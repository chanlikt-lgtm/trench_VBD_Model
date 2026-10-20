#!/usr/bin/env python3
"""Step C: matrix-free symmetric saddle-point MINRES polish of the L1 QP.

Symmetric KKT form (no factorization of K_FF):
    [ K_FF  C_A' ] [ u   ]   [ g     ]
    [ C_A    0   ] [ lam'] = [ -cD_A ]
with lam = -lam' >= 0.
"""
from __future__ import annotations
import json
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

import p2_audit as A
from p2_bernstein import B
from p2_nested_consistency import assemble_p2

M = A.M
V_D = A.V_DRAIN


def build():
    z = np.load("p2_L1_p2mesh.npz")
    nodes, tris, mat = z["nodes"], z["tris"], z["mat"]
    used = np.unique(tris.ravel())
    o2n = np.full(len(nodes), -1, dtype=int)
    o2n[used] = np.arange(len(used))
    nodes = nodes[used]; tris = o2n[tris]
    x = nodes[:, 0]
    cxs = nodes[tris][:, :3, 0].mean(axis=1)
    si = mat == M.SI
    body = si & (cxs <= A.BODY_BOTTOM + 1e-12)
    charged = si & ~body
    K, f = assemble_p2(nodes, tris, mat, charged)
    fixed = np.zeros(len(nodes), bool); val = np.zeros(len(nodes))
    si_nodes = np.unique(tris[si])
    gate = np.unique(tris[mat == M.GATE]); shield = np.unique(tris[mat == M.SHIELD])
    src = si_nodes[x[si_nodes] <= 1e-12]; drn = np.flatnonzero(x >= M.BOTTOM_DEPTH - 1e-9)
    for s in (np.unique(tris[body]), gate, shield, src):
        fixed[s] = True
    val[fixed] = V_D
    fixed[drn] = True; val[drn] = 0.0
    idx = np.flatnonzero(~fixed)
    cons = np.flatnonzero(charged); nC = 6*len(cons)
    C = sp.coo_matrix((np.tile(B.ravel(), len(cons)),
                       (np.repeat(np.arange(nC), 6), np.tile(tris[cons], (1, 6)).ravel())),
                      shape=(nC, len(nodes))).tocsr()
    C_F = C[:, idx].tocsc()
    cD = C[:, fixed] @ val[fixed]
    Kff = K[idx][:, idx].tocsc()
    g = -f[idx] - (K[idx][:, fixed] @ val[fixed])
    chk = np.load("p2_L1_osqp_checkpoint.npz")
    xf = chk["x"]; y = chk["y"]; sK = float(chk["sK"])
    lam_osqp = -sK*y
    return dict(nodes=nodes, tris=tris, mat=mat, si=si, fixed=fixed, val=val,
                idx=idx, C=C, C_F=C_F, cD=cD, K=K, Kff=Kff, g=g, f=f,
                xf=xf, lam=lam_osqp, nfree=len(idx))


def kkt_solve(d, W, x0u, x0l, maxiter=50000, rtol=1e-13):
    Kff, C_F = d["Kff"], d["C_F"]
    CW = C_F[W]
    n = d["nfree"]; m = len(W)
    def mv(v):
        u = v[:n]; l = v[n:]
        return np.concatenate([Kff @ u + CW.T @ l, CW @ u])
    Aop = spla.LinearOperator((n+m, n+m), matvec=mv, dtype=float)
    dK = Kff.diagonal(); dC = np.asarray(CW.multiply(CW).sum(axis=1)).ravel()
    dK = np.where(dK > 0, dK, 1.0); dC = np.where(dC > 0, dC, 1.0)
    def pmv(v):
        return np.concatenate([v[:n]/dK, v[n:]/dC])
    Mop = spla.LinearOperator((n+m, n+m), matvec=pmv, dtype=float)
    b = np.concatenate([d["g"], -d["cD"][W]])
    x0 = np.concatenate([x0u, x0l])
    sol, info = spla.minres(Aop, b, x0=x0, M=Mop, rtol=rtol, maxiter=maxiter)
    res = mv(sol) - b
    return sol[:n], sol[n:], int(info), float(np.linalg.norm(res)/max(np.linalg.norm(b), 1e-30))


def main():
    d = build()
    lam_vec = d["lam"].copy()
    u = d["xf"].copy()
    for it in range(8):
        W = np.flatnonzero(lam_vec > 1e-8)
        u, l, info, kkt_rel = kkt_solve(d, W, u, -lam_vec[W])
        lam_new = -l
        slack = d["C_F"] @ u + d["cD"]
        lam_vec = np.zeros(len(d["cD"]))
        lam_vec[W] = lam_new
        viol = np.flatnonzero(slack < -1e-9)
        neg = np.flatnonzero(lam_vec < -1e-9)
        print(f"  iter {it}: |W|={len(W)} info={info} kkt_rel={kkt_rel:.3e} "
              f"min slack={slack.min():+.3e} viol={len(viol)} neg={len(neg)}",
              flush=True)
        if viol.size:
            lam_vec[viol] = np.maximum(lam_vec[viol], 1e-6)
        if neg.size:
            lam_vec[neg] = 0.0
        if not viol.size and not neg.size:
            break
    W = np.flatnonzero(lam_vec > 0.0)
    l = -lam_vec[W]

    lam_full = np.zeros(len(d["cD"])); lam_full[W] = -l
    r = d["Kff"] @ u - d["g"]
    stat = r - d["C_F"].T @ lam_full
    slack = d["C_F"] @ u + d["cD"]
    comp = float(np.max(np.abs(lam_full*slack)))
    rel = float(np.linalg.norm(stat)/max(np.linalg.norm(r), 1e-30))
    u_full = d["val"].copy(); u_full[d["idx"]] = u
    overshoot = float(-u_full.min())
    si_set = set(np.unique(d["tris"][d["si"]]).tolist())
    ox = [int(i) for i in np.unique(d["tris"][d["mat"] == M.OX]) if int(i) not in si_set]
    ox_free = np.intersect1d(np.array(ox, dtype=int), d["idx"])
    Ku = d["K"] @ u_full + d["f"]
    ox_res = float(np.abs(Ku[ox_free]).max()) if len(ox_free) else 0.0
    print(f"[C] final |W|={len(W)}  stationarity={rel:.3e}  comp={comp:.3e}")
    print(f"    min slack={slack.min():+.3e}  lam min={lam_full.min():+.3e}  "
          f"oxide={ox_res:.3e}  overshoot={overshoot:+.3e}")
    ok = (slack.min() >= -1e-9 and lam_full.min() >= -1e-9 and rel < 1e-8
          and comp < 1e-10 and ox_res < 1e-9 and overshoot < 1e-9)
    print(f"    GATE: {'PASS' if ok else 'FAIL'}")
    np.savez_compressed("p2_L1_polished.npz", nodes=d["nodes"], tris=d["tris"],
                        mat=d["mat"], u=u_full, lam=lam_full, active=W)
    json.dump(dict(active=int(len(W)), stationarity=rel, complementarity=comp,
                   min_slack=float(slack.min()), lam_min=float(lam_full.min()),
                   oxide=ox_res, overshoot=overshoot, gate=bool(ok)),
              open("p2_L1_polished.json", "w"), indent=2)
    print("wrote p2_L1_polished.npz / .json")


if __name__ == "__main__":
    main()
