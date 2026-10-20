#!/usr/bin/env python3
"""L1 scaled/preconditioned KKT polish (Option a).

Symmetric KKT (F7 multiplier convention):
    [ K_FF   -C_A' ] [ u   ]   [ g     ]
    [ -C_A     0   ] [ lam ] = [ cD_A  ]
with lam >= 0.  Diagonal congruence scaling D H D is symmetric, so MINRES
applies.  Acceptance is computed on the ORIGINAL unscaled problem.
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
    return dict(nodes=nodes, tris=tris, mat=mat, si=si, fixed=fixed, val=val,
                idx=idx, C_F=C_F, cD=cD, K=K, Kff=Kff, g=g, f=f, nfree=len(idx))


def main():
    d = build()
    pol = np.load("p2_L1_polished.npz")
    W = pol["active"].astype(int)
    u = pol["u"][d["idx"]].copy()
    lam_full = pol["lam"].copy()
    lam = lam_full[W].copy()
    Kff, C_F, cD, g = d["Kff"], d["C_F"], d["cD"], d["g"]
    CW = C_F[W]
    n = d["nfree"]; m = len(W)

    # diagonal congruence scaling
    dK = np.maximum(Kff.diagonal(), 1e-12)
    du = 1.0/np.sqrt(dK)
    CWu = CW.multiply(du[None, :])
    rn = np.sqrt(np.asarray(CWu.multiply(CWu).sum(axis=1)).ravel())
    rn = np.maximum(rn, 1e-12)
    dl = 1.0/rn
    D = np.concatenate([du, dl])

    def Hmv(zt):
        z = D*zt
        uu, ll = z[:n], z[n:]
        return D*np.concatenate([Kff@uu - CW.T@ll, -CW@uu])

    Hop = spla.LinearOperator((n+m, n+m), matvec=Hmv, dtype=float)
    b = np.concatenate([g, cD[W]])
    bt = D*b
    zt0 = np.concatenate([u/du, lam/dl])
    print(f"[scaled] |W|={m}  du range {du.min():.3e}..{du.max():.3e}  "
          f"dl range {dl.min():.3e}..{dl.max():.3e}", flush=True)

    zt, info = spla.minres(Hop, bt, x0=zt0, rtol=1e-14, maxiter=200000)
    z = D*zt
    u = z[:n].copy(); lam = z[n:].copy()
    print(f"[scaled] minres info={info}", flush=True)

    def unscaled():
        r = Kff@u - g
        stat = r - CW.T@lam
        slack = C_F@u + cD
        return r, stat, slack

    for k in range(3):
        r, stat, slack = unscaled()
        rel = float(np.linalg.norm(stat)/max(np.linalg.norm(r), 1e-30))
        print(f"  refine {k}: stat_rel={rel:.3e} min slack={slack.min():+.3e}", flush=True)
        if rel < 5e-10:
            break
        res = np.concatenate([-stat, -slack[W]])
        rt = D*res
        zt_corr, info = spla.minres(Hop, rt, rtol=1e-14, maxiter=200000)
        zc = D*zt_corr
        u += zc[:n]; lam += zc[n:]

    r, stat, slack = unscaled()
    rel = float(np.linalg.norm(stat)/max(np.linalg.norm(r), 1e-30))
    lam_full = np.zeros(len(cD)); lam_full[W] = lam
    comp = float(np.max(np.abs(lam_full*slack)))
    u_full = d["val"].copy(); u_full[d["idx"]] = u
    overshoot = float(-u_full.min())
    si_set = set(np.unique(d["tris"][d["si"]]).tolist())
    ox = [int(i) for i in np.unique(d["tris"][d["mat"] == M.OX]) if int(i) not in si_set]
    ox_free = np.intersect1d(np.array(ox, dtype=int), d["idx"])
    Ku = d["K"] @ u_full + d["f"]
    ox_res = float(np.abs(Ku[ox_free]).max()) if len(ox_free) else 0.0
    viol = int((slack < -1e-9).sum())
    print(f"[scaled] final stat_rel={rel:.3e} comp={comp:.3e} min slack={slack.min():+.3e} "
          f"viol={viol} lam_min={lam_full.min():+.3e} oxide={ox_res:.3e} "
          f"overshoot={overshoot:+.3e}")
    ok = (rel < 1e-8 and comp < 1e-10 and slack.min() >= -1e-9 and
          lam_full.min() >= -1e-9 and ox_res < 1e-9 and overshoot < 1e-9)
    print(f"  GATE: {'PASS' if ok else 'FAIL'}")
    np.savez_compressed("p2_L1_scaled.npz", nodes=d["nodes"], tris=d["tris"],
                        mat=d["mat"], u=u_full, lam=lam_full, active=W)
    json.dump(dict(active=int(m), stationarity=rel, complementarity=comp,
                   min_slack=float(slack.min()), violations=viol,
                   lam_min=float(lam_full.min()), oxide=ox_res,
                   overshoot=overshoot, gate=bool(ok),
                   minres_info=int(info)), open("p2_L1_scaled.json", "w"), indent=2)
    print("wrote p2_L1_scaled.npz / .json")


if __name__ == "__main__":
    main()
