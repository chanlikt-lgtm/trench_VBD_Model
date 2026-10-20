#!/usr/bin/env python3
"""L1b Step B/C: scaled matrix-free MINRES KKT polish from the OSQP checkpoint.

Same convention as p2_L1_scaled.py (symmetric KKT, lam >= 0, acceptance
computed on the ORIGINAL unscaled system), plus bulk active-set updates so a
bad OSQP active set cannot poison the polish.
"""
from __future__ import annotations
import json
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

import p2_audit as A
from p2_bernstein import B
from p2_nested_consistency import assemble_p2
from p2_p3_polish import independent

M = A.M
V_D = A.V_DRAIN


def build():
    z = np.load("p2_L1b_p2mesh.npz")
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


def polish(d, W, u, lam, tag=""):
    """Run scaled MINRES + residual correction on the given active set."""
    Kff, C_F, cD, g = d["Kff"], d["C_F"], d["cD"], d["g"]
    n = d["nfree"]; m = len(W)
    CW = C_F[W]

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
    zt0 = np.concatenate([u/du, lam/dl])
    zt, info = spla.minres(Hop, D*b, x0=zt0, rtol=1e-14, maxiter=200000)
    z = D*zt
    u = z[:n].copy(); lam = z[n:].copy()

    for k in range(6):
        r = Kff@u - g
        stat = r - CW.T@lam
        slackW = CW@u + cD[W]
        rel = float(np.linalg.norm(stat)/max(np.linalg.norm(r), 1e-30))
        print(f"  {tag} refine {k}: stat_rel={rel:.3e} min slackW={slackW.min():+.3e}",
              flush=True)
        if rel < 5e-10:
            break
        res = np.concatenate([-stat, -slackW])
        zt_corr, _ = spla.minres(Hop, D*res, rtol=1e-14, maxiter=200000)
        zc = D*zt_corr
        u += zc[:n]; lam += zc[n:]
    return u, lam, int(info)


def gate(d, W, u, lam):
    Kff, C_F, cD, g = d["Kff"], d["C_F"], d["cD"], d["g"]
    CW = C_F[W]
    r = Kff@u - g
    stat = r - CW.T@lam
    rel = float(np.linalg.norm(stat)/max(np.linalg.norm(r), 1e-30))
    slack = C_F@u + cD
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
    return dict(rel=rel, slack=slack, lam_full=lam_full, comp=comp,
                ox=ox_res, viol=viol, overshoot=overshoot, u_full=u_full,
                lam_min=float(lam_full.min()))


def main():
    d = build()
    chk = np.load("p2_L1b_osqp_checkpoint.npz")
    u = chk["x"].copy()
    print(f"[L1b-B] checkpoint status={chk['status']} iters={int(chk['iters'])} "
          f"free={d['nfree']}", flush=True)

    Kff, C_F, cD = d["Kff"], d["C_F"], d["cD"]
    slack0 = C_F@u + cD
    W = np.flatnonzero(slack0 < 1e-7)
    W = W[independent(C_F[W])]
    lam = np.zeros(len(W))
    print(f"[L1b-B] initial |W|={len(W)}", flush=True)

    best = None
    for outer in range(6):
        u, lam, info = polish(d, W, u, lam, tag=f"[{outer}]")
        gres = gate(d, W, u, lam)
        print(f"  outer {outer}: |W|={len(W)} stat={gres['rel']:.3e} "
              f"min slack={gres['slack'].min():+.3e} lam_min={gres['lam_min']:+.3e} "
              f"viol={gres['viol']}", flush=True)
        ok = (gres["rel"] < 1e-8 and gres["comp"] < 1e-10 and
              gres["slack"].min() >= -1e-9 and gres["lam_min"] >= -1e-9 and
              gres["ox"] < 1e-9 and gres["overshoot"] < 1e-9)
        if best is None or gres["rel"] < best[0]["rel"]:
            best = (gres, W.copy(), u.copy(), lam.copy())
        if ok:
            print("  GATE: PASS", flush=True)
            break
        # bulk active-set update
        add = np.flatnonzero(gres["slack"] < -1e-12)
        drop = W[gres["lam_full"][W] < -1e-9]
        if len(add) == 0 and len(drop) == 0:
            print("  active set stable, no further bulk update", flush=True)
            break
        if len(add):
            W = np.union1d(W, add)
        if len(drop):
            W = np.setdiff1d(W, drop)
        W = W[independent(C_F[W])]
        lam = np.zeros(len(W))

    gres, W, u, lam = best
    W = np.asarray(W, dtype=int)
    ok = (gres["rel"] < 1e-8 and gres["comp"] < 1e-10 and
          gres["slack"].min() >= -1e-9 and gres["lam_min"] >= -1e-9 and
          gres["ox"] < 1e-9 and gres["overshoot"] < 1e-9)
    print(f"[L1b] final |W|={len(W)} stat={gres['rel']:.3e} comp={gres['comp']:.3e} "
          f"min slack={gres['slack'].min():+.3e} viol={gres['viol']} "
          f"lam_min={gres['lam_min']:+.3e} oxide={gres['ox']:.3e} "
          f"overshoot={gres['overshoot']:+.3e}", flush=True)
    print(f"  GATE: {'PASS' if ok else 'FAIL'}", flush=True)

    np.savez_compressed("p2_L1b_scaled.npz", nodes=d["nodes"], tris=d["tris"],
                        mat=d["mat"], u=gres["u_full"], lam=gres["lam_full"],
                        active=W)
    json.dump(dict(active=int(len(W)), stationarity=gres["rel"],
                   complementarity=gres["comp"],
                   min_slack=float(gres["slack"].min()),
                   violations=int(gres["viol"]),
                   lam_min=float(gres["lam_min"]), oxide=gres["ox"],
                   overshoot=gres["overshoot"], gate=bool(ok),
                   minres_info=int(info)), open("p2_L1b_scaled.json", "w"), indent=2)
    print("wrote p2_L1b_scaled.npz / .json")


if __name__ == "__main__":
    main()
