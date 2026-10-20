#!/usr/bin/env python3
"""P2.3 obstacle solve on r4, P1-matched model, with KKT polish.

Same construction as the r3 solve, using the corrected operator and the
vectorised assembly. One assembly, one OSQP solve, rank-reduced Schur polish.
"""
from __future__ import annotations
import json
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
import scipy.linalg as sla
import osqp

import p2_audit as A
from p2_bernstein import B
from p2_nested_consistency import assemble_p2
from p2_p3_polish import independent, schur, finish

M = A.M
V_D = A.V_DRAIN


def main():
    z = np.load("p2_r4_fixed.npz")
    nodes, tris, mat = z["nodes"], z["tris"], z["mat"]
    x = nodes[:, 0]
    cxs = nodes[tris][:, :3, 0].mean(axis=1)
    si = mat == M.SI
    body = si & (cxs <= A.BODY_BOTTOM + 1e-12)
    charged = si & ~body

    K, f = assemble_p2(nodes, tris, mat, charged)

    fixed = np.zeros(len(nodes), bool)
    val = np.zeros(len(nodes))
    si_nodes = np.unique(tris[si])
    gate = np.unique(tris[mat == M.GATE])
    shield = np.unique(tris[mat == M.SHIELD])
    src = si_nodes[x[si_nodes] <= 1e-12]
    drn = np.flatnonzero(x >= M.BOTTOM_DEPTH - 1e-9)
    for s in (np.unique(tris[body]), gate, shield, src):
        fixed[s] = True
    val[fixed] = V_D
    fixed[drn] = True
    val[drn] = 0.0
    free = ~fixed
    idx = np.flatnonzero(free)

    cons = np.flatnonzero(charged)
    nC = 6*len(cons)
    C = sp.coo_matrix((np.tile(B.ravel(), len(cons)),
                       (np.repeat(np.arange(nC), 6), np.tile(tris[cons], (1, 6)).ravel())),
                      shape=(nC, len(nodes))).tocsr()
    C_F = C[:, idx].tocsc()
    cD = C[:, fixed] @ val[fixed]
    Kff = K[idx][:, idx].tocsc()
    g = -f[idx] - (K[idx][:, fixed] @ val[fixed])
    sK = float(Kff.diagonal().max())
    Kff_s = (Kff/sK).tocsc()
    g_s = g/sK

    print(f"[r4] dofs={len(nodes)} free={len(idx)} constraints={nC} "
          f"charged={len(cons)} body={int(body.sum())}", flush=True)

    prob = osqp.OSQP()
    prob.setup(P=Kff_s, q=-g_s, A=C_F, l=-cD, u=np.full(len(cD), np.inf),
               verbose=False, eps_abs=1e-8, eps_rel=1e-8, max_iter=200000,
               polishing=True, scaling=10)
    res = prob.solve()
    print(f"[r4] OSQP {res.info.status} iters={res.info.iter} "
          f"obj={res.info.obj_val:.6e}", flush=True)
    u0 = res.x

    # polish
    Ksol = spla.factorized(Kff)
    kg = Ksol(g)
    W = np.flatnonzero(C_F @ u0 + cD < 1e-7)
    W = W[independent(C_F[W])]
    print(f"[r4] independent active={len(W)}", flush=True)
    for it in range(8):
        Cr = C_F[W]
        S = schur(Ksol, Cr, cD[W], kg)
        rhs = -cD[W] - Cr @ kg
        lam = sla.solve(S, rhs, assume_a="sym")
        u = finish(Ksol, Cr, kg, lam)
        slack = C_F @ u + cD
        print(f"  iter {it}: |W|={len(W)} min slack={slack.min():+.3e} "
              f"lam min={lam.min():+.3e}", flush=True)
        if slack.min() < -1e-9:
            j = int(np.argmin(slack))
            if j in W:
                break
            W = np.union1d(W, [j]); W = W[independent(C_F[W])]; continue
        if lam.min() < -1e-9:
            j = W[int(np.argmin(lam))]
            W = np.setdiff1d(W, [j]); W = W[independent(C_F[W])]; continue
        break

    uf = u
    u_full = val.copy(); u_full[idx] = uf
    phi = V_D - u_full
    lam_full = np.zeros(len(cD)); lam_full[W] = lam
    r = Kff @ uf - g
    stat = r - C_F.T @ lam_full
    slack = C_F @ uf + cD
    comp = np.max(np.abs(lam_full*slack))
    area = np.abs(A.element_geom(nodes, tris)[2])
    ue = u_full[tris]
    fully = charged & (np.abs(ue) < 1e-8).all(axis=1)
    print(f"[r4] stationarity={np.linalg.norm(stat)/max(np.linalg.norm(r),1e-30):.3e} "
          f"min slack={slack.min():+.3e} lam min={lam_full.min():+.3e} "
          f"comp={comp:.3e}", flush=True)
    print(f"[r4] phi max={phi.max():.9f} overshoot={phi.max()-V_D:+.3e} "
          f"fully neutral area={area[fully].sum():.9f} elements={int(fully.sum())}",
          flush=True)

    np.savez_compressed("p2_p3_r4_polished.npz", nodes=nodes, tris=tris, mat=mat,
                        phi=phi, u=u_full, lam=lam_full, active=W)
    json.dump(dict(status=str(res.info.status), iters=int(res.info.iter),
                   free=int(len(idx)), constraints=int(nC), active=int(len(W)),
                   stationarity=float(np.linalg.norm(stat)/max(np.linalg.norm(r), 1e-30)),
                   min_slack=float(slack.min()), lam_min=float(lam_full.min()),
                   lam_max=float(lam_full.max()), complementarity=float(comp),
                   overshoot=float(phi.max()-V_D),
                   neutral_area=float(area[fully].sum())),
              open("p2_p3_r4.json", "w"), indent=2)
    print("wrote p2_p3_r4_polished.npz / p2_p3_r4.json")


if __name__ == "__main__":
    main()
