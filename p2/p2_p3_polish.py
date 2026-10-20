#!/usr/bin/env python3
"""P2.3 active-set KKT polish (rank-checked, memory-bounded).

1. Start from the OSQP solution and its active Bernstein rows.
2. Reduce the active rows to a maximal independent subset via pivoted QR of
   the Gram matrix C_A C_A^T (the saddle-point matrix need not be invertible
   if the rows are dependent).
3. Solve the equality-constrained KKT system by the Schur complement, with one
   cached factorization of K_FF and batched solves (no n_free x r dense block).
4. Iterate add/drop on violated primal rows / negative multipliers.
5. Audit stationarity, primal, dual, complementarity, oxide equilibrium and
   the interface residual AFTER activation.
"""
from __future__ import annotations
import json
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
import scipy.linalg as sla

import p2_audit as A
import p2_interface as PI
from p2_bernstein import B

M = A.M
V_D = A.V_DRAIN


def build():
    z = np.load("p2_r3_fixed.npz")
    nodes, tris, mat = z["nodes"], z["tris"], z["mat"]
    x = nodes[:, 0]
    K, _ = A.assemble(nodes, tris, mat)
    area = np.abs(A.element_geom(nodes, tris)[2])

    cxs = nodes[tris][:, :3, 0].mean(axis=1)
    si = mat == M.SI
    body = si & (cxs <= A.BODY_BOTTOM + 1e-12)
    charged = si & ~body

    f = np.zeros(len(nodes))
    for e in np.flatnonzero(charged):
        f[tris[e][3:]] += A.RHO*area[e]/3.0

    fixed = np.zeros(len(nodes), bool)
    val = np.zeros(len(nodes))
    si_nodes = np.unique(tris[si])
    body_nodes = np.unique(tris[body])
    gate = np.unique(tris[mat == M.GATE])
    shield = np.unique(tris[mat == M.SHIELD])
    src = si_nodes[x[si_nodes] <= 1e-12]
    drn = np.flatnonzero(x >= M.BOTTOM_DEPTH - 1e-9)
    for s in (body_nodes, gate, shield, src):
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
    return dict(nodes=nodes, tris=tris, mat=mat, area=area, si=si, body=body,
                charged=charged, fixed=fixed, val=val, free=free, idx=idx,
                cons=cons, C=C, C_F=C_F, cD=cD, K=K, Kff=Kff, g=g, f=f)


def independent(CW, tol=1e-10):
    G = (CW @ CW.T).toarray()
    _, R, piv = sla.qr(G, pivoting=True)
    d = np.abs(np.diag(R))
    r = int((d > tol*d.max()).sum())
    return np.array(piv[:r])


def schur(Ksol, Cr, cD_W, kg, batch=64):
    r = Cr.shape[0]
    S = np.empty((r, r))
    for s in range(0, r, batch):
        e = min(s+batch, r)
        Xb = Ksol(Cr[s:e].T.toarray())
        S[s:e, :] = (Cr @ Xb).T
    return S


def finish(Ksol, Cr, kg, lam, batch=64):
    u = kg.copy()
    r = Cr.shape[0]
    for s in range(0, r, batch):
        e = min(s+batch, r)
        Xb = Ksol(Cr[s:e].T.toarray())
        u += Xb @ lam[s:e]
    return u


def main():
    d = build()
    z0 = np.load("p2_p3_r3.npz")
    u0 = z0["u"][d["idx"]]
    C_F, cD, Kff, g = d["C_F"], d["cD"], d["Kff"], d["g"]

    Ksol = spla.factorized(Kff)
    kg = Ksol(g)

    W = np.flatnonzero(C_F @ u0 + cD < 1e-7)
    print(f"[polish] initial active={len(W)}")
    W = W[independent(C_F[W])]
    print(f"[polish] independent active={len(W)}")

    for it in range(8):
        Cr = C_F[W]
        S = schur(Ksol, Cr, cD[W], kg)
        rhs = -cD[W] - Cr @ kg
        lam = sla.solve(S, rhs, assume_a="sym")
        u = finish(Ksol, Cr, kg, lam)
        slack = C_F @ u + cD
        print(f"  iter {it}: |W|={len(W)} min slack={slack.min():+.3e} "
              f"lam min={lam.min():+.3e}")
        if slack.min() < -1e-9:
            j = int(np.argmin(slack))
            if j in W:
                break
            W = np.union1d(W, [j])
            W = W[independent(C_F[W])]
            continue
        if lam.min() < -1e-9:
            j = W[int(np.argmin(lam))]
            W = np.setdiff1d(W, [j])
            W = W[independent(C_F[W])]
            continue
        break
    else:
        print("  WARNING: polish did not converge")

    uf = u
    u_full = d["val"].copy()
    u_full[d["idx"]] = uf
    phi = V_D - u_full
    lam_full = np.zeros(len(cD)); lam_full[W] = lam
    r = Kff @ uf - g
    stat = r - C_F.T @ lam_full
    slack = C_F @ uf + cD
    comp = np.max(np.abs(lam_full*slack))

    print(f"\n[polish] final |W|={len(W)}")
    print(f"  stationarity rel = {np.linalg.norm(stat)/max(np.linalg.norm(r),1e-30):.3e}")
    print(f"  min slack = {slack.min():+.3e}   lam min = {lam_full.min():+.3e}   "
          f"lam max = {lam_full.max():.3e}")
    print(f"  complementarity = {comp:.3e}")
    print(f"  phi max = {phi.max():.9f}   overshoot = {phi.max()-V_D:+.3e}")

    # oxide equilibrium on pure-oxide free rows (positions in the free vector)
    ox_nodes = [i for i in np.unique(d["tris"][d["mat"] == M.OX])
                if int(i) not in set(np.unique(d["tris"][d["si"]]).tolist())]
    pos = np.searchsorted(d["idx"], np.array(ox_nodes, dtype=int))
    pos = pos[pos < len(d["idx"])]
    print(f"  pure-oxide free rows max|K u| = {np.abs(r[pos]).max():.3e}")

    recs = PI.interface_flux(d["nodes"], d["tris"], d["mat"], phi)
    w = np.array([q["ds"] for q in recs])
    ds = np.array([q["d_si"] for q in recs])
    do = np.array([q["d_ox"] for q in recs])
    Rab = float(np.sum(w*np.abs(ds-do))/np.sum(w*(np.abs(ds)+np.abs(do))/2))
    print(f"  interface R_abs after activation = {Rab:.6f} (n={len(recs)})")

    np.savez_compressed("p2_p3_r3_polished.npz", nodes=d["nodes"], tris=d["tris"],
                        mat=d["mat"], phi=phi, u=u_full, lam=lam_full, active=W)
    json.dump(dict(active=int(len(W)),
                   stationarity=float(np.linalg.norm(stat)/max(np.linalg.norm(r), 1e-30)),
                   min_slack=float(slack.min()), lam_min=float(lam_full.min()),
                   lam_max=float(lam_full.max()), complementarity=float(comp),
                   overshoot=float(phi.max()-V_D), R_abs_after=Rab,
                   oxide_res=float(np.abs(r[pos]).max())),
              open("p2_p3_polish.json", "w"), indent=2)
    print("wrote p2_p3_r3_polished.npz / p2_p3_polish.json")


if __name__ == "__main__":
    main()
