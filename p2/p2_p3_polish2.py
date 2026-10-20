#!/usr/bin/env python3
"""P2.3 KKT polish — vectorised assembly version.

Uses assemble_p2 from p2_nested_consistency (vectorised, corrected QW) instead
of the serial loop, then the same rank-reduced Schur-complement solve.
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
from p2_nested_consistency import assemble_p2
from p2_p3_polish import independent, schur, finish

M = A.M
V_D = A.V_DRAIN


def build():
    z = np.load("p2_r3_fixed.npz")
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
    return dict(nodes=nodes, tris=tris, mat=mat, si=si, body=body, fixed=fixed,
                val=val, idx=idx, cons=cons, C=C, C_F=C_F, cD=cD, K=K, Kff=Kff, g=g)


def main():
    d = build()
    z0 = np.load("p2_p3_r3.npz")
    u0 = z0["u"][d["idx"]]
    C_F, cD, Kff, g = d["C_F"], d["cD"], d["Kff"], d["g"]

    Ksol = spla.factorized(Kff)
    kg = Ksol(g)

    W = np.flatnonzero(C_F @ u0 + cD < 1e-7)
    W = W[independent(C_F[W])]
    print(f"[polish] initial independent active={len(W)}", flush=True)

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
                print("    most violated row already active; stopping")
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
    print(f"  min slack = {slack.min():+.3e}  lam min = {lam_full.min():+.3e}  "
          f"lam max = {lam_full.max():.3e}")
    print(f"  complementarity = {comp:.3e}")
    print(f"  phi max = {phi.max():.9f}  overshoot = {phi.max()-V_D:+.3e}")

    # save the polished solution before any post-processing so a failure in the
    # audit cannot lose it
    np.savez_compressed("p2_p3_r3_polished.npz", nodes=d["nodes"], tris=d["tris"],
                        mat=d["mat"], phi=phi, u=u_full, lam=lam_full, active=W)
    print("  saved p2_p3_r3_polished.npz", flush=True)

    si_set = set(np.unique(d["tris"][d["si"]]).tolist())
    ox_nodes = [int(i) for i in np.unique(d["tris"][d["mat"] == M.OX])
                if int(i) not in si_set]
    pos = np.searchsorted(d["idx"], np.array(ox_nodes, dtype=int))
    pos = pos[pos < len(d["idx"])]
    print(f"  pure-oxide free rows max|K u| = {np.abs(r[pos]).max():.3e}", flush=True)

    recs = PI.interface_flux(d["nodes"], d["tris"], d["mat"], phi)
    w = np.array([q["ds"] for q in recs])
    ds = np.array([q["d_si"] for q in recs])
    do = np.array([q["d_ox"] for q in recs])
    Rab = float(np.sum(w*np.abs(ds-do))/np.sum(w*(np.abs(ds)+np.abs(do))/2))
    print(f"  interface R_abs after activation = {Rab:.6f} (n={len(recs)})", flush=True)

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
