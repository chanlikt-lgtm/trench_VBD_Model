#!/usr/bin/env python3
"""P1 obstacle solution on the P2 mesh's vertex triangles.

Same geometry, materials, P1-matched body (uncharged, all its nodes pinned at
u_D = V_D), same Dirichlet sets and u = V_D - phi convention as the P2 solve.
This gives a P1 obstacle field on identical geometry for a common zero-set
comparison.
"""
from __future__ import annotations
import json
import sys
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

import p2_audit as A

M = A.M
V_D = A.V_DRAIN
EPS = {M.SI: A.EPS_SI, M.OX: A.EPS_OX, M.GATE: A.EPS_SI, M.SHIELD: A.EPS_SI}


def main():
    FN = sys.argv[1] if len(sys.argv) > 1 else "p2_r3_fixed.npz"
    TAG = sys.argv[2] if len(sys.argv) > 2 else ""
    z = np.load(FN)
    nodes6, tris6, mat = z["nodes"], z["tris"], z["mat"]

    tv = tris6[:, :3]
    used = np.unique(tv.ravel())
    o2n = np.full(len(nodes6), -1, dtype=int)
    o2n[used] = np.arange(len(used))
    nodes = nodes6[used]
    tris = o2n[tv]
    x = nodes[:, 0]

    cxs = nodes6[tris6][:, :3, 0].mean(axis=1)
    si = mat == M.SI
    body = si & (cxs <= A.BODY_BOTTOM + 1e-12)
    charged = si & ~body

    p = nodes[tris]
    e1 = p[:, 1] - p[:, 0]; e2 = p[:, 2] - p[:, 0]
    det = e1[:, 0]*e2[:, 1] - e2[:, 0]*e1[:, 1]
    area = 0.5*np.abs(det)
    b = np.stack([p[:, 1, 1]-p[:, 2, 1], p[:, 2, 1]-p[:, 0, 1], p[:, 0, 1]-p[:, 1, 1]], 1)
    c = np.stack([p[:, 2, 0]-p[:, 1, 0], p[:, 0, 0]-p[:, 2, 0], p[:, 1, 0]-p[:, 0, 0]], 1)
    eps = np.array([EPS[m] for m in mat])
    rows = np.zeros((len(tris), 9), dtype=int)
    cols = np.zeros((len(tris), 9), dtype=int)
    vals = np.zeros((len(tris), 9))
    for e in range(len(tris)):
        ke = eps[e]*(np.outer(b[e], b[e]) + np.outer(c[e], c[e]))/(4*area[e])
        ids = tris[e]
        rows[e] = np.repeat(ids, 3)
        cols[e] = np.tile(ids, 3)
        vals[e] = ke.ravel()
    K = sp.coo_matrix((vals.ravel(), (rows.ravel(), cols.ravel())),
                      shape=(len(nodes), len(nodes))).tocsr()

    f = np.zeros(len(nodes))
    for e in np.flatnonzero(charged):
        f[tris[e]] += A.RHO*area[e]/3.0

    fixed = np.zeros(len(nodes), bool)
    val = np.zeros(len(nodes))
    grounded = np.zeros(len(nodes), bool)
    for e in np.flatnonzero(body | (mat == M.GATE) | (mat == M.SHIELD)):
        grounded[tris[e]] = True
    si_nodes = np.unique(tris[charged])
    src = si_nodes[x[si_nodes] <= 1e-12]
    drn = np.flatnonzero(x >= M.BOTTOM_DEPTH - 1e-9)
    fixed[grounded] = True; val[grounded] = V_D
    fixed[src] = True; val[src] = V_D
    fixed[drn] = True; val[drn] = 0.0

    free = ~fixed
    rhs = -f - K @ val
    idx = np.flatnonzero(free)
    u = val.copy()
    u[idx] = spla.spsolve(K[idx][:, idx].tocsc(), rhs[idx])
    active = np.zeros(len(nodes), bool)
    active[idx] = (u[idx] < 0) & np.isin(idx, si_nodes)

    for it in range(50):
        free2 = free & ~active
        u2 = val.copy()
        u2[active] = 0.0
        j = np.flatnonzero(free2)
        if j.size:
            u2[j] = spla.spsolve(K[j][:, j].tocsc(), rhs[j])
        # physical obstacle residual for K u = -f with Dirichlet elimination;
        # u2 already contains the prescribed values, so do NOT subtract rhs
        r = K @ u2 + f
        r[fixed] = 0.0
        dg = K.diagonal(); dg = np.where(dg > 0, dg, 1.0)
        new = np.zeros(len(nodes), bool)
        cand = si_nodes[free[si_nodes]]
        new[cand] = (u2[cand] - r[cand]/dg[cand]) <= 0
        if np.array_equal(new, active):
            break
        active = new
    else:
        print("WARNING: P1 complementarity did not converge")

    phi = V_D - u2
    cand = si_nodes[free[si_nodes]]
    res = (K @ u2 + f)[cand]
    act = active[cand]
    print(f"[P1obs] KKT: min u={u2[cand].min():+.3e}  "
          f"min residual(active)={res[act].min() if act.any() else 0.0:+.3e}  "
          f"max |residual|(inactive)={np.abs(res[~act]).max() if (~act).any() else 0.0:.3e}  "
          f"comp={np.abs(u2[cand][act]*res[act]).max() if act.any() else 0.0:.3e}")
    act_el = np.zeros(len(tris), bool)
    act_el[charged] = active[tris[charged]].any(axis=1)
    print(f"[P1obs] nodes={len(nodes)} tris={len(tris)} charged={int(charged.sum())}")
    print(f"[P1obs] active nodes={int(active.sum())}  iterations={it+1}")
    print(f"[P1obs] active area={area[act_el].sum():.9f} um^2  "
          f"charged area={area[charged].sum():.9f}")
    print(f"[P1obs] phi max={phi.max():.9f}  min u={u2.min():.3e}")

    np.savez_compressed(f"p1_on_p2mesh_obstacle{TAG}.npz", nodes=nodes, tris=tris,
                        mat=mat, u=u2, phi=phi, active=active)
    json.dump(dict(active_nodes=int(active.sum()),
                   active_area=float(area[act_el].sum()),
                   charged_area=float(area[charged].sum()),
                   iterations=int(it+1), phi_max=float(phi.max()),
                   min_u=float(u2.min())),
              open(f"p1_on_p2mesh_obstacle{TAG}.json", "w"), indent=2)
    print(f"wrote p1_on_p2mesh_obstacle{TAG}.npz / .json")


if __name__ == "__main__":
    main()
