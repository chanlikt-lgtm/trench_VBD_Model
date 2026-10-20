#!/usr/bin/env python3
"""P1 linear elements on the P2 mesh's own vertex triangles.

This removes the geometry/mesh difference from the comparison: same nodes, same
triangles, same material tags, same body treatment (uncharged + all nodes
pinned), same Dirichlet sets. Only the element order and load distribution
differ from the P2 variant.
"""
from __future__ import annotations
import json
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

import p2_audit as A

M = A.M
V_D = A.V_DRAIN
EPS = {M.SI: A.EPS_SI, M.OX: A.EPS_OX, M.GATE: A.EPS_SI, M.SHIELD: A.EPS_SI}


def main():
    z = np.load("p2_r3_solution.npz")
    nodes6, tris6, mat = z["nodes"], z["tris"], z["mat"]
    x6 = nodes6[:, 0]

    # vertex-only mesh
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
    e1 = p[:, 1] - p[:, 0]
    e2 = p[:, 2] - p[:, 0]
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
    si_nodes = np.unique(tris[si])
    src = si_nodes[x[si_nodes] <= 1e-12]
    drn = np.flatnonzero(x >= M.BOTTOM_DEPTH - 1e-9)

    fixed[grounded] = True; val[grounded] = V_D
    fixed[src] = True; val[src] = V_D
    fixed[drn] = True; val[drn] = 0.0

    idx = np.flatnonzero(~fixed)
    u = val.copy()
    u[idx] = spla.spsolve(K[idx][:, idx].tocsc(), (-f - K @ val)[idx])
    phi = V_D - u
    free_si = si_nodes[~fixed[si_nodes]]
    print(f"[P1onP2] nodes={len(nodes)} tris={len(tris)} charged area={area[charged].sum():.9f} "
          f"charge={A.RHO*area[charged].sum():.6f}")
    print(f"[P1onP2] max(phi)-V_D = {(phi[si_nodes]-V_D).max():+.9f}  "
          f"min u free Si = {u[free_si].min():+.9f} at "
          f"(x={nodes[free_si[np.argmin(u[free_si])],0]:.5f}, "
          f"y={nodes[free_si[np.argmin(u[free_si])],1]:.5f})")

    np.savez_compressed("p1_on_p2mesh.npz", nodes=nodes, tris=tris, mat=mat,
                        phi=phi, u=u, fixed=fixed)

    # probe comparison against the P2 variant on identical geometry
    v = np.load("p2_r3_variant.npz")
    nodes2, tris2, phi2 = v["nodes"], v["tris"], v["phi"]
    xs = np.linspace(2.25, 2.415, 12)
    ys = np.array([-0.30, -0.20, -0.10, 0.0, 0.10, 0.20, 0.30])
    pts = np.array([[a, bb] for a in xs for bb in ys])

    def locate(nodes_, tris_, pts_):
        Aa = nodes_[tris_[:, 0]]; Bb = nodes_[tris_[:, 1]]; Cc = nodes_[tris_[:, 2]]
        det_ = ((Bb[:, 0]-Aa[:, 0])*(Cc[:, 1]-Aa[:, 1]) -
                (Cc[:, 0]-Aa[:, 0])*(Bb[:, 1]-Aa[:, 1]))
        el = np.full(len(pts_), -1, int); bar = np.zeros((len(pts_), 2))
        for k, q in enumerate(pts_):
            L1 = ((Bb[:, 1]-Cc[:, 1])*(q[0]-Cc[:, 0]) +
                  (Cc[:, 0]-Bb[:, 0])*(q[1]-Cc[:, 1]))/det_
            L2 = ((Cc[:, 1]-Aa[:, 1])*(q[0]-Cc[:, 0]) +
                  (Aa[:, 0]-Cc[:, 0])*(q[1]-Cc[:, 1]))/det_
            L3 = 1-L1-L2
            m = (L1 >= -1e-9) & (L2 >= -1e-9) & (L3 >= -1e-9)
            if m.any():
                j = int(np.flatnonzero(m)[0])
                el[k] = j; bar[k] = (L1[j], L2[j])
        return el, bar

    e1_, b1_ = locate(nodes, tris, pts)
    e2_, b2_ = locate(nodes2, tris2, pts)
    p1v = np.full(len(pts), np.nan); p2v = np.full(len(pts), np.nan)
    for k in range(len(pts)):
        if e1_[k] >= 0:
            L1, L2 = b1_[k]
            ids = tris[e1_[k]]
            p1v[k] = L1*phi[ids[0]] + L2*phi[ids[1]] + (1-L1-L2)*phi[ids[2]]
        if e2_[k] >= 0:
            L1, L2 = b2_[k]
            p2v[k] = float(np.dot(phi2[tris2[e2_[k]]], A.shape_vals(L1, L2)))
    ok = ~np.isnan(p1v) & ~np.isnan(p2v)
    d = p2v[ok] - p1v[ok]
    print(f"[same-mesh] probes={len(pts)} matched={int(ok.sum())}")
    print(f"  dphi (P2 - P1lin): max|.|={np.abs(d).max():.6e} rms={np.sqrt(np.mean(d*d)):.6e} "
          f"mean={d.mean():+.6e}")
    print(f"  near drain (x>2.38): max|.|={np.abs(d[np.repeat(xs, len(ys))[ok] > 2.38]).max():.6e}"
          if (np.repeat(xs, len(ys))[ok] > 2.38).any() else "")
    json.dump(dict(p1onp2_overshoot=float((phi[si_nodes]-V_D).max()),
                   p1onp2_min_u=float(u[free_si].min()),
                   probes=int(ok.sum()), max_abs=float(np.abs(d).max()),
                   rms=float(np.sqrt(np.mean(d*d))), mean=float(d.mean())),
              open("p1_on_p2mesh.json", "w"), indent=2)


if __name__ == "__main__":
    main()
