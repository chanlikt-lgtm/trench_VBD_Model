#!/usr/bin/env python3
"""Wall statistics excluding the gate-corner box by LOCATION, so all seven
high-mismatch samples are removed regardless of pinning state."""
from __future__ import annotations
import numpy as np
import p2_interface as PI
import p2_audit as A

EPS_SI, EPS_OX = A.EPS_SI, A.EPS_OX


def collect(fn):
    z = np.load(fn)
    nodes, tris, mat, phi = z["nodes"], z["tris"], z["mat"], z["phi"]
    _, _, _, g1, g2, g3 = A.element_geom(nodes, tris)
    rows = []
    for ed in PI.build_edge_map(nodes, tris, mat):
        a, b, m = ed["nodes"]
        p0, pm, p1 = nodes[a], nodes[m], nodes[b]
        for t, w in zip(PI.G1, PI.GW):
            pt, dpt = PI.edge_points(p0, pm, p1, t)
            ds = np.linalg.norm(dpt)*w
            tn = dpt/np.linalg.norm(dpt)
            nrm = np.array([tn[1], -tn[0]])
            cen = nodes[tris[ed["si"]][:3]].mean(axis=0)
            if np.dot(cen - pt, nrm) > 0:
                nrm = -nrm
            vals = []
            for e in (ed["si"], ed["ox"]):
                P = nodes[tris[e]][:3]
                Mx = np.array([[P[1,0]-P[0,0], P[2,0]-P[0,0]],
                               [P[1,1]-P[0,1], P[2,1]-P[0,1]]])
                st = np.linalg.solve(Mx, pt - P[0])
                G = A.shape_grads(st[0], st[1], g1[e], g2[e], g3[e])
                vals.append(np.array([G[:,0] @ phi[tris[e]], G[:,1] @ phi[tris[e]]]))
            rows.append(dict(x=float(pt[0]), y=float(pt[1]), ds=float(ds),
                             d_si=EPS_SI*float(vals[0] @ nrm),
                             d_ox=EPS_OX*float(vals[1] @ nrm)))
    return rows


def R(rows):
    w = np.array([r["ds"] for r in rows]); ds = np.array([r["d_si"] for r in rows])
    do = np.array([r["d_ox"] for r in rows])
    return float(np.sum(w*np.abs(ds-do))/np.sum(w*(np.abs(ds)+np.abs(do))/2))


def main():
    for tag, fn in (("r3", "p2_r3_solution.npz"), ("r4", "p2_r4_solution.npz")):
        rows = collect(fn)
        wall = [r for r in rows if r["x"] <= 2.04]
        box = [r for r in wall if abs(r["x"]-0.54) <= 0.03 and abs(r["y"]) >= 0.34]
        rest = [r for r in wall if not (abs(r["x"]-0.54) <= 0.03 and abs(r["y"]) >= 0.34)]
        print(f"[{tag}] wall n={len(wall)} R={R(wall):.6f} | "
              f"corner-box n={len(box)} | wall excl. box n={len(rest)} R={R(rest):.6f}")
        xs = [r["x"] for r in box]
        if xs:
            print(f"      box x range {min(xs):.4f}..{max(xs):.4f}")


if __name__ == "__main__":
    main()
