#!/usr/bin/env python3
"""Split the arc region at the x=2.15 arc/cap junction.

The frozen profile is tangent-continuous at x=2.05 (2.4 deg kink) but has a
79.5 deg tangent discontinuity at x=2.15: the arc arrives vertical, the cap
leaves at atan(0.185) = 10.5 deg. That corner is where the chord normal is
meaningless. Report arc-proper and corner separately.
"""
from __future__ import annotations
import json
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
    if not rows:
        return None
    w = np.array([r["ds"] for r in rows])
    ds = np.array([r["d_si"] for r in rows])
    do = np.array([r["d_ox"] for r in rows])
    return float(np.sum(w*np.abs(ds-do)) / np.sum(w*(np.abs(ds)+np.abs(do))/2))


def main():
    out = {}
    for tag, fn in (("r3", "p2_r3_solution.npz"), ("r4", "p2_r4_solution.npz")):
        rows = collect(fn)
        arc = [r for r in rows if r["x"] > 2.04]
        proper = [r for r in arc if r["x"] <= 2.14]
        corner = [r for r in arc if r["x"] > 2.14]
        vals = dict(arc=R(arc), arc_proper=R(proper), corner=R(corner),
                    n_arc=len(arc), n_proper=len(proper), n_corner=len(corner))
        out[tag] = vals
        print(f"[{tag}] arc n={vals['n_arc']} R_abs={vals['arc']:.6f} | "
              f"proper n={vals['n_proper']} R_abs={vals['arc_proper']:.6f} | "
              f"corner(x>2.14) n={vals['n_corner']} R_abs={vals['corner']:.6f}")

    # corner angle from the frozen constants
    import math
    slope_cap = (0.575 - 0.525) / (2.42 - 2.15)
    print(f"\n[frozen profile] arc tangent at x=2.15: 90.0 deg (vertical)")
    print(f"[frozen profile] cap slope: {math.degrees(math.atan(slope_cap)):.2f} deg")
    print(f"[frozen profile] corner discontinuity: "
          f"{90.0 - math.degrees(math.atan(slope_cap)):.2f} deg")
    slope_taper = (0.425 - 0.34) / 2.05
    print(f"[frozen profile] taper slope at x=2.05: {math.degrees(math.atan(slope_taper)):.2f} deg")
    print(f"[frozen profile] junction kink at x=2.05: "
          f"{math.degrees(math.atan(slope_taper)):.2f} deg")

    with open("p2_arc_split.json", "w") as fh:
        json.dump(out, fh, indent=2)


if __name__ == "__main__":
    main()
