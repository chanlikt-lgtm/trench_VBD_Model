#!/usr/bin/env python3
"""Test the hypothesis for the r4 wall tail: interface quadrature points whose
Si-side normal flux is identically zero sit on Si elements that are entirely
Dirichlet-pinned (an equipotential contact), where the no-free-charge interface
condition is replaced by electrode charge balance.

Replicates the frozen reconstruction path exactly (same edge map, same 2-point
rule, same normal orientation) and additionally records the adjacent element.
"""
from __future__ import annotations
import json
import numpy as np
import p2_interface as PI
import p2_audit as A

EPS_SI, EPS_OX = A.EPS_SI, A.EPS_OX


def main():
    with open("p2_nearzero.json") as fh:
        dref = json.load(fh)["reference"]["d_ref"]

    z = np.load("p2_r4_solution.npz")
    nodes, tris, mat, phi, fixed = z["nodes"], z["tris"], z["mat"], z["phi"], z["fixed"]
    p, det, area, g1, g2, g3 = A.element_geom(nodes, tris)
    edges = PI.build_edge_map(nodes, tris, mat)

    rows = []
    for ed in edges:
        a, b, m = ed["nodes"]
        p0, pm, p1 = nodes[a], nodes[m], nodes[b]
        for t, w in zip(PI.G1, PI.GW):
            pt, dpt = PI.edge_points(p0, pm, p1, t)
            tn = dpt / np.linalg.norm(dpt)
            nrm = np.array([tn[1], -tn[0]])
            cen_si = nodes[tris[ed["si"]][:3]].mean(axis=0)
            if np.dot(cen_si - pt, nrm) > 0:
                nrm = -nrm
            vals = []
            for e in (ed["si"], ed["ox"]):
                P = nodes[tris[e]][:3]
                Mx = np.array([[P[1, 0]-P[0, 0], P[2, 0]-P[0, 0]],
                               [P[1, 1]-P[0, 1], P[2, 1]-P[0, 1]]])
                st = np.linalg.solve(Mx, pt - P[0])
                vals.append(PI.grad_at(nodes, tris, phi, e, st[0], st[1], g1, g2, g3))
            d_si = EPS_SI * float(vals[0] @ nrm)
            d_ox = EPS_OX * float(vals[1] @ nrm)
            rows.append((pt[0], pt[1], d_si, d_ox, ed["si"],
                         "arc" if pt[0] > 2.04 else "wall"))

    x = np.array([r[0] for r in rows]); y = np.array([r[1] for r in rows])
    ds = np.array([r[2] for r in rows]); do = np.array([r[3] for r in rows])
    si_el = np.array([r[4] for r in rows])
    reg = np.array([r[5] for r in rows])
    dscale = (np.abs(ds) + np.abs(do)) / 2
    keep = dscale >= 1e-3 * dref
    rel = np.abs(ds - do) / np.maximum(2 * dscale, 1e-30)

    bad = keep & (rel > 0.99) & (reg == "wall")
    zero_si = np.abs(ds) < 1e-9

    print(f"replicated samples={len(rows)}  bad={int(bad.sum())}  "
          f"d_si==0: {int((zero_si & (reg=='wall')).sum())} wall")

    # Are the Si elements behind the bad points fully pinned?
    all_pinned = np.array([bool(np.all(fixed[tris[e]])) for e in si_el])
    any_pinned = np.array([bool(np.any(fixed[tris[e]])) for e in si_el])
    print(f"bad points: si element fully Dirichlet-pinned = {int(all_pinned[bad].sum())}/{int(bad.sum())}")
    print(f"bad points: si element partly pinned          = {int(any_pinned[bad].sum())}/{int(bad.sum())}")
    print(f"bad points: phi range on si element = "
          f"{phi[tris[si_el[bad]]].min():.3e} .. {phi[tris[si_el[bad]]].max():.3e}")

    print(f"bad x range {x[bad].min():.4f}..{x[bad].max():.4f}  "
          f"y range {y[bad].min():.4f}..{y[bad].max():.4f}")
    # where does the body-contact Dirichlet block end?
    si_nodes = np.unique(tris[mat == PI.M.SI])
    bod = si_nodes[(nodes[si_nodes, 0] <= A.BODY_BOTTOM + 1e-12)]
    print(f"body-contact Si nodes: {len(bod)}, max x = {nodes[bod,0].max():.4f} "
          f"(BODY_BOTTOM={A.BODY_BOTTOM})")

    # how many bad points have an all-pinned Si element AND phi exactly 0 there
    z_el = np.array([np.all(phi[tris[e]] == 0.0) for e in si_el])
    print(f"bad points with si element all-zero potential = {int(z_el[bad].sum())}/{int(bad.sum())}")

    # same accounting on r3 for comparison
    z3 = np.load("p2_r3_solution.npz")
    n3, t3, m3, ph3, fx3 = z3["nodes"], z3["tris"], z3["mat"], z3["phi"], z3["fixed"]
    rows3 = []
    p3, det3, a3, h1, h2, h3 = A.element_geom(n3, t3)
    for ed in PI.build_edge_map(n3, t3, m3):
        a, b, m = ed["nodes"]
        p0, pm, p1 = n3[a], n3[m], n3[b]
        for t, w in zip(PI.G1, PI.GW):
            pt, dpt = PI.edge_points(p0, pm, p1, t)
            tn = dpt / np.linalg.norm(dpt)
            nrm = np.array([tn[1], -tn[0]])
            cen = n3[t3[ed["si"]][:3]].mean(axis=0)
            if np.dot(cen - pt, nrm) > 0:
                nrm = -nrm
            vals = []
            for e in (ed["si"], ed["ox"]):
                P = n3[t3[e]][:3]
                Mx = np.array([[P[1,0]-P[0,0], P[2,0]-P[0,0]],
                               [P[1,1]-P[0,1], P[2,1]-P[0,1]]])
                st = np.linalg.solve(Mx, pt - P[0])
                vals.append(PI.grad_at(n3, t3, ph3, e, st[0], st[1], h1, h2, h3))
            rows3.append((pt[0], EPS_SI*float(vals[0] @ nrm),
                          EPS_OX*float(vals[1] @ nrm), ed["si"],
                          "arc" if pt[0] > 2.04 else "wall"))
    ds3 = np.array([r[1] for r in rows3]); do3 = np.array([r[2] for r in rows3])
    e3 = np.array([r[3] for r in rows3]); rg3 = np.array([r[4] for r in rows3])
    ds3c = (np.abs(ds3)+np.abs(do3))/2
    k3 = ds3c >= 1e-3*dref
    r3rel = np.abs(ds3-do3)/np.maximum(2*ds3c, 1e-30)
    b3 = k3 & (r3rel > 0.99) & (rg3 == "wall")
    pin3 = np.array([bool(np.all(fx3[t3[e]])) for e in e3])
    print(f"[r3] wall bad={int(b3.sum())} of kept {int((k3&(rg3=='wall')).sum())}, "
          f"fully pinned={int(pin3[b3].sum())}")


if __name__ == "__main__":
    main()
