#!/usr/bin/env python3
"""Fixed-probe electric-field comparison on the existing r3/r4 solutions.

E = -grad(phi), reconstructed with each field's own interpolant and the same
probe coordinates on both meshes. Three groups:
  regular depleted silicon (away from front and corner)
  near the neutral front
  near but not at the 79.51-degree corner (2.15, +/-0.525)

Also reports u = V_D - phi and whether the probe is neutral (u <= 1e-8).
"""
from __future__ import annotations
import json
import numpy as np
import p2_audit as A

M = A.M
V_D = A.V_DRAIN

PROBES = [
    ("regular", 1.60, 0.00), ("regular", 1.80, 0.10), ("regular", 1.80, -0.10),
    ("regular", 2.00, 0.00), ("regular", 2.05, 0.20), ("regular", 1.90, -0.25),
    ("front", 2.33, 0.00), ("front", 2.35, 0.00), ("front", 2.36, 0.08),
    ("front", 2.34, -0.08), ("front", 2.37, 0.12),
    ("corner", 2.20, 0.40), ("corner", 2.22, 0.45), ("corner", 2.18, -0.42),
    ("corner", 2.25, 0.35), ("corner", 2.12, 0.45),
]


def load_p2(fn):
    z = np.load(fn)
    nodes, tris, mat = z["nodes"], z["tris"], z["mat"]
    phi = z["phi"]
    cxs = nodes[tris][:, :3, 0].mean(axis=1)
    charged = (mat == M.SI) & (cxs > A.BODY_BOTTOM + 1e-12)
    els = np.flatnonzero(charged)
    _, det, area, g1, g2, g3 = A.element_geom(nodes, tris)
    return dict(nodes=nodes, tris=tris, phi=phi, els=els,
                g1=g1, g2=g2, g3=g3)


def load_p1(fn):
    z = np.load(fn)
    nodes, tris, phi = z["nodes"], z["tris"], z["phi"]
    p = nodes[tris]
    det = ((p[:, 1, 0]-p[:, 0, 0])*(p[:, 2, 1]-p[:, 0, 1]) -
           (p[:, 2, 0]-p[:, 0, 0])*(p[:, 1, 1]-p[:, 0, 1]))
    b = np.stack([p[:, 1, 1]-p[:, 2, 1], p[:, 2, 1]-p[:, 0, 1], p[:, 0, 1]-p[:, 1, 1]], 1)
    c = np.stack([p[:, 2, 0]-p[:, 1, 0], p[:, 0, 0]-p[:, 2, 0], p[:, 1, 0]-p[:, 0, 0]], 1)
    return dict(nodes=nodes, tris=tris, phi=phi, det=det, b=b, c=c)


def probe_p2(d, x0, y0):
    nodes, tris, els = d["nodes"], d["tris"], d["els"]
    Aa = nodes[tris[els, 0]]; Bb = nodes[tris[els, 1]]; Cc = nodes[tris[els, 2]]
    det = ((Bb[:, 0]-Aa[:, 0])*(Cc[:, 1]-Aa[:, 1]) -
           (Cc[:, 0]-Aa[:, 0])*(Bb[:, 1]-Aa[:, 1]))
    L1 = ((Bb[:, 1]-Cc[:, 1])*(x0-Cc[:, 0]) + (Cc[:, 0]-Bb[:, 0])*(y0-Cc[:, 1]))/det
    L2 = ((Cc[:, 1]-Aa[:, 1])*(x0-Cc[:, 0]) + (Aa[:, 0]-Cc[:, 0])*(y0-Cc[:, 1]))/det
    L3 = 1-L1-L2
    m = (L1 >= -1e-9) & (L2 >= -1e-9) & (L3 >= -1e-9)
    if not m.any():
        return None
    j = int(np.flatnonzero(m)[0]); e = int(els[j])
    G = A.shape_grads(L1[j], L2[j], d["g1"][e], d["g2"][e], d["g3"][e])
    grad = np.array([G[:, 0] @ d["phi"][tris[e]], G[:, 1] @ d["phi"][tris[e]]])
    u = V_D - float(np.dot(d["phi"][tris[e]], A.shape_vals(L1[j], L2[j])))
    return dict(E=-grad, u=u)


def probe_p1(d, x0, y0):
    nodes, tris = d["nodes"], d["tris"]
    Aa = nodes[tris[:, 0]]; Bb = nodes[tris[:, 1]]; Cc = nodes[tris[:, 2]]
    det = ((Bb[:, 0]-Aa[:, 0])*(Cc[:, 1]-Aa[:, 1]) -
           (Cc[:, 0]-Aa[:, 0])*(Bb[:, 1]-Aa[:, 1]))
    L1 = ((Bb[:, 1]-Cc[:, 1])*(x0-Cc[:, 0]) + (Cc[:, 0]-Bb[:, 0])*(y0-Cc[:, 1]))/det
    L2 = ((Cc[:, 1]-Aa[:, 1])*(x0-Cc[:, 0]) + (Aa[:, 0]-Cc[:, 0])*(y0-Cc[:, 1]))/det
    L3 = 1-L1-L2
    m = (L1 >= -1e-9) & (L2 >= -1e-9) & (L3 >= -1e-9)
    if not m.any():
        return None
    j = int(np.flatnonzero(m)[0]); e = int(j)
    ph = d["phi"][tris[e]]
    grad = np.array([d["b"][e] @ ph, d["c"][e] @ ph])/d["det"][e]
    u = V_D - float(L1[j]*ph[0] + L2[j]*ph[1] + (1-L1[j]-L2[j])*ph[2])
    return dict(E=-grad, u=u)


def main():
    r3p2 = load_p2("p2_p3_r3_polished.npz")
    r4p2 = load_p2("p2_p3_r4_polished.npz")
    r3p1 = load_p1("p1_on_p2mesh_obstacle.npz")
    r4p1 = load_p1("p1_on_p2mesh_obstacle_r4.npz")

    rows = []
    print(f"{'grp':8s} {'x':>5s} {'y':>6s} | "
          f"{'P1r3':>10s} {'P2r3':>10s} {'P1r4':>10s} {'P2r4':>10s} | "
          f"{'u_r3':>9s} {'u_r4':>9s}")
    for grp, x0, y0 in PROBES:
        a = probe_p1(r3p1, x0, y0); b = probe_p2(r3p2, x0, y0)
        c = probe_p1(r4p1, x0, y0); d = probe_p2(r4p2, x0, y0)
        if None in (a, b, c, d):
            print(f"{grp:8s} {x0:5.2f} {y0:6.2f}   outside one of the meshes")
            continue
        E = [np.linalg.norm(z["E"]) for z in (a, b, c, d)]
        rows.append(dict(group=grp, x=x0, y=y0,
                         P1r3=E[0], P2r3=E[1], P1r4=E[2], P2r4=E[3],
                         u_r3=b["u"], u_r4=d["u"]))
        print(f"{grp:8s} {x0:5.2f} {y0:6.2f} | "
              f"{E[0]:10.4f} {E[1]:10.4f} {E[2]:10.4f} {E[3]:10.4f} | "
              f"{b['u']:9.4f} {d['u']:9.4f}")

    print("\n[changes r3 -> r4 in |E|]")
    for r in rows:
        print(f"  {r['group']:8s} ({r['x']:.2f},{r['y']:+.2f}): "
              f"P1 {r['P1r4']-r['P1r3']:+.5f}  P2 {r['P2r4']-r['P2r3']:+.5f}  "
              f"P2-P1(r4) {r['P2r4']-r['P1r4']:+.5f}")

    json.dump(rows, open("p2_field_probes.json", "w"), indent=2)
    print("\nwrote p2_field_probes.json")


if __name__ == "__main__":
    main()
