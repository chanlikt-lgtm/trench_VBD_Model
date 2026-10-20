#!/usr/bin/env python3
"""Two follow-ups on the closure checks.

1. det sign per material: is the mixed sign an orientation convention rather
   than degenerate geometry?
2. tangent-angle statistics on arc edges that lie strictly INSIDE the arc,
   excluding edges that touch the x=2.05 or x=2.15 junctions where the analytic
   circle formula does not apply.
"""
from __future__ import annotations
import numpy as np
import p2_interface as PI
import p2_audit as A

R, XC, A_ARC_END = 0.10, 2.05, 0.525
NAMES = {A.M.SI: "Si", A.M.OX: "SiO2", A.M.GATE: "Gate", A.M.SHIELD: "Shield"}


def main():
    z = np.load("p2_r4_solution.npz")
    nodes, tris, mat = z["nodes"], z["tris"], z["mat"]
    _, det, _, _, _, _ = A.element_geom(nodes, tris)

    print("[det sign by material]")
    for code, name in NAMES.items():
        m = mat == code
        if not m.any():
            continue
        print(f"  {name:6s} n={int(m.sum()):6d}  det>0: {int((det[m] > 0).sum()):6d}  "
              f"det<0: {int((det[m] < 0).sum()):6d}  min|det|={np.abs(det[m]).min():.3e}")

    # tangent error on strictly interior arc edges
    ang, ang_all = [], []
    for ed in PI.build_edge_map(nodes, tris, mat):
        a, b, m = ed["nodes"]
        p0, pm, p1 = nodes[a], nodes[m], nodes[b]
        if not (XC < 0.5*(p0[0]+p1[0]) < XC+R):
            continue
        on_junction = any(abs(p[0]-XC) < 1e-7 or abs(p[0]-(XC+R)) < 1e-7 for p in (p0, p1))
        for t in PI.G1:
            pt, dpt = PI.edge_points(p0, pm, p1, t)
            tn = dpt/np.linalg.norm(dpt)
            rad = np.array([pt[0]-XC, abs(pt[1])-A_ARC_END])
            tan = np.array([-rad[1], rad[0]]); tan /= np.linalg.norm(tan)
            e = np.degrees(np.arccos(np.clip(abs(np.dot(tn, tan)), -1, 1)))
            ang_all.append(e)
            if not on_junction:
                ang.append(e)

    ang = np.array(ang); ang_all = np.array(ang_all)
    print(f"\n[arc tangent error] all arc-tagged Gauss points n={len(ang_all)}: "
          f"med={np.median(ang_all):.3f} p95={np.percentile(ang_all,95):.3f} max={ang_all.max():.3f} deg")
    print(f"[arc tangent error] excluding junction edges n={len(ang)}: "
          f"med={np.median(ang):.3f} p95={np.percentile(ang,95):.3f} max={ang.max():.3f} deg")
    if ang.any():
        print(f"  implied normal-direction fractional error sin(theta): "
              f"med={np.sin(np.radians(np.median(ang))):.4f} max={np.sin(np.radians(ang.max())):.4f}")


if __name__ == "__main__":
    main()
