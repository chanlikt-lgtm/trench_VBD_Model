#!/usr/bin/env python3
"""Tangent error restricted to Gauss points that actually lie on the analytic
circle (| |rad| - R | < 5e-4). Points on the cap or the taper are outside the
circle's domain and are reported separately, not silently averaged in.
"""
from __future__ import annotations
import numpy as np
import p2_interface as PI
import p2_audit as A

R, XC, A_ARC_END = 0.10, 2.05, 0.525


def main():
    z = np.load("p2_r4_solution.npz")
    nodes, tris, mat = z["nodes"], z["tris"], z["mat"]

    on_arc, off_arc = [], []
    expected = []
    for ed in PI.build_edge_map(nodes, tris, mat):
        a, b, m = ed["nodes"]
        p0, pm, p1 = nodes[a], nodes[m], nodes[b]
        if not (XC < 0.5*(p0[0]+p1[0]) < XC+R):
            continue
        c = np.linalg.norm(p1-p0)
        theta_edge = 2*np.arcsin(min(c/(2*R), 1.0))
        for t in PI.G1:
            pt, dpt = PI.edge_points(p0, pm, p1, t)
            tn = dpt/np.linalg.norm(dpt)
            rad = np.array([pt[0]-XC, abs(pt[1])-A_ARC_END])
            rn = np.linalg.norm(rad)
            tan = np.array([-rad[1], rad[0]])/rn
            e = np.degrees(np.arccos(np.clip(abs(np.dot(tn, tan)), -1, 1)))
            if abs(rn - R) < 5e-4:
                on_arc.append(e)
                expected.append(np.degrees(0.289*theta_edge))
            else:
                off_arc.append((abs(rn-R), e))

    on_arc = np.array(on_arc); expected = np.array(expected)
    print(f"on-arc Gauss points n={len(on_arc)}")
    print(f"  measured tangent error: med={np.median(on_arc):.3f} p95={np.percentile(on_arc,95):.3f} max={on_arc.max():.3f} deg")
    print(f"  chord-geometry expectation med={np.median(expected):.3f} deg "
          f"(0.289*theta_edge, theta=c/R)")
    print(f"  ratio measured/expected med={np.median(on_arc/np.maximum(expected,1e-12)):.3f}")
    if off_arc:
        d = np.array([o[0] for o in off_arc]); e = np.array([o[1] for o in off_arc])
        print(f"off-arc points n={len(off_arc)}  | |rad|-R | med={np.median(d):.3e} "
              f"max={d.max():.3e}  their angle med={np.median(e):.2f} max={e.max():.2f} deg")


if __name__ == "__main__":
    main()
