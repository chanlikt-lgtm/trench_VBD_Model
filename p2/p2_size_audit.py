#!/usr/bin/env python3
"""Compare realised element-size fields across L1 / L1b / intermediate.

h is reported as the equilateral-equivalent side length h = sqrt(4A/sqrt(3)).
"""
from __future__ import annotations
import json
import numpy as np

BOXES = {
    "A_core": (2.335, 2.355, -0.05, 0.05),
    "A_full": (2.32, 2.37, -0.13, 0.13),
    "B_core": (2.400, 2.415, -0.10, 0.10),
    "B_full": (2.395, 2.420, -0.24, 0.24),
    "far":    (1.00, 1.50, -0.60, 0.60),
    "mid_far": (2.10, 2.25, -0.60, 0.60),
}
MESHES = {
    "L1":  "p2_L1_p2mesh.npz",
    "L1b": "p2_L1b_p2mesh.npz",
    "mid": "p2_rmid_p2mesh.npz",
}


def stats(path):
    z = np.load(path)
    nodes, tris = z["nodes"], z["tris"]
    p = nodes[tris][:, :3]
    det = ((p[:, 1, 0]-p[:, 0, 0])*(p[:, 2, 1]-p[:, 0, 1]) -
           (p[:, 2, 0]-p[:, 0, 0])*(p[:, 1, 1]-p[:, 0, 1]))
    area = 0.5*np.abs(det)
    h = np.sqrt(4*area/np.sqrt(3))
    cx = p[:, :, 0].mean(1); cy = p[:, :, 1].mean(1)
    out = {}
    for k, (x0, x1, y0, y1) in BOXES.items():
        m = (cx >= x0) & (cx <= x1) & (cy >= y0) & (cy <= y1)
        if m.sum() < 5:
            out[k] = None
            continue
        out[k] = dict(n=int(m.sum()),
                      p10=float(np.percentile(h[m], 10)),
                      med=float(np.median(h[m])),
                      p90=float(np.percentile(h[m], 90)))
    return out, len(nodes), len(tris)


def main():
    res = {}
    for name, path in MESHES.items():
        try:
            res[name], nn, ne = stats(path)
        except FileNotFoundError:
            print(f"{name}: missing {path}")
            continue
        print(f"\n=== {name}  nodes={nn} elems={ne} ===")
        for k, v in res[name].items():
            if v is None:
                print(f"  {k:8s}: (too few)")
            else:
                print(f"  {k:8s}: n={v['n']:7d}  h p10/med/p90 = "
                      f"{v['p10']:.5f} / {v['med']:.5f} / {v['p90']:.5f}")
    print("\n=== contrast ratios (med in region / med far) ===")
    for name in res:
        far = res[name].get("far")
        if not far:
            continue
        row = []
        for k in ("A_core", "A_full", "B_core", "B_full"):
            v = res[name].get(k)
            row.append(f"{k}={v['med']/far['med']:.3f}" if v else f"{k}=n/a")
        print(f"  {name:4s} far_h={far['med']:.5f}  " + "  ".join(row))
    json.dump(res, open("p2_size_audit.json", "w"), indent=2)


if __name__ == "__main__":
    main()

