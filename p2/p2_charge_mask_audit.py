#!/usr/bin/env python3
"""Audit the charged-Si mask on L1 vs L1b.

charged = (mat == SI) & (element centroid x > BODY_BOTTOM)

The centroid test is a whole-element classification, so elements straddling
x = BODY_BOTTOM are assigned entirely to one side. This script measures that
error exactly by clipping, and checks whether the L1/L1b charged-area
difference is explained by it rather than by a moved mask.
"""
from __future__ import annotations
import json
import numpy as np

import p2_audit as A

M = A.M
XB = A.BODY_BOTTOM


def clip_area(pts, xb):
    """Area of the part of a convex polygon with x >= xb."""
    out = []
    n = len(pts)
    for i in range(n):
        cur = pts[i]; nxt = pts[(i+1) % n]
        cin = cur[0] >= xb; nin = nxt[0] >= xb
        if cin:
            out.append(cur)
        if cin != nin:
            t = (xb - cur[0])/(nxt[0] - cur[0])
            out.append((xb, cur[1] + t*(nxt[1] - cur[1])))
    if len(out) < 3:
        return 0.0
    a = 0.0
    for i in range(len(out)):
        x1, y1 = out[i]; x2, y2 = out[(i+1) % len(out)]
        a += x1*y2 - x2*y1
    return abs(a)/2


def audit(path, label):
    z = np.load(path)
    nodes, tris, mat = z["nodes"], z["tris"], z["mat"]
    si = mat == M.SI
    t = tris[si]
    p = nodes[t][:, :3]
    det = ((p[:, 1, 0]-p[:, 0, 0])*(p[:, 2, 1]-p[:, 0, 1]) -
           (p[:, 2, 0]-p[:, 0, 0])*(p[:, 1, 1]-p[:, 0, 1]))
    area = 0.5*np.abs(det)
    cxs = p[:, :, 0].mean(1)
    charged = cxs > XB
    xmin = p[:, :, 0].min(1); xmax = p[:, :, 0].max(1)
    straddle = (xmin <= XB) & (xmax > XB)

    exact_ch = 0.0
    for i in np.flatnonzero(straddle):
        exact_ch += clip_area([tuple(q) for q in p[i]], XB)
    exact_ch += area[cxs <= XB].sum()*0  # non-straddling below line contribute 0
    exact_ch += area[~straddle & (xmin > XB)].sum()

    classified = area[charged].sum()
    res = dict(
        label=label, si_elements=int(si.sum()), si_area=float(area.sum()),
        charged_elements=int(charged.sum()), charged_area=float(classified),
        straddle_elements=int(straddle.sum()),
        straddle_area=float(area[straddle].sum()),
        straddle_above=float(area[straddle & charged].sum()),
        straddle_below=float(area[straddle & ~charged].sum()),
        exact_charged=float(exact_ch), exact_body=float(area.sum() - exact_ch),
        # signed classification error: classified minus exact
        error=float(classified - exact_ch),
        body_bottom=float(XB),
    )
    print(f"\n=== {label} ===")
    print(f"  Si elements {res['si_elements']:7d}  area {res['si_area']:.12f}")
    print(f"  charged (centroid test) {res['charged_elements']:7d}  "
          f"area {res['charged_area']:.12f}")
    print(f"  straddling x={XB}: {res['straddle_elements']:5d} elements, "
          f"area {res['straddle_area']:.3e} "
          f"(above {res['straddle_above']:.3e} / below {res['straddle_below']:.3e})")
    print(f"  exact charged {res['exact_charged']:.12f}  exact body {res['exact_body']:.12f}")
    print(f"  classification error (classified - exact) = {res['error']:+.3e}")
    return res


def main():
    a = audit("p2_L1_p2mesh.npz", "L1")
    b = audit("p2_L1b_p2mesh.npz", "L1b")
    print("\n=== comparison ===")
    print(f"  Si area          L1 {a['si_area']:.12f}  L1b {b['si_area']:.12f}  "
          f"diff {b['si_area']-a['si_area']:+.3e}")
    print(f"  classified ch.   L1 {a['charged_area']:.12f}  L1b {b['charged_area']:.12f}  "
          f"diff {b['charged_area']-a['charged_area']:+.3e}")
    print(f"  exact charged    L1 {a['exact_charged']:.12f}  L1b {b['exact_charged']:.12f}  "
          f"diff {b['exact_charged']-a['exact_charged']:+.3e}")
    print(f"  mask error       L1 {a['error']:+.3e}  L1b {b['error']:+.3e}  "
          f"diff {b['error']-a['error']:+.3e}")
    obs = b["charged_area"] - a["charged_area"]
    exp = b["error"] - a["error"]
    print(f"\n  observed charged-area change  {obs:+.6e}")
    print(f"  predicted from mask error     {exp:+.6e}")
    print(f"  residual                      {obs-exp:+.3e}")
    json.dump(dict(L1=a, L1b=b, observed_delta=obs, predicted_delta=exp,
                   residual=obs-exp), open("p2_charge_mask_audit.json", "w"), indent=2)
    print("wrote p2_charge_mask_audit.json")


if __name__ == "__main__":
    main()
