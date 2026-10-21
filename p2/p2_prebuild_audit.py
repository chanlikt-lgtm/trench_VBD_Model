#!/usr/bin/env python3
"""Pre-build audit for the controlled area-convergence study.

1. Realised element size of the elements that straddle the neutral boundary
   (u = eps), per mesh — the quantity that should control area convergence,
   as opposed to sizes inside the analysis windows.
2. Element-level neutral area binned along x, to locate where the L1b -> L1
   increase sits inside the remainder region.
"""
from __future__ import annotations
import json
import numpy as np

import p2_audit as A
from p2_bernstein import B

M = A.M
EPS = 1e-8
REG = {"A": (2.32, 2.37, -0.13, 0.13), "B": (2.395, 2.420, -0.24, 0.24)}


def in_region(cx, cy, r):
    return (cx >= r[0]) & (cx <= r[1]) & (cy >= r[2]) & (cy <= r[3])


def audit(fn):
    z = np.load(fn)
    nodes, tris, mat, u = z["nodes"], z["tris"], z["mat"], z["u"]
    area = np.abs(A.element_geom(nodes, tris)[2])
    P = nodes[tris][:, :3]
    cx = P[:, :, 0].mean(1); cy = P[:, :, 1].mean(1)
    charged = (mat == M.SI) & (cx > A.BODY_BOTTOM + 1e-12)
    h = np.sqrt(4.0*area/np.sqrt(3.0))

    c = (B @ u[tris].T).T - EPS          # (nel, 6) Bernstein coefficients
    allneg = (c <= 0).all(1)
    allpos = (c > 0).all(1)
    mixed = charged & ~allneg & ~allpos
    neutral = charged & allneg

    reg = np.where(in_region(cx, cy, REG["A"]), "A",
          np.where(in_region(cx, cy, REG["B"]), "B", "REM"))

    out = {"boundary": {}, "neutral_area_by_band": {}}
    for tag, m in (("A", reg == "A"), ("B", reg == "B"), ("REM", reg == "REM")):
        sel = mixed & m
        if sel.sum() < 3:
            out["boundary"][tag] = None
            continue
        hh = h[sel]
        out["boundary"][tag] = dict(n=int(sel.sum()),
                                    p10=float(np.percentile(hh, 10)),
                                    med=float(np.median(hh)),
                                    p90=float(np.percentile(hh, 90)),
                                    area=float(area[sel].sum()))

    edges = np.arange(0.50, 2.60, 0.05)
    bins = []
    for a, b in zip(edges[:-1], edges[1:]):
        m = neutral & (cx >= a) & (cx < b) & (reg == "REM")
        bins.append((float(a), float(b), float(area[m].sum())))
    out["neutral_area_by_band"]["REM"] = bins
    out["totals"] = dict(mixed=int(mixed.sum()), neutral=int(neutral.sum()),
                         neutral_area=float(area[neutral].sum()))
    return out


def main():
    meshes = [("mid", "p2_p3_rmid_polished.npz"),
              ("L1b", "p2_L1b_scaled.npz"),
              ("M2", "p2_M2_scaled.npz"),
              ("L1", "p2_L1_scaled.npz")]
    res = {}
    for tag, fn in meshes:
        r = audit(fn)
        res[tag] = r
        print(f"\n[{tag}] boundary elements: "
              + "  ".join(f"{k}:{v['n'] if v else '-'} "
                          f"(h_med {v['med']:.5f})" if v else f"{k}:-"
                          for k, v in r["boundary"].items()))
    print("\nREM element-level neutral area by x band (L1 - L1b)")
    b1 = res["L1"]["neutral_area_by_band"]["REM"]
    b0 = res["L1b"]["neutral_area_by_band"]["REM"]
    rows = []
    for (a, e, va), (_, _, vb) in zip(b1, b0):
        rows.append((a, e, vb, va, va - vb))
    rows.sort(key=lambda t: -abs(t[4]))
    print(f"{'x band':>14s} {'L1b':>12s} {'L1':>12s} {'delta':>12s}")
    for a, e, vb, va, d in rows[:10]:
        print(f"{a:5.2f}-{e:4.2f}   {vb:12.3e} {va:12.3e} {d:+12.3e}")
    tot = sum(r[4] for r in rows)
    print(f"REM total delta (element-level): {tot:+.3e}")
    srt = sorted(rows, key=lambda t: -t[4])
    print(f"largest positive band: x {srt[0][0]:.2f}-{srt[0][1]:.2f}  {srt[0][4]:+.3e}"
          f"  ({100*srt[0][4]/tot:.0f}% of REM delta)")
    json.dump(res, open("p2_prebuild_audit.json", "w"), indent=2)
    print("\nwrote p2_prebuild_audit.json")


if __name__ == "__main__":
    main()
