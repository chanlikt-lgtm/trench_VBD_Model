#!/usr/bin/env python3
"""Spatially decomposed Bernstein bounds for the P2 neutral area.

Regions are fixed physical rectangles in silicon, independent of mesh
connectivity:
    A  : 2.32 <= x <= 2.37,  |y| <= 0.13
    B  : 2.395 <= x <= 2.420, |y| <= 0.24
    REM: charged Si outside A and B

The subdivision carries physical coordinates alongside the parent-frame
barycentrics, so a sub-triangle can be classified as inside / outside /
straddling the region. Lower bounds count only sub-triangles fully inside the
region whose Bernstein coefficients are all <= 0; upper bounds additionally
carry every straddling or unresolved sub-triangle at full area. Both are
therefore valid one-sided bounds.
"""
from __future__ import annotations
import json
import numpy as np

import p2_audit as A
from p2_bernstein import B
from p2_area_bounds import eval_bern, sub_coeff, AREA_TOL, MAX_DEPTH, EPS

M = A.M

REGIONS = {
    "A":   (2.32, 2.37, -0.13, 0.13),
    "B":   (2.395, 2.420, -0.24, 0.24),
}


def classify(ph, rect):
    x0, x1, y0, y1 = rect
    P = np.asarray(ph, dtype=float).reshape(-1, 2)
    xs = P[:, 0]; ys = P[:, 1]
    if xs.max() <= x0 or xs.min() >= x1 or ys.max() <= y0 or ys.min() >= y1:
        return -1
    if xs.min() >= x0 and xs.max() <= x1 and ys.min() >= y0 and ys.max() <= y1:
        return 1
    return 0


def mid3(P, Q):
    return ((P[0]+Q[0])/2, (P[1]+Q[1])/2, (P[2]+Q[2])/2)


def mid2(P, Q):
    return ((P[0]+Q[0])/2, (P[1]+Q[1])/2)


def recurse(pc, bar, ph, area, depth, rect, acc):
    rel = classify(ph, rect)
    if rel < 0:
        return
    sc = sub_coeff(pc, bar[0], bar[1], bar[2])
    if rel > 0:
        if all(ci <= 0.0 for ci in sc):
            acc[0] += area
            return
        if all(ci > 0.0 for ci in sc):
            return
    if area < AREA_TOL or depth >= MAX_DEPTH:
        acc[1] += area
        return
    M01b = mid3(bar[0], bar[1]); M12b = mid3(bar[1], bar[2]); M20b = mid3(bar[2], bar[0])
    m01 = mid2(ph[0], ph[1]); m12 = mid2(ph[1], ph[2]); m20 = mid2(ph[2], ph[0])
    kids_b = ((bar[0], M01b, M20b), (M01b, bar[1], M12b),
              (M20b, M12b, bar[2]), (M01b, M12b, M20b))
    kids_p = ((ph[0], m01, m20), (m01, ph[1], m12),
              (m20, m12, ph[2]), (m01, m12, m20))
    for kb, kp in zip(kids_b, kids_p):
        recurse(pc, kb, kp, area/4.0, depth+1, rect, acc)


def region_bounds(nodes, tris, mat, u, rect):
    area = np.abs(A.element_geom(nodes, tris)[2])
    cxs = nodes[tris][:, :3, 0].mean(axis=1)
    charged = (mat == M.SI) & (cxs > A.BODY_BOTTOM + 1e-12)
    lower = 0.0; unresolved = 0.0
    for e in np.flatnonzero(charged):
        P = nodes[tris[e]][:, :3]
        if classify(P, rect) < 0:
            continue
        c = B @ u[tris[e]] - EPS
        acc = [0.0, 0.0]
        bar0 = ((1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0))
        recurse(c, bar0, (tuple(P[0]), tuple(P[1]), tuple(P[2])),
                area[e], 0, rect, acc)
        lower += acc[0]; unresolved += acc[1]
    return lower, lower + unresolved


def remainder_bounds(nodes, tris, mat, u):
    """Charged Si outside A and B, by masked whole-domain recursion."""
    area = np.abs(A.element_geom(nodes, tris)[2])
    cxs = nodes[tris][:, :3, 0].mean(axis=1)
    charged = (mat == M.SI) & (cxs > A.BODY_BOTTOM + 1e-12)
    lower = 0.0; unresolved = 0.0
    for e in np.flatnonzero(charged):
        P = nodes[tris[e]][:, :3]
        c = B @ u[tris[e]] - EPS
        acc = [0.0, 0.0]
        bar0 = ((1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0))

        def rec(pc, bar, ph, ar, depth):
            ca = classify(ph, REGIONS["A"]); cb = classify(ph, REGIONS["B"])
            if ca > 0 or cb > 0:
                return                      # belongs to A or B
            straddle = (ca == 0) or (cb == 0)
            sc = sub_coeff(pc, bar[0], bar[1], bar[2])
            if not straddle:
                if all(ci <= 0.0 for ci in sc):
                    acc[0] += ar; return
                if all(ci > 0.0 for ci in sc):
                    return
            if ar < AREA_TOL or depth >= MAX_DEPTH:
                acc[1] += ar; return
            M01b = mid3(bar[0], bar[1]); M12b = mid3(bar[1], bar[2]); M20b = mid3(bar[2], bar[0])
            m01 = mid2(ph[0], ph[1]); m12 = mid2(ph[1], ph[2]); m20 = mid2(ph[2], ph[0])
            for kb, kp in (((bar[0], M01b, M20b), (ph[0], m01, m20)),
                           ((M01b, bar[1], M12b), (m01, ph[1], m12)),
                           ((M20b, M12b, bar[2]), (m20, m12, ph[2])),
                           ((M01b, M12b, M20b), (m01, m12, m20))):
                rec(pc, kb, kp, ar/4.0, depth+1)

        rec(c, bar0, (tuple(P[0]), tuple(P[1]), tuple(P[2])), area[e], 0)
        lower += acc[0]; unresolved += acc[1]
    return lower, lower + unresolved


def main():
    meshes = [("mid", "p2_p3_rmid_polished.npz"),
              ("L1b", "p2_L1b_scaled.npz"),
              ("M2",  "p2_M2_scaled.npz"),
              ("L1",  "p2_L1_scaled.npz")]
    out = {}
    for tag, fn in meshes:
        z = np.load(fn)
        nodes, tris, mat, u = z["nodes"], z["tris"], z["mat"], z["u"]
        rec = {}
        for name, rect in REGIONS.items():
            lo, up = region_bounds(nodes, tris, mat, u, rect)
            rec[name] = [lo, up]
            print(f"[{tag}] {name:4s} lower={lo:.12f} upper={up:.12f} gap={up-lo:.2e}")
        lo, up = remainder_bounds(nodes, tris, mat, u)
        rec["REM"] = [lo, up]
        print(f"[{tag}] REM  lower={lo:.12f} upper={up:.12f} gap={up-lo:.2e}")
        s_lo = sum(v[0] for v in rec.values()); s_up = sum(v[1] for v in rec.values())
        rec["SUM"] = [s_lo, s_up]
        print(f"[{tag}] SUM  lower={s_lo:.12f} upper={s_up:.12f}")
        out[tag] = rec
    json.dump(out, open("p2_area_regions.json", "w"), indent=2)
    print("\nwrote p2_area_regions.json")


if __name__ == "__main__":
    main()
