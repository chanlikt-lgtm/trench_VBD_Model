#!/usr/bin/env python3
"""CP-L2-A6: fixed geometric arc-adjacent regional neutral-area diagnostic.

Region R_arc = { silicon points within 0.01 um of either canonical arc circle }.
Bounds come from the same Bernstein subdivision used for the whole domain: a
sub-triangle is assigned only when all three vertices are inside (rel > 0) or all
outside (rel < 0); straddlers are subdivided to AREA_TOL and counted in the
upper bound only. Complement bounds are derived from the certified whole-domain
enclosure by set subtraction.
"""
from __future__ import annotations
import json
import time
import numpy as np

import p2_audit as A
import p2_area_bounds as AB
from p2_bernstein import B
from p2_area_bounds import sub_coeff, AREA_TOL, MAX_DEPTH

M = A.M
R_ARC = abs(M.A_BOT - M.A_ARC_END)
TOL_R = 0.01


def arc_dist(pts):
    cx = pts[:, 0]; cy = pts[:, 1]
    d1 = np.abs(np.hypot(cx - M.TIP_START, cy - M.A_ARC_END) - R_ARC)
    d2 = np.abs(np.hypot(cx - M.TIP_START, cy + M.A_ARC_END) - R_ARC)
    return np.minimum(d1, d2)


def mid2(P, Q):
    return ((P[0]+Q[0])/2, (P[1]+Q[1])/2)


def mid3(P, Q):
    return ((P[0]+Q[0])/2, (P[1]+Q[1])/2, (P[2]+Q[2])/2)


def recurse(pc, bar, ph, area, depth, acc, inside_region):
    d = arc_dist(np.asarray(ph))
    if inside_region:
        rel_in = bool((d <= TOL_R).all())
        rel_out = bool((d > TOL_R).all())
    else:
        rel_in = bool((d > TOL_R).all())
        rel_out = bool((d <= TOL_R).all())
    if rel_out:
        return
    sc = sub_coeff(pc, bar[0], bar[1], bar[2])
    if rel_in:
        if all(c <= 0.0 for c in sc):
            acc[0] += area; return
        if all(c > 0.0 for c in sc):
            return
    if area < AB.AREA_TOL or depth >= MAX_DEPTH:
        acc[1] += area; return
    m01b = mid3(bar[0], bar[1]); m12b = mid3(bar[1], bar[2]); m20b = mid3(bar[2], bar[0])
    m01 = mid2(ph[0], ph[1]); m12 = mid2(ph[1], ph[2]); m20 = mid2(ph[2], ph[0])
    kids = (((bar[0], m01b, m20b), (ph[0], m01, m20)),
            ((m01b, bar[1], m12b), (m01, ph[1], m12)),
            ((m20b, m12b, bar[2]), (m20, m12, ph[2])),
            ((m01b, m12b, m20b), (m01, m12, m20)))
    for kb, kp in kids:
        recurse(pc, kb, kp, area/4.0, depth+1, acc, inside_region)


def region_bounds(nodes, tris, mat, u, inside_region):
    area = np.abs(A.element_geom(nodes, tris)[2])
    P = nodes[tris][:, :3]
    cx = P[:, :, 0].mean(1); cy = P[:, :, 1].mean(1)
    charged = (mat == M.SI) & (cx > A.BODY_BOTTOM + 1e-12)
    dmin = np.minimum(arc_dist(P[:, 0]), np.minimum(arc_dist(P[:, 1]), arc_dist(P[:, 2])))
    dmax = np.maximum(arc_dist(P[:, 0]), np.maximum(arc_dist(P[:, 1]), arc_dist(P[:, 2])))
    if inside_region:
        cand = charged & (dmin <= TOL_R)
    else:
        cand = charged & (dmax > TOL_R)
    lo = 0.0; un = 0.0
    bar0 = ((1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0))
    for e in np.flatnonzero(cand):
        c = B @ u[tris[e]] - AB.EPS
        acc = [0.0, 0.0]
        recurse(c, bar0, (tuple(P[e, 0]), tuple(P[e, 1]), tuple(P[e, 2])),
                area[e], 0, acc, inside_region)
        lo += acc[0]; un += acc[1]
    return lo, lo + un, int(cand.sum())


def main():
    levels = [("L0", "../p2_3_spatial_S0/CP-M2-C_L0_polished.npz",
               json.load(open("../p2_3_spatial_S0/A_S0-L0.json"))),
              ("L1", "../p2_3_spatial_S0/CP-L1_L1_polished.npz",
               json.load(open("../p2_3_spatial_S0/A_S0-L1.json"))),
              ("L2", "../p2_3_spatial_S0/CP-L2_L2_polished.npz",
               json.load(open("../p2_3_spatial_S0/A_S0-L2.json")))]
    AB.AREA_TOL = 1e-13
    out = dict(region="R_arc: |dist to canonical arc| <= 0.01 um",
               radius=R_ARC, tol=TOL_R, AREA_TOL=AB.AREA_TOL, levels={})
    for tag, path, whole in levels:
        z = np.load(path)
        t0 = time.time()
        alo, aup, nc = region_bounds(z["nodes"], z["tris"], z["mat"], z["u"], True)
        dt = time.time() - t0
        wl, wu = whole["lower"], whole["upper"]
        clo, cup = wl - aup, wu - alo
        out["levels"][tag] = dict(arc_lower=float(alo), arc_upper=float(aup),
                                  arc_width=float(aup-alo), candidates=nc,
                                  seconds=round(dt, 1),
                                  whole=[wl, wu],
                                  comp_lower=float(clo), comp_upper=float(cup),
                                  comp_note="derived by set subtraction, conservative")
        print(f"[{tag}] arc [{alo:.15f}, {aup:.15f}] width {aup-alo:.3e} "
              f"cand {nc} ({dt:.0f}s) | comp [{clo:.12f}, {cup:.12f}]", flush=True)
    json.dump(out, open("../p2_3_spatial_S0/CP-L2-A6_arc_region.json", "w"), indent=2)
    print("wrote CP-L2-A6_arc_region.json")


if __name__ == "__main__":
    main()
