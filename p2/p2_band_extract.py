#!/usr/bin/env python3
"""Extract the frozen refinement band from the certified L1b solution.

Rule (documented, reproducible):
  * contour the L1b P2 field at u = 1e-8 on the silicon sub-triangulation;
  * split the contour into connected polylines;
  * DROP any component whose arc length is below L_MIN = 0.02 um
    (numerical slivers at region edges);
  * keep every component at or above L_MIN, in decreasing length order;
  * resample the kept polylines at DS = 0.002 um.

The kept components define the band for all levels; the band is frozen here and
must not be re-derived from later solutions.
"""
from __future__ import annotations
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.tri as mtri

import p2_audit as A

M = A.M
EPS = 1e-8
L_MIN = 0.02
DS = 0.002


def resample(poly, ds):
    d = np.linalg.norm(np.diff(poly, axis=0), axis=1)
    s = np.concatenate([[0.0], np.cumsum(d)])
    if s[-1] <= 0:
        return poly[:1]
    n = max(int(s[-1]/ds) + 1, 2)
    si = np.linspace(0.0, s[-1], n)
    return np.column_stack([np.interp(si, s, poly[:, 0]),
                            np.interp(si, s, poly[:, 1])])


def main():
    z = np.load("p2_L1b_scaled.npz")
    nodes, tris, mat, u = z["nodes"], z["tris"], z["mat"], z["u"]
    si = mat == M.SI
    t3 = tris[si][:, :3]
    tri = mtri.Triangulation(nodes[:, 0], nodes[:, 1], t3)
    fig, ax = plt.subplots()
    cs = ax.tricontour(tri, u, levels=[EPS])
    segs = [np.asarray(s) for s in cs.allsegs[0] if len(s) >= 2]
    plt.close(fig)

    comps = []
    for s in segs:
        L = float(np.linalg.norm(np.diff(s, axis=0), axis=1).sum())
        comps.append((L, s))
    comps.sort(key=lambda t: -t[0])
    kept = [(L, s) for L, s in comps if L >= L_MIN]
    dropped = [(L, s) for L, s in comps if L < L_MIN]

    print(f"components: {len(comps)}  kept {len(kept)}  dropped {len(dropped)}")
    for i, (L, s) in enumerate(kept):
        print(f"  keep {i}: length {L:.4f} um  points {len(s)}  "
              f"x[{s[:,0].min():.3f},{s[:,0].max():.3f}] "
              f"y[{s[:,1].min():.3f},{s[:,1].max():.3f}]")
    for L, s in dropped:
        print(f"  drop  : length {L:.4f} um  points {len(s)}")

    pts = np.vstack([resample(s, DS) for _, s in kept])
    np.savez_compressed("p2_band_contour.npz", pts=pts,
                        lengths=np.array([L for L, _ in kept]),
                        dropped=np.array([L for L, _ in dropped]))
    json.dump(dict(n_components=len(comps), kept=[float(L) for L, _ in kept],
                   dropped=[float(L) for L, _ in dropped],
                   total_kept_length=float(sum(L for L, _ in kept)),
                   l_min=L_MIN, ds=DS, eps=EPS,
                   n_points=int(len(pts))),
              open("p2_band_contour.json", "w"), indent=2)
    print(f"\nband points {len(pts)}  total length {sum(L for L,_ in kept):.3f} um")
    print("wrote p2_band_contour.npz / .json")


if __name__ == "__main__":
    main()
