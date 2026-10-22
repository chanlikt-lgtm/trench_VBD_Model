#!/usr/bin/env python3
"""L1 post-processing: bounded neutral area and root-based front cuts."""
from __future__ import annotations
import json
import numpy as np

import p2_area_bounds as AB
from p2_exact_metrics import cut_intervals
from p2_topology_clean import clean

M = AB.M


def main():
    z = np.load("p2_B1X_scaled.npz")
    nodes, tris, mat, u = z["nodes"], z["tris"], z["mat"], z["u"]
    lo, up, nch, nm = AB.bounds(nodes, tris, mat, u)
    print(f"[B1X area] charged={nch} mixed={nm}  lower={lo:.12f}  upper={up:.12f} "
          f"gap={up-lo:.3e}")

    out = {"area": dict(lower=lo, upper=up, gap=up-lo, charged=nch, mixed=nm),
           "cuts": {}}
    for x0 in (2.35, 2.38, 2.39, 2.40, 2.41):
        out["cuts"][str(x0)] = {}
        for eps in (1e-8, 1e-6, 1e-5):
            iv = cut_intervals(nodes, tris, mat, u, x0, eps, True)
            main, short = clean(iv)
            out["cuts"][str(x0)][str(eps)] = dict(main=main, short=len(short))
            print(f"  x={x0:.2f} eps={eps:g}: main={[(round(a,4),round(b,4)) for a,b in main]}"
                  f"  short={len(short)}")

    json.dump(out, open("p2_B1X_post.json", "w"), indent=2)
    print("wrote p2_B1X_post.json")


if __name__ == "__main__":
    main()
