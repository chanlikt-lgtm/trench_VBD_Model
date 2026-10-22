#!/usr/bin/env python3
"""AREA_TOL = 1e-13 postprocessing audit for all four certified solutions.

Postprocessing only: no OSQP, no polish. The module-level tolerance in
p2_area_bounds is raised in place before each call, so the certified solution
archives are read and never modified.
"""
from __future__ import annotations
import json
import time
import numpy as np

import p2_area_bounds as AB

TOL = 1e-13
SOLUTIONS = [("B0", "p2_B0_scaled.npz"),
             ("B1", "p2_B1_scaled.npz"),
             ("B0X", "p2_B0X_scaled.npz"),
             ("B1X", "p2_B1X_scaled.npz")]


def main():
    AB.AREA_TOL = TOL
    print(f"AREA_TOL = {AB.AREA_TOL:g}  MAX_DEPTH = {AB.MAX_DEPTH}", flush=True)
    out = {}
    for tag, fn in SOLUTIONS:
        z = np.load(fn)
        t0 = time.time()
        lo, up, nch, nm = AB.bounds(z["nodes"], z["tris"], z["mat"], z["u"])
        dt = time.time() - t0
        out[tag] = dict(lower=float(lo), upper=float(up), width=float(up-lo),
                        charged=int(nch), mixed=int(nm), seconds=round(dt, 1))
        print(f"[{tag}] lower={lo:.15f} upper={up:.15f} width={up-lo:.3e} "
              f"mixed={nm} ({dt:.0f}s)", flush=True)
    json.dump(out, open("p2_area_a13.json", "w"), indent=2)
    print("wrote p2_area_a13.json")


if __name__ == "__main__":
    main()
