#!/usr/bin/env python3
"""Read-only P1/P2 equivalence audit.

P1 archives store element blocks per material: Si, SiO2, Gate, Shield.
Recompute P1's own element_material() and report:
  - geometric silicon area before the METAL conversions
  - Si area after conversion (what P1 actually charges)
  - body-region area removed from charge support
  - integrated donor charge implied by each

No frozen result is modified.
"""
from __future__ import annotations
import sys
import numpy as np

sys.path.insert(0, "..")
import conforming_fem as CF

BLOCKS = {"Si": CF.SI, "SiO2": CF.OX, "Gate": CF.METAL, "Shield": CF.METAL}


def load_blocks(fn):
    d = np.load(fn)
    nodes = d["nodes"]
    tris, mat = [], []
    for k, v in BLOCKS.items():
        if k in d.files:
            arr = np.asarray(d[k], dtype=int)
            tris.append(arr)
            mat.append(np.full(len(arr), v))
    return nodes, np.vstack(tris), np.concatenate(mat)


def areas(nodes, tris):
    p = nodes[tris]
    det = ((p[:, 1, 0]-p[:, 0, 0])*(p[:, 2, 1]-p[:, 0, 1]) -
           (p[:, 2, 0]-p[:, 0, 0])*(p[:, 1, 1]-p[:, 0, 1]))
    return 0.5*np.abs(det)


def main():
    for fn in ("../r3_mesh.npz", "../r5_mesh.npz"):
        try:
            nodes, tris, block_mat = load_blocks(fn)
        except FileNotFoundError:
            print(f"{fn}: missing")
            continue
        if tris.shape[1] != 3:
            print(f"{fn}: {tris.shape[1]}-node elements, skipped")
            continue
        a = areas(nodes, tris)
        cx = nodes[tris, 0].mean(axis=1)
        cy = nodes[tris, 1].mean(axis=1)
        half = CF.half_width(cx, 0.10)

        raw_si = np.abs(cy) < half
        mat = CF.element_material(nodes, tris, 0.10)
        post_si = mat == CF.SI
        body = (cx <= CF.BODY_BOTTOM) & (np.abs(cy) < half)
        gate_metal = (cx >= CF.GATE_TOP) & (cx <= CF.GATE_BOTTOM) & (np.abs(cy) >= half + CF.OXIDE_T)
        shield_metal = (cx >= CF.SHIELD_TOP) & (cx <= CF.SHIELD_BOTTOM) & (np.abs(cy) >= half + CF.OXIDE_T)

        print(f"\n{fn}: elements={len(tris)} total area={a.sum():.6f}")
        print(f"  geometric Si (|y|<a)      = {a[raw_si].sum():.9f} um^2")
        print(f"  Si after METAL conversion = {a[post_si].sum():.9f} um^2")
        print(f"  body block (x<=0.54)      = {a[body].sum():.9f} um^2")
        print(f"  gate metal block          = {a[gate_metal].sum():.9f} um^2")
        print(f"  shield metal block        = {a[shield_metal].sum():.9f} um^2")
        print(f"  removed from Si           = {a[raw_si].sum() - a[post_si].sum():.9f} um^2")
        print(f"  donor charge (post)       = {CF.RHO_SI*a[post_si].sum():.6f}")
        print(f"  donor charge (geometric)  = {CF.RHO_SI*a[raw_si].sum():.6f}")

    print("\n[reference] P1 archived si_area      = 1.954558180854431")
    print("[reference] P2 r3 Si area (measured) = 1.954550")


if __name__ == "__main__":
    main()
