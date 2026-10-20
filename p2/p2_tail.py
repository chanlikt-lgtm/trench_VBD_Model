#!/usr/bin/env python3
"""Inspect the r4 wall tail: pairs whose relative mismatch is ~1.0 after the
near-zero exclusion, and count interface edges the matcher refused to pair."""
from __future__ import annotations
import json
import numpy as np
import p2_interface as PI

M = PI.M


def main():
    with open("p2_nearzero.json") as fh:
        dref = json.load(fh)["reference"]["d_ref"]

    z = np.load("p2_r4_solution.npz")
    nodes, tris, mat, phi = z["nodes"], z["tris"], z["mat"], z["phi"]
    recs = PI.interface_flux(nodes, tris, mat, phi)

    w = np.array([q["ds"] for q in recs])
    ds = np.array([q["d_si"] for q in recs])
    do = np.array([q["d_ox"] for q in recs])
    dscale = (np.abs(ds) + np.abs(do)) / 2
    rel = np.abs(ds - do) / np.maximum(2 * dscale, 1e-30)
    keep = dscale >= 1e-3 * dref

    print(f"samples={len(recs)} kept={keep.sum()} near_zero={(~keep).sum()}")
    for reg in ("wall", "arc"):
        m = np.array([q["region"] == reg for q in recs]) & keep
        bad = m & (rel > 0.99)
        print(f"{reg}: kept={m.sum()} rel>0.99: {bad.sum()} "
              f"({100*bad.sum()/max(m.sum(),1):.2f}%)")
        if bad.any():
            print(f"   d_si min/med/max {np.min(np.abs(ds[bad])):.3e} / "
                  f"{np.median(np.abs(ds[bad])):.3e} / {np.max(np.abs(ds[bad])):.3e}")
            print(f"   d_ox min/med/max {np.min(np.abs(do[bad])):.3e} / "
                  f"{np.median(np.abs(do[bad])):.3e} / {np.max(np.abs(do[bad])):.3e}")
            print(f"   exact zeros: si={int((ds[bad]==0).sum())} ox={int((do[bad]==0).sum())}")
            # coordinate spread of the bad set
            xs = np.array([q["x"] for q, b in zip(recs, bad) if b])
            ys = np.array([q["y"] for q, b in zip(recs, bad) if b])
            print(f"   x range {xs.min():.3f}..{xs.max():.3f}  y range {ys.min():.3f}..{ys.max():.3f}")

    # edge matcher accounting on both grids
    for tag, fn in (("r3", "p2_r3_solution.npz"), ("r4", "p2_r4_solution.npz")):
        zz = np.load(fn)
        tn, tt, tm = zz["nodes"], zz["tris"], zz["mat"]
        key = {}
        for e in range(len(tt)):
            t = tt[e]
            for (a, b, mm) in ((0, 1, 3), (1, 2, 4), (2, 0, 5)):
                k = (min(int(t[a]), int(t[b])), max(int(t[a]), int(t[b])))
                key.setdefault(k, []).append(e)
        n3 = sum(1 for v in key.values() if len(v) > 2)
        mixed = sum(1 for v in key.values()
                    if len(v) == 2 and {tm[v[0]], tm[v[1]]} == {M.SI, M.OX})
        one = sum(1 for v in key.values()
                  if len(v) == 1 and tm[v[0]] in (M.SI, M.OX))
        print(f"[{tag}] edges total={len(key)} nonmanifold(>2)={n3} "
              f"si-si or ox-ox pairs={sum(1 for v in key.values() if len(v)==2 and tm[v[0]]==tm[v[1]] and tm[v[0]] in (M.SI,M.OX))} "
              f"mixed={mixed}")


if __name__ == "__main__":
    main()
