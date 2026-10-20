#!/usr/bin/env python3
"""Distribution of opposite-sign pairs on r3: how many are genuinely above the
near-zero band, and where they sit."""
from __future__ import annotations
import json
import numpy as np
import p2_interface as PI

dref = json.load(open("p2_nearzero.json"))["reference"]["d_ref"]
z = np.load("p2_r3_solution.npz")
recs = PI.interface_flux(z["nodes"], z["tris"], z["mat"], z["phi"])
ds = np.array([q["d_si"] for q in recs])
do = np.array([q["d_ox"] for q in recs])
dscale = (np.abs(ds) + np.abs(do)) / 2
opp = np.sign(ds) != np.sign(do)

print(f"D_ref={dref:.6f}  opposite-sign pairs={opp.sum()}")
print(f"  dscale min={dscale[opp].min():.3e} med={np.median(dscale[opp]):.3e} "
      f"p95={np.percentile(dscale[opp],95):.3e} max={dscale[opp].max():.3e}")
for th in (1e-3, 1e-2, 1e-1):
    n = int((dscale[opp] > th*dref).sum())
    print(f"  above {th:g}*D_ref ({th*dref:.3f}): {n}")
above = opp & (dscale > 1e-3*dref)
xs = np.array([q["x"] for q in recs])[above]
ys = np.array([q["y"] for q in recs])[above]
dsi_a = ds[above]; dox_a = do[above]
print(f"  above-threshold set n={above.sum()}")
if above.any():
    print(f"    x {xs.min():.3f}..{xs.max():.3f}  y {ys.min():.3f}..{ys.max():.3f}")
    print(f"    |d_si| min/med/max {np.abs(dsi_a).min():.3e}/{np.median(np.abs(dsi_a)):.3e}/{np.abs(dsi_a).max():.3e}")
    print(f"    |d_ox| min/med/max {np.abs(dox_a).min():.3e}/{np.median(np.abs(dox_a)):.3e}/{np.abs(dox_a).max():.3e}")
    print(f"    rel mismatch med={np.median(np.abs(dsi_a-dox_a)/(np.abs(dsi_a)+np.abs(dox_a))):.4f}")
