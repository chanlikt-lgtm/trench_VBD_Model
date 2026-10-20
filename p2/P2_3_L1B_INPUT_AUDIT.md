# P2.3 — L1b Input Audit (mesh, realised sizes, charge mask, donor load)

**Record date:** 2026-10-20
**Status:** L1b mesh and inputs audited; QP solved; polishing in progress.
L1 remains the F7-certified reference. P2.3 open; P1 frozen; Stage 3 suspended.

## 1. Mesh

| Metric | L1b |
|---|---:|
| Nodes | 426,380 |
| Elements | 212,505 |
| Si elements / area | 158,667 / 1.954544960806 µm² |
| SiO2 / Gate / Shield elements | 31,728 / 6,448 / 15,662 |
| Minimum triangle area | 2.909e-06 |
| Maximum aspect ratio | 2.52 |
| Minimum interior angle | 36.09° |
| `H_LOCAL` requested | 0.0055 µm |
| Background nominal | 0.0085 µm, tip 0.0034 µm |

Windows unchanged from L1: A = `2.32 ≤ x ≤ 2.37, |y| ≤ 0.13`;
B = `2.395 ≤ x ≤ 2.420, |y| ≤ 0.24`; transition `Thickness = 0.012`.

## 2. Realised element sizes

Measured as the equilateral-equivalent side `h = sqrt(4A/√3)` from element
areas, not from Gmsh targets.

| Region | Intermediate | L1 | L1b | L1b / intermediate |
|---|---:|---:|---:|---:|
| A core (median h) | 0.00780 | 0.00399 | 0.00524 | 0.672 |
| B core (median h) | 0.00807 | 0.00399 | 0.00527 | 0.653 |

The nominal 0.0085 background is not realised near the front: the base
grading is already finer there (0.0078 at the windows, ~0.0044 near
x = 1.9). Contrast is therefore quoted against the intermediate mesh in the
same box, not against a far-field median.

Both L1 and L1b spread the window size outward through `Thickness = 0.012`,
so neither has a sharp window edge. L1b consequently tests **local element
size**, not an isolated change in transition sharpness; a difference in
neutral area or front position cannot be attributed uniquely to refinement
contrast.

## 3. Charge mask audit (L1 vs L1b)

`charged = (mat == SI) & (element centroid x > 0.54)`. Exact comparison by
Sutherland–Hodgman clipping of each Si element against `x ≥ 0.54`.

| Quantity | L1 | L1b |
|---|---:|---:|
| Exact Si area (µm²) | 1.954544960806 | 1.954544960806 |
| Exact charged area (µm²) | 1.575254229099 | 1.575254229099 |
| Centroid-classified charged area (µm²) | 1.575290810623 | 1.575280762537 |
| Charged-area overcount (µm²) | +0.000036582 | +0.000026533 |
| Relative overcount | +2.32e-05 | +1.68e-05 |
| Straddling elements at x = 0.54 | 220 | 213 |
| Straddle area (µm²) | 4.756e-03 | 4.629e-03 |
| Straddle above / below line (µm²) | 2.501e-03 / 2.255e-03 | 2.315e-03 / 2.314e-03 |

Observed charged-area change `-1.004809e-05 µm²` equals the change in
classification error (`+2.653e-05` minus `+3.658e-05`) with residual
`-2.2e-16`. The physical donor region and silicon geometry are therefore
identical, to machine precision, across the two meshes.

## 4. Interpretation (revised)

Sharing a classification rule does **not** make the comparison unbiased:

> Both meshes use the same centroid-classification rule, and their
> charged-area difference is fully explained by its mesh-dependent
> integration error. The physical donor region is unchanged, but the
> resulting discrete load difference is not guaranteed to cancel in the
> L1/L1b solution comparison.

The discrete charge source is not identical across meshes. In
`assemble_p2` the donor load is applied as a whole-element lump,

```
f[tris[e][3:]] += RHO * area[e] / 3.0
```

on the three P2 mid-side nodes of each charged element. Elements crossing
`x = 0.54` therefore contribute their full area regardless of how much of
them lies on the body side. L1 and L1b carry different source-integration
errors, these enter the assembled load vector, and in an obstacle problem a
small load perturbation can move an active-set boundary. The comparison is
between two discretisations of the same implemented formulation, not between
two identical discrete problems.

## 5. Terminology

The charge mask and the Bernstein constraint mask are conceptually distinct.
A centroid-charged triangle may extend into the body-side region even though
its centroid does not lie there. In the current code both masks are generated
from the same boolean `charged`, so they coincide; they are documented
separately here so that a future change to one is not silently inherited by
the other.

## 6. Deferred control experiment

If L1b continues to show unexpected mesh sensitivity, the appropriate
control is **exact partial-element donor-load integration** at `x = 0.54`:
integrating the quadratic P2 shape functions over the clipped charged
portion of each straddling element, rather than correcting only its scalar
charged area. That would separate donor-mask quadrature error from effects of
local mesh refinement. It is a separate controlled experiment and must not be
introduced as a silent change to the L1b formulation.

## 7. State at time of writing

L1b OSQP: `solved`, 14,625 iterations, objective `-1.014015e+05`,
primal residual `1.385e-08`, dual residual `1.663e-09`; checkpoint
`p2_L1b_osqp_checkpoint.npz` written and reload-verified
(`x(350190,)`, `y(862794,)`, `idx(350190,)`).

Free DOFs 350,190; constraints 862,794; charged elements 143,799;
body elements 14,868.

F7 certification, area bounds and front cuts are **pending** and are not
reported here.
