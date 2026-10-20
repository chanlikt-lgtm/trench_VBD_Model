# P2.3 — KKT polish and neutral-geometry audit

**Date:** 19 October 2026
**Status:** polish at machine precision; neutral set one drain-connected component of half-width ≈ 0.20 µm; interface flux re-audited after activation; P2.3 remains open pending one refinement; P1 frozen; Stage 3 suspended.

## 1. KKT polish

Active rows reduced from 2,229 to **646 independent** rows by pivoted QR of the Gram matrix (the record's warning about dependent active rows was correct).

| Quantity | Before polish | After polish |
|---|---:|---:|
| Stationarity (relative) | 1.933e-05 | **6.687e-11** |
| Min slack | −3.459e-10 | −1.871e-16 |
| λ min / max | −2.804e-17 / 3.074e-04 | **0.000e+00** / 1.258e-01 |
| Complementarity | 1.383e-12 | **8.581e-18** |
| Overshoot | +1.080e-11 | +0.000e+00 |
| Pure-oxide residual | — | 1.805e-12 |

## 2. Neutral-set classes

| Class | Elements | Area (µm²) |
|---|---:|---:|
| Charged silicon | 50,037 | 1.575120976 |
| Fully neutral (all six nodes \|u\| < 1e-8) | 291 | **0.020917675** |
| Any zero node | 480 | 0.035064617 |
| Drain-adjacent | 171 | 0.013474421 |
| Fully neutral ∧ drain-adjacent | 59 | 0.004657445 |
| Fully neutral ∧ not drain-adjacent | 232 | 0.016260230 |

## 3. Components and drain connectivity

**One** connected component: 291 elements, area 0.020917675 µm²,
x [2.3320, 2.4200], **y [−0.2009, +0.2009]**, touching the drain.

The earlier ±0.575 transverse bound came from the any-node class (drain-adjacent elements with a fixed u = 0 node). The fully neutral region is a genuine central plug of half-width ≈ 0.20 µm, not a full-width boundary artifact.

Drain-plane silicon nodes: 167, all with |u| < 1e-8 (fraction 1.0).

## 4. Front intersections

| x (µm) | Neutral interval in y (µm) |
|---|---|
| 2.35 | (−0.062, +0.054) and a small second interval (0.072, 0.074) |
| 2.38 | (−0.134, +0.130) |
| 2.39 | (−0.148, +0.156) |
| 2.40 | (−0.176, +0.178) |
| 2.41 | (−0.194, +0.194) |

The plug widens downstream. Upper/lower asymmetry is ≈ 0.008 µm at 2.38–2.39; maximum potential asymmetry over the sampled region is 9.867e-04 V at (2.34, 0.06).

## 5. Interface flux after activation

| Region | n | R_abs | R_sgn | R_abs before obstacle (r3 linear) |
|---|---:|---:|---:|---:|
| Wall | 896 | 0.015503 | 0.008270 | 0.014307 |
| Arc proper | 90 | 0.014669 | 0.013878 | 0.014412 |
| Cap + corner | 158 | 0.032182 | 0.029529 | 0.042254 |
| All | 1,144 | 0.017264 | — | 0.017702 |

Activation slightly increases wall and arc mismatch and reduces the cap/corner mismatch.

## 6. Comparison with P1

| Quantity | P1 | P2 |
|---|---:|---:|
| Neutral/active area | 0.020330952 µm² (70 active nodes) | 0.020917675 µm² (291 fully neutral elements) |
| Upstream extent | front ≈ 2.35 µm | 2.332 µm |
| Half-width | — | 0.2009 µm |

Area difference +2.9 %, under different element-level definitions. P1 and P2 zero sets have not yet been compared with a single common physical definition.

## 7. Remaining gates

1. P2 zero-set comparison with P1 using one consistent physical criterion.
2. One P2 refinement (r4) with the obstacle to test neutral-area and front convergence.
3. Electric-field comparison at fixed probes.

P2.3 remains open until these pass.
