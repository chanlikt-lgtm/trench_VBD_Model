# P2.3 — r3 → r4 spatial convergence under a common zero-set definition

**Date:** 19 October 2026
**Status:** P2.3 converging but not yet certified; neutral area grows +2.13 % between r3 and r4 in *both* formulations; interface mismatch −25 %; symmetry improves 3.8×; P1 frozen; Stage 3 suspended.

## 1. Common zero-set definition

Neutral region = `{u ≤ ε}` with u = V_D − φ, sampled on an identical 231-point barycentric grid per charged silicon element, each field using its own interpolation (linear for P1, quadratic for P2).

| ε | r3 P1 | r3 P2 | r3 P2/P1 | r4 P1 | r4 P2 | r4 P2/P1 |
|---|---:|---:|---:|---:|---:|---:|
| 1e-8 | 0.021862911 | 0.021584154 | 0.9872 (−1.28 %) | 0.022327871 | 0.022044317 | 0.9873 (−1.27 %) |
| 1e-6 | 0.021862911 | 0.021597463 | 0.9879 | 0.022327871 | 0.022105673 | 0.9900 |
| 1e-4 | 0.022157553 | 0.022583365 | 1.0192 | 0.022593534 | 0.023021518 | 1.0189 |

Both formulations grow by **+2.13 %** from r3 to r4; the P1/P2 ratio is stable at **−1.27 %**. The area is not yet constant, so the continuum limit is not certified; the stable ratio indicates a common discretization trend rather than a formulation error.

## 2. KKT residuals

| | r3 | r4 |
|---|---:|---:|
| Free DOFs | 123,752 | 240,606 |
| Constraints | 300,222 | 583,674 |
| Independent active rows | 646 | 1,219 |
| Stationarity | 6.687e-11 | 1.482e-10 |
| Min slack | −1.871e-16 | −2.336e-16 |
| λ min | 0.000e+00 | 0.000e+00 |
| Complementarity | 8.581e-18 | 5.868e-18 |
| Overshoot | 0 | 0 |

## 3. Neutral topology

| | r3 | r4 |
|---|---|---|
| Fully neutral elements | 291 | 571 |
| Element-level area | 0.020917675 | 0.021562124 |
| Components | 1 | 1 |
| x extent | 2.3320–2.4200 | 2.3328–2.4200 |
| y half-width | 0.2009 | 0.2150 |
| Drain-connected | yes | yes |

## 4. Front intersections (common criterion, ε = 1e-8)

| x | r3 P1 | r3 P2 | r4 P1 | r4 P2 |
|---|---|---|---|---|
| 2.35 | (−0.062, 0.084) | (−0.062, 0.054)+(0.072, 0.074) | (−0.068, 0.070) | (−0.068, 0.070) |
| 2.38 | (−0.134, 0.130) | (−0.134, 0.130) | (−0.134, 0.132) | (−0.134, 0.132) |
| 2.39 | (−0.148, 0.156) | (−0.148, 0.156) | (−0.154, 0.154) | (−0.154, 0.154) |
| 2.40 | (−0.176, 0.178) | (−0.176, 0.178) | (−0.172, 0.178) | (−0.172, 0.168) |
| 2.41 | (−0.194, 0.194) | (−0.194, 0.194) | (−0.188, 0.188)+(0.204, 0.204) | same as P1 |

The r3 discrepancy at x = 2.35 is resolved at r4; P1 and P2 now agree at four of five cuts. The sliver at x = 2.41 appears in both formulations at r4 and needs one more resolution.

## 5. Interface flux after activation

| Region | r3 | r4 | change |
|---|---:|---:|---:|
| Wall | 0.015503 | 0.011660 | −24.8 % |
| Arc proper | 0.014669 | 0.010883 | −25.8 % |
| Cap + corner | 0.032182 | 0.024996 | −22.3 % |

## 6. Symmetry

Maximum |u(x,y) − u(x,−y)|: **9.867e-04 V** at r3 → **2.589e-04 V** at r4 (3.8× improvement).

## 7. Assessment

- Converging: interface mismatch −25 %, symmetry 3.8× better, P1/P2 front agreement improved, stable P1/P2 area ratio.
- Not certified: neutral area +2.13 % between r3 and r4 in both formulations; the front moves by ≤ 0.006 µm; one more resolution is required to extrapolate the limit and to interpret the x = 2.41 sliver.
- Electric-field comparison at fixed probes has not yet been performed.

P2.3 remains open.

## 8. Three-resolution comparison (r3, r4, intermediate h=0.0085)

Intermediate solve: 404,129 DOFs, 328,393 free, 797,268 constraints, 132,878 charged elements. OSQP solved in 15,225 iterations (obj −9.819705e+04); 1,666 independent active rows.

| KKT | r3 | r4 | intermediate |
|---|---:|---:|---:|
| Stationarity | 6.687e-11 | 1.482e-10 | 2.219e-10 |
| Min slack | −1.871e-16 | −2.336e-16 | −4.298e-16 |
| λ min | 0 | 0 | 0 |
| Complementarity | 8.581e-18 | 5.868e-18 | 5.163e-18 |
| Overshoot | 0 | 0 | 0 |

### Common-criterion neutral area (ε = 1e-8 V)

| Mesh | P2 area (µm²) | Increment | Increment ratio |
|---|---:|---:|---:|
| r3 | 0.021584154 | — | — |
| r4 | 0.022044317 | +0.000460163 | — |
| intermediate | 0.022200727 | +0.000156410 | **0.340** |

The increment falls by a factor **2.94**. Element-level fully neutral area behaves identically: 0.020917675 → 0.021562124 → 0.021781461, increments +0.000644449 and +0.000219337 (ratio 0.340).

### Topology

| | r3 | r4 | intermediate |
|---|---:|---:|---:|
| Components | 1 | 1 | 1 |
| x extent | 2.3320–2.4200 | 2.3328–2.4200 | 2.3343–2.4200 |
| y half-width | 0.2009 | 0.2150 | 0.2114 |

### Front intersections (ε = 1e-8 V)

| x | r3 | r4 | intermediate |
|---|---|---|---|
| 2.35 | (−0.062, 0.054)+(0.072, 0.074) | (−0.068, 0.070) | (−0.074, 0.068) |
| 2.38 | (−0.134, 0.130) | (−0.134, 0.132) | (−0.130, 0.132) |
| 2.39 | (−0.148, 0.156) | (−0.154, 0.154) | (−0.156, 0.158) |
| 2.40 | (−0.176, 0.178) | (−0.172, 0.168) | (−0.180, 0.174) |
| 2.41 | (−0.194, 0.194) | (−0.188, 0.188)+(0.204, 0.204) | (−0.196,−0.192)+(−0.188,0.190)+(0.194,0.198) |

### Symmetry and interface

| | r3 | r4 | intermediate |
|---|---:|---:|---:|
| Max asymmetry (V) | 9.867e-04 | 2.589e-04 | 2.742e-04 |
| Wall R_abs | 0.015503 | 0.011660 | 0.010813 |
| Arc R_abs | 0.014669 | 0.010883 | 0.009235 |
| Cap R_abs | 0.032182 | 0.024996 | 0.023428 |

### Assessment

- Neutral-area increment falls 2.94×, consistent with convergence; the residual increment is 0.000156 µm² (0.70 % of the area). No asymptotic order is certified.
- Front positions at x = 2.38, 2.39, 2.40 are stable within 0.004–0.006 µm.
- The upstream edge at x = 2.35 continues to move outward (−0.062 → −0.068 → −0.074) and is not stabilized.
- The x = 2.41 sliver persists and fragments: one interval at r3, two at r4, three at the intermediate resolution.
- Symmetry is stable at ~2.6–2.7e-4 V; all interface mismatches continue to decrease.

P2.3 is converging in area and interface, but the upstream tip and the x = 2.41 feature are not certified. Targeted local refinement of those features is indicated before any full r5 solve.

## 9. Post-audit corrections (F1-F3)

### F1 — P1 comparator residual (HIGH, fixed)

The active-set residual was `K u - rhs` with `rhs = -f - K u_D`, double-counting the Dirichlet terms. Corrected to `r = K u + f` on free rows. Re-run on the existing meshes:

| | r3 | r4 |
|---|---:|---:|
| active nodes | 152 | 297 |
| min u | 0 | 0 |
| min residual (active) | +1.503e-02 | +2.331e-03 |
| max abs residual (inactive) | 1.701e-12 | 1.632e-12 |
| complementarity | 0 | 0 |

The active sets and areas were unchanged; the KKT audit is now explicit and passes.

### F2 — controlled-error neutral area (equal-area quadrature)

The 231-point lattice over-weights edges and vertices. Replaced with equal-area sub-triangle centroids at K = 8, 16, 32.

| Mesh | K=8 | K=16 | K=32 | K-trend |
|---|---:|---:|---:|---:|
| r3 | 0.020917675 | 0.020917675 | 0.020917675 | exact element union |
| r4 | 0.021562124 | 0.021569662 | 0.021570700 | +7.5e-6, +1.0e-6 |
| intermediate | 0.021781461 | 0.021781624 | 0.021788063 | +1.6e-7, +6.4e-6 |

Increments (K=32): r3→r4 = 0.000653025; r4→mid = 0.000217363; **ratio 0.333**. The 3× reduction survives the estimator correction. The lattice values were biased high by ~3 %, but the increment ratio was nearly unchanged (0.340 → 0.333) because the bias largely cancels between meshes.

P1 exact clipping: r3 0.021196729; r4 0.021845237. P1/P2 ratios: **0.98684** and **0.98743** (−1.32 %, −1.26 %), stable and close to the previous −1.27 %.

### F3 — analytic front cuts

Vertical cuts were replaced by exact element/line intersection with bisection of u(t)=eps. Main intervals:

| x | r3 | r4 | intermediate |
|---|---|---|---|
| 2.35 | (−0.0624, 0.0542) | (−0.0697, 0.0714) | (−0.0754, 0.0686) |
| 2.38 | (−0.1355, 0.1313) | (−0.1355, 0.1330) | (−0.1314, 0.1331) |
| 2.39 | (−0.1490, 0.1576) | (−0.1544, 0.1550) | (−0.1564, 0.1593) |
| 2.40 | (−0.1763, 0.1795) | (−0.1729, 0.1696) | (−0.1816, 0.1746) |
| 2.41 | (−0.1945, 0.1958) | (−0.1892, 0.1899)+(0.2037, 0.2050) | (−0.1895, 0.1903)+(0.1922, 0.1981) |

The x = 2.35 lower edge still drifts outward; the x = 2.41 sliver persists. The extractor emits small boundary artefacts that should be filtered before the final record.

## 10. F2 — provable Bernstein area bounds (closed)

Each charged element's `u - eps` is expressed in Bernstein form. Convex-hull classification (all coefficients <= 0 -> neutral; all > 0 -> outside; mixed -> subdivide) gives

```
LOWER <= A_true <= UPPER,   UPPER - LOWER = unresolved area
```

Result at eps = 1e-8 V, subdivision tolerance 1e-11 um^2:

| Mesh | LOWER | UPPER | gap | mixed elements |
|---|---:|---:|---:|---:|
| r3 | 0.020926837673 | 0.020930251232 | 3.414e-06 | 189 |
| r4 | 0.021574089791 | 0.021578224392 | 4.135e-06 | 261 |
| intermediate | 0.021789286073 | 0.021793581617 | 4.296e-06 | 310 |

**Correction.** An earlier version of this section reported lower = upper with gap 0 (0.021318701188 / 0.021851690577 / 0.022032403652). That was a coordinate-frame bug: during recursion the sub-triangle Bernstein coefficients were evaluated at parent-frame nodes. After the fix the bounds are genuine intervals with a positive gap. The manufactured linear test now returns [0.375, 0.375 + 1.907e-06], containing the exact 0.375.

Increment intervals: Delta(r3->r4) in [0.000643839, 0.000651387]; Delta(r4->mid) in [0.000211062, 0.000219492]; ratio in [0.324, 0.341], midpoint **0.332**. The ~3x reduction is robust within these bounds.

P1 exact clipping: 0.021196729 (r3), 0.021845237 (r4). P2/P1 = **0.9873–0.9874** (r3) and **0.9876–0.9878** (r4): stable near 0.987, i.e. −1.2 to −1.3 %.

The gaps (~4e-6) are 0.6–1.9 % of the increments. A tighter AREA_TOL (1e-13) would reduce them below 1e-6 at roughly 10x the runtime; not required for the convergence comparison.

### Remaining audit items

F4 (maps use vertex-only contours: qualitative), F5 (patch test must call `interface_flux`), F7 (hard KKT/status gates), F9 (edge-connected components), F11 (streamline vectors), F6/F8 (quarantine legacy assembler and `conforming_fem.solve_complementarity`). The frozen P1 archive was not regenerated through the legacy routine; its provenance still needs an independent audit.
