# P2.2 — Bernstein constraint layer, manufactured active test, and P1/P2 equivalence audit

**Date:** 19 October 2026
**Status:** P2.2 unit tests pass · manufactured active-obstacle test passes · inactive trench KKT audit passes · **P1/P2 physical equivalence FAILS** · P2.3 blocked · P2.1 closed · P1 frozen · Stage 3 suspended

## 1. Bernstein layer — unit tests

Transform: `b_v = u_v`, `b_ij = 2u_ij − (u_i+u_j)/2`, node order v0,v1,v2,m01,m12,m20.

| Test | Result |
|---|---:|
| Transform exactness, 200 random quadratic triangles | 7.105e-15 |
| Min interior value for random `b >= 0` | +1.5e-02 |
| Nodal non-negativity counterexample, `b = [1,1,0,-1,-0.5,-0.5]` | min −0.249854 (−1/4) |
| Sufficient-not-necessary, `b = [1,1,1,-1,1,-1]` | min +2.3e-09, function (2L1−1)² |

Constraint matrix on r3 silicon: 333,660 rows × 152,489 columns; 12 structural non-zeros per element = 667,320.

## 2. Inactive physical-trench census

| | r3 | r4 |
|---|---:|---:|
| Negative Bernstein coefficients | 0 | 0 |
| min b | +0.000e+00 | +0.000e+00 |
| min u on free Si nodes | +0.07973584 V | +0.06657532 V |
| R_free | 9.94e-15 | 1.30e-14 |

λ = 0 satisfies stationarity and complementarity trivially; feasibility holds. The constraint is inactive at both resolutions.

## 3. Manufactured active-obstacle test (OSQP 1.1.3)

Two-material patch: oxide left, silicon right; Dirichlet u=1 at x=0, u=0 at x=1; negative silicon source S=60; Bernstein constraints on silicon only.

| Check | Result |
|---|---:|
| Unconstrained min b | −7.652244e-02 |
| Negative coefficients / affected elements | 320 / 64 of 64 |
| QP status / iterations | solved / 150 |
| Primal min slack | −1.735e-18 |
| Active constraints | 256 / 384 |
| Stationarity, ‖Px+q+Aᵀy‖/‖r‖ | 3.339e-14 |
| λ = −y, min / max / negative count | −0.000e+00 / 7.813e-02 / 0 |
| Complementarity max |λ·slack| | 4.948e-20 |
| Pure-oxide free rows max|K u| | 1.648e-14 (119 nodes) |
| Dense interior sampling, 64 Si elements | −8.919e-25 |

This is the required nonempty-active-set validation: multipliers are nonnegative, stationarity holds to 3e-14, complementarity to 5e-20, oxide equilibrium to 2e-14, and sampled interior values are nonnegative.

## 4. P1/P2 equivalence audit — the two problems are not the same

Recomputed P1's own `element_material()` and `area[mat==SI]` on the archived P1 meshes (read-only):

| Quantity | P1 r3 | P1 r5 | P2 r3 |
|---|---:|---:|---:|
| Geometric Si (|y| < a) | 1.964870600 | 1.964706279 | — |
| Si after METAL conversion (charged) | **1.585441860** | **1.585451240** | **1.954550** |
| Body block x ≤ 0.54 | 0.379428740 | 0.379255038 | (charged Si) |
| Removed from Si | 0.379428740 | 0.379255038 | 0 |
| Donor charge, post-conversion | **1979.525384** | **1979.537096** | **2440.380108** |

**P2 carries 0.369108 µm² more charged silicon and 460.85 more donor-charge units than P1 — 23.3 % more.**

P1 `assemble` gives METAL the silicon permittivity (`eps = EPS_OX if OX else EPS_SI`) but skips the donor load (`if mat[e] != SI: continue`). P2 keeps the body slab as charged silicon and pins its nodes. These are different volume charge supports and different effective boundary implementations.

### Archived P1 normalization is inconsistent

The frozen record reports `silicon_area_um2 = 1.954558180854431`. P1's own post-conversion charged area on the archived r3 mesh is **1.585442**. The archived number matches the *pre-conversion* geometric area (1.964871) to 0.5 % and P2's total (1.954550) to 8 ppm. Therefore the archived `active_silicon_area_um2 / silicon_area_um2` fraction uses the wrong denominator:

- reported: 0.020330952 / 1.954558 = 1.04 %
- corrected: 0.020330952 / 1.585442 = **1.28 %**

This is a reporting correction only; the frozen solution is untouched.

## 5. Consequence

P1 and P2 do not solve the same discrete problem, so their activation states cannot be compared as a discretization study. P2's inactive margin (+0.0797 V at r3) and P1's active overshoot (−0.2125 V) differ by ~0.29 V, but the charge-support difference runs the wrong way to explain it (P2 has *more* charge yet *less* overshoot). The remaining candidate is the element formulation itself (P1 linear with vertex loads vs P2 quadratic with midside loads), plus any P1 complementarity-specific treatment.

## 6. Next steps

1. Build a **comparison variant** of P2 with the body slab uncharged (matching P1's charge support), clearly labelled, and compare unconstrained fields at fixed physical probes. Do not alter the frozen P1 baseline.
2. Re-run the frozen P1 linear (unconstrained) solve and save its field so the comparison is pointwise, not just the recorded scalar overshoot.
3. Only after the variants agree on unconstrained electrostatics, revisit activation and the neutral set.
4. No P2 trench obstacle solve and no further P2 refinement until then.
