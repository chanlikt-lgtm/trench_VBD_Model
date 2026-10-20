# P2 quadrature defect — factor-2 stiffness error

**Date:** 19 October 2026
**Status:** P2 linear results invalidated; P2.1 closure suspended; P2.2 census invalidated; P1 correct and unchanged; P2.3 blocked; Stage 3 suspended.

## 1. Detection

Nested-space identity on the r3 mesh, before Dirichlet elimination:

| Check | Result |
|---|---:|
| ‖f_P1 − Tᵀ f_P2‖ / ‖f_P1‖ | 1.098e-16 |
| ‖K_P1 − Tᵀ K_P2 T‖ / ‖K_P1‖ | 5.000e-01 |
| ‖K_P1 − ½ Tᵀ K_P2 T‖ / ‖½ Tᵀ K_P2 T‖ | **2.359e-16** |
| diagonal ratio median | **2.000000000** (p05 = p95) |

The load transfer is exact. The operator is exactly twice the correct one.

## 2. Root cause

The P2 stiffness is assembled as

```python
ke += QW[q] * abs(det[m]) * eps[m] * (G @ G.T)
```

with `QW = [1/3, 1/3, 1/3]`, so `sum(QW) = 1.0`. The reference triangle has area **0.5**, so a rule that integrates 1 exactly must have weights summing to **0.5**. The assembled element integral of the constant 1 is therefore `sum(QW)·|det| = 2A` instead of `A`: every P2 stiffness entry is doubled.

The physical load `f = ρA/3` on midside nodes is correct, so every P2 linear solve to date solved

```
2 K φ = f
```

i.e. the charge contribution was halved relative to the Dirichlet data. That is exactly the signature observed: P1 overshoots, P2 did not.

The earlier patch tests could not catch this because both the stiffness and the test load used the same over-scaled quadrature, so the scaling cancelled.

## 3. Corrected solutions (K/2, same meshes, same BCs, same load)

| | r3 | r4 |
|---|---:|---:|
| R_free | 1.304e-14 | 1.678e-14 |
| phi range | [0, 34.490612] | [0, 34.490945] |
| max(phi) − V_D | **+0.216403 V** | **+0.216737 V** |
| min u free Si | −0.216403 V | −0.216737 V |
| location of min u | (2.35372, 0.00354) | (2.35177, −0.00039) |
| min Bernstein b | −2.174e-01 | −2.175e-01 |
| negative coefficients | 3,598 | 7,015 |
| affected elements | 678 / 55,610 | 1,272 / 107,973 |

P1 for comparison: +0.220543 V on the same r3 mesh, +0.232581 V on its own archive mesh, minimum at (2.35115, −0.00261).

The corrected P2 overshoot and the location of the minimum now agree with P1 to within the discretization difference. The Bernstein constraint is genuinely violated, so the obstacle is active and P2.3 becomes meaningful.

## 4. What is invalidated

- P2 linear solutions r3 and r4 (solved `2Kφ = f`).
- P2.1 interface audit numbers: R_abs 0.019182 → 0.015473, wall 0.011079 → 0.009343, arc 0.033922 → 0.026612. All were computed from the wrong solution. **P2.1 closure is suspended pending re-run.**
- P2.2 inactive-trench census: min free u +0.0797/+0.0666 V and "zero violations" are wrong. The corrected census shows 3,598/7,015 violations.
- P1/P2 activation margins and the probe comparison in `P1_P2_EQUIVALENCE_STUDY.md` (the P2 side is invalid; the P1 side stands).

## 5. What remains valid

- The Bernstein transform and its unit tests (algebraic, independent of K).
- The manufactured active-obstacle KKT validation (scale-invariant; the QP was self-consistent).
- The P1 reference and its assembly (it is the reference in the identity, and Tᵀ(½K₂)T = K₁ to 2.4e-16).
- The nested-space test itself, which found the defect.

## 6. Fix and re-run plan

1. Change the P2 quadrature to weights summing to 0.5 (`QW/2`, or multiply the weight by `|det|/2 = A`) in `p2_audit.assemble` and every copy (`p2_diag.py`, `p2_obstacle_test.py`, `p2_obstacle_osqp.py`, `p2_corrected_check.py`).
2. Re-verify `Tᵀ K₂ T = K₁` and `Tᵀ f₂ = f₁` to machine precision.
3. Re-run P2 linear solves r3/r4, the interface audit, and the Bernstein census.
4. Re-validate P2.1 on the corrected solutions; do not assume the previous convergence conclusion survives.
5. Then run the P2.3 obstacle QP, which now has a nonempty active set.
6. Update the P2.1 and P2.2 records; keep P1 frozen.

## 7. Re-validation after the fix

Correction applied: `QW = [1/6, 1/6, 1/6]` in `p2_audit.py`; every dependent script imports it.

### 7.1 Nested identities

| | r3 | r4 |
|---|---:|---:|
| ‖K1 − TᵀK2T‖/‖K1‖ | **2.359e-16** | **2.347e-16** |
| max abs difference | 3.553e-14 | 4.263e-14 |
| ‖f1 − Tᵀf2‖/‖f1‖ | 1.098e-16 | 1.095e-16 |

The formulation now matches to machine precision.

### 7.2 Manufactured patch tests (r3)

| Check | Result |
|---|---:|
| ∫N/A, three Si elements | [0, 0, 0, 1/3, 1/3, 1/3] |
| six-DOF graph | 1 component, 0 unanchored |
| K_FF zero diagonals / rows without off-diagonal / asymmetry | 0 / 0 / 0 |
| D1, phi = x²−y², zero source | 1.776e-14 (rel 3.97e-13) |
| D2, phi = x²+y², source −4ε | 2.156e-14 (rel 4.80e-13) |

### 7.3 Corrected linear solves

| | r3 | r4 |
|---|---:|---:|
| R_free | 1.304e-14 | 1.678e-14 |
| phi max | 34.490612 | 34.490945 |
| overshoot | **+0.216403 V** | **+0.216737 V** |
| min u free Si | −0.216403 V | −0.216737 V |
| min-u location | (2.353717, 0.003543) | (2.351772, −0.000390) |
| min Bernstein b | −2.173578e-01 | −2.174655e-01 |
| negative coefficients | 3,598 | 7,015 |
| affected elements | 678 / 55,610 | 1,272 / 107,973 |

P1 same-mesh: +0.220543 V at (2.35115, −0.00261). The corrected P2 overshoot and location now agree with P1.

### 7.4 Corrected interface audit (frozen reconstruction, unchanged)

| | r3 | r4 | ratio | p |
|---|---:|---:|---:|---:|
| R_abs all | 0.017702 | 0.014568 | 0.823 | 0.55 |
| R_abs wall | 0.014307 | 0.012101 | 0.846 | 0.50 |
| R_abs arc | 0.030967 | 0.024254 | 0.783 | 0.72 |

Both wall and arc decrease independently, as before, with revised magnitudes. The full classified audit (electrode-adjacent, near-zero, arc proper) has not yet been re-run on the corrected fields.

### 7.5 Superseded evidence preserved

Original solutions remain at `p2_r3_solution.npz`, `p2_r4_solution.npz`; corrected solutions are written separately as `p2_r3_fixed.npz`, `p2_r4_fixed.npz`.

### 7.6 Remaining re-validation

1. Re-run the classified interface audit (electrode-adjacent / near-zero / arc proper) on the corrected solutions.
2. Repeat the same-mesh P1/P2 probe comparison with the corrected P2 variant.
3. Only then consider P2.3.

## 8. Final re-validation audits (completed)

### 8.1 Classified interface audit, corrected solutions

| r3 | n | R_abs | R_signed | med | p95 | max | near-zero | electrode |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| wall all | 896 | 0.014307 | 0.005725 | — | — | — | 120 | 168 |
| wall interface | 728 | 0.010914 | 0.009153 | 0.00425 | 0.01229 | 1.0000 | 0 | — |
| arc proper | 90 | 0.014412 | 0.013621 | 0.00712 | 0.01597 | 0.0223 | 0 | 0 |
| cap+corner | 158 | 0.042254 | 0.008647 | 0.00747 | 0.09023 | 0.2915 | 0 | 0 |
| all | 1,144 | 0.017702 | 0.004636 | — | — | — | 120 | 168 |

| r4 | n | R_abs | R_signed | med | p95 | max | near-zero | electrode |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| wall all | 1,246 | 0.012101 | 0.002523 | — | — | — | 168 | 232 |
| wall interface | 1,014 | 0.009510 | 0.005133 | 0.00317 | 0.00913 | 1.0000 | 0 | — |
| arc proper | 124 | 0.010703 | 0.010068 | 0.00501 | 0.01185 | 0.0200 | 0 | 0 |
| cap+corner | 214 | 0.033470 | 0.006866 | 0.00616 | 0.05710 | 0.2742 | 0 | 0 |
| all | 1,584 | 0.014568 | 0.002008 | — | — | — | 168 | 232 |

Interface-class convergence: wall 0.010914 → 0.009510 (ratio 0.871, p ≈ 0.40); arc proper 0.014412 → 0.010703 (0.743, p ≈ 0.84); cap+corner 0.042254 → 0.033470 (0.792, p ≈ 0.66). Both wall and arc decrease independently. All near-zero pairs are on the wall.

### 8.2 Same-mesh P1/P2 probe comparison, corrected

| group | n | rms dphi | max abs | mean |
|---|---:|---:|---:|---:|
| near drain (x > 2.38) | 21 | **0.002642 V** | 0.003836 V | +0.002537 V |
| other archived probes | 63 | 0.007524 V | 0.013345 V | +0.007119 V |
| all 84 archived | 84 | **0.006649 V** | 0.013345 V | +0.005974 V |
| corner neighbourhood (2.15, ±0.525) | 72 | **0.038869 V** | 0.062123 V | +0.037873 V |

The superseded rms was ≈ 2.06 V. The corrected same-mesh difference is 6.6 mV overall (×310 smaller), 2.6 mV near the drain, and 38.9 mV around the 79.51° corner — the corner remains the largest residual, consistent with reduced regularity but small in absolute terms.

Variant rebuild (matched physics, own meshes): P1 overshoot +0.232581 V; P2 variant +0.223689 V (3.8 % lower); probe rms 0.0187 V.

### 8.3 Disposition

- P2.1 variational, manufactured, linear-solve and interface checks: **pass on the corrected operator**.
- Interface classification and matched-field comparison: **pass**.
- P2.1 can be formally **reclosed**; P2.2 mathematical validation stands; the physical obstacle is **active**.
- P2.3 remains gated on one explicit model choice: whether the first physical solve uses the corrected standard P2 model (body slab charged silicon, pinned) or the P1-matched variant (body uncharged, equipotential).
