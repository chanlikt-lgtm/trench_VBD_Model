# P2.1 — Interface convergence, near-zero classification, and second resolution

**Date:** 19 October 2026
**Status:** P2.1 interface gate **met on the available evidence**; P2.2/P2.3 not started; P1 frozen; Stage 3 suspended.
**Diagnostic:** frozen — `p2_interface.py` (edge-matched 2-point Gauss, per-element P2 gradients, normal from the quadratic edge tangent, oriented Si→oxide by the Si element centroid), `p2_nearzero.py` (classification), `p2_final_table.py` (two-resolution table).

## 1. Count reconciliation (896 vs 980)

The earlier "810 same-sign" was a **fraction printed as a count**. Both classifications come from the same matched pairs:

| Set | Same sign | Opposite sign | Sum |
|---|---:|---:|---:|
| wall | 726 | 170 | **896** |
| arc | 248 | 0 | **248** |
| all | 974 | 170 | **1144** |

Same-sign fraction on the wall is 0.81027; 0.81027 × 896 = 726. No double counting, no differing subsets, no zero-handling difference. (`zero_pairs = 168` counts individual *sides* that are exactly zero and is not a pair count.)

## 2. Near-zero classification (r3)

Reference scale frozen: **D_ref = 282.983677**, defined as the median of D_scale = (|D_Si| + |D_ox|)/2 over all 1144 samples. The same D_ref is used on r4 so classification is identical.

| threshold | cutoff | near-zero | wall | arc | rel median | rel p95 | R_abs |
|---|---:|---:|---:|---:|---:|---:|---:|
| 1e-4 | 0.0283 | 103 | 103 | 0 | 0.00591 | 1.000 | 0.019182 |
| 1e-3 | 0.2830 | 127 | 127 | 0 | 0.00574 | 0.1338 | 0.019182 |
| 1e-2 | 2.8298 | 151 | 151 | 0 | 0.00528 | 0.0295 | 0.019182 |

R_abs is bit-identical at every threshold: all samples stay in the integral; only the relative distribution excludes near-zero pairs.

**Opposite-sign pairs:** 170, all on the wall. D_scale min 2.161e-6 (7.6e-9 D_ref), median 4.727e-3, p95 8.45, max 33.95. 127/170 lie below the 1e-3 cutoff; the 43 above it are all at x 0.409–0.559, |y| ≤ 0.363, with |D_Si| median 0.000 (exact zeros) and |D_ox| median 4.434 — i.e. they are the contact-adjacent class of §4, not genuine flux disagreement. Sign disagreement is confined to essentially zero flux plus the electrode footprint.

## 3. Second resolution (r4) — identical geometry, BCs, physics

| | r3 | r4 |
|---|---:|---:|
| h / h_tip | 0.014 / 0.0055 | 0.010 / 0.004 |
| elements | 75,970 | 147,629 |
| DOFs / free | 152,489 / 123,867 | 296,014 / 240,762 |
| R_free | 9.94e-15 | 1.30e-14 |
| phi max | 34.274208536 | 34.274208536 |
| interface samples | 1,144 | 1,584 |

### Interface class (electrode-adjacent excluded), frozen definitions

| metric | r3 | r4 | ratio | observed p |
|---|---:|---:|---:|---:|
| wall R_abs | 0.011079 | 0.009343 | 0.843 | 0.51 |
| wall median rel | 0.00468 | 0.00337 | 0.720 | 0.98 |
| wall p95 rel | 0.01346 | 0.00980 | 0.728 | 0.94 |
| arc R_abs | 0.033922 | 0.026612 | 0.785 | 0.72 |
| arc median rel | 0.00711 | 0.00468 | 0.658 | 1.24 |
| arc p95 rel | 0.04445 | 0.03031 | 0.682 | 1.14 |

Observed p uses h ratio 0.714 and one refinement step only: descriptive, not a certified order.

One-sided means are stable under refinement: wall interface −311.18/−314.33 → −310.61/−312.64 (0.2 %/0.5 %); arc −389.79/−387.94 → −390.51/−388.90 (0.2 %/0.2 %). The reconstruction itself is converged; the pointwise difference is what shrinks.

## 4. Electrode-adjacent class (separated, not discarded)

At r4, 60 of 1069 kept wall samples have rel > 0.99. **55 have an all-Dirichlet-pinned Si element and all-zero potential there**; x 0.4135–0.5637, y ±0.363, straddling the body-contact boundary at x = 0.54. Where the Si element is equipotential, grad(phi) = 0 identically, so eps_Si grad(phi).n = 0 while the oxide side carries flux: that flux terminates on the electrode, and the no-free-charge condition does not apply. r3 shows the same effect (43 bad, 41 fully pinned).

Consequently the electrode class returns R_abs = 2.000000 exactly — the algebraic maximum, since the numerator is |D_ox| and the denominator |D_ox|/2. This is a class definition artefact, not a solver result.

| class (r4 wall) | n | R_abs | median | p95 | mean D_Si | mean D_ox |
|---|---:|---:|---:|---:|---:|---:|
| all | 1246 | 0.011482 | 0.00385 | 1.000 | −252.78 | −254.08 |
| no-free-charge interface | 1014 | 0.009343 | 0.00337 | 0.00980 | −310.61 | −312.64 |
| electrode-adjacent | 232 | 2.000000 | 1.00000 | 1.000 | 0.00 | 1.85 |

## 5. Open items

1. Interface-class max remains 1.0000 at both resolutions (2 samples r3, 5 r4) — Si elements only *partially* pinned, at the contact edge. They need an edge-proximity classification rather than the all-pinned test.
2. Only one refinement step; the rate is ~0.5–1.0, which cannot distinguish O(h) from O(h^1.5).
3. ~~The quadratic-edge tangent and element Jacobian on the rounded tip have invariance evidence but no independent geometric check yet.~~ **Closed — see §6.**

## 6. Closure checks (requested 12:48)

### 6.1 r3 wall baseline: 0.013879 vs 0.011079 is a region change

Same 2-point quadrature, same |D| weights, same D_ref. Only the sample class differs.

| r3 wall | N | D | R_abs | n |
|---|---:|---:|---:|---:|
| all | 11.755640 | 847.029194 | 0.013879 | 896 |
| interface (electrode excluded) | 9.370778 | 845.836763 | 0.011079 | 728 |
| electrode-adjacent | 2.384862 | 1.192431 | 2.000000 | 168 |

N and D add exactly: 9.370778 + 2.384862 = 11.755640; 845.836763 + 1.192431 = 847.029194. r4 behaves identically (9.752656 = 7.927546 + 1.825111; 849.404011 = 848.491455 + 0.912555). Arc baseline unchanged throughout at 0.033922.

### 6.2 The seven high-mismatch samples

| tag | x | y | pinned/6 | D_Si | D_ox |
|---|---:|---:|---:|---:|---:|
| r3 | 0.5591 | −0.3632 | 0 | −2.3015 | +2.4976 |
| r3 | 0.5474 | +0.3627 | 1 | −12.333 | +55.564 |
| r4 | 0.5554 | −0.3630 | 0 | −2.2056 | +6.5559 |
| r4 | 0.5637 | +0.3634 | 0 | −0.6976 | +0.6198 |
| r4 | 0.5386 | −0.3623 | 5 | −9.0608 | +47.368 |
| r4 | 0.5338 | +0.3621 | 5 | −52.645 | +44.531 |
| r4 | 0.5386 | +0.3623 | 5 | −52.645 | +58.522 |

All lie within 0.024 µm of x = 0.54 at y ≈ ±a(0.54) = ±0.362: the gate-bottom edge (GATE_BOTTOM = 0.54) and the body-contact boundary (BODY_BOTTOM = 0.54) meeting the interface. Classification: **contact-corner**.

Correction to §4: rel = 1.0 is *not* evidence of a zero side. For opposite signs, |a−b|/(|a|+|b|) = 1 identically. Four of these seven have no pinned node and are ordinary opposite-sign pairs at the corner; three have a 5/6-pinned, near-equipotential Si element.

Excluding all-pinned **and** partially-pinned samples:

| | r3 | r4 | ratio | p |
|---|---:|---:|---:|---:|
| wall no-free-charge | 0.010368 (n=724) | 0.007757 (n=1010) | 0.748 | 0.86 |

### 6.3 Tip geometry

| check | result |
|---|---|
| C1 midside − chord midpoint | **0.000e+00** for every edge → straight-sided (`Mesh.SecondOrderLinear = 1`) |
| C2 162 arc nodes vs circle | max abs(r−R) = 1.93e-05 µm; mean 9.52e-06 µm |
| C3 tangent error, arc Gauss points | median 1.125° (chord geometry); max 88.4° only at the x = 2.15 junction |
| C4 det sign | Si 107,973 all negative; SiO₂ 23,378 all positive; gate/shield split 50/50; min abs(det) = 6.69e-06 |
| C5 grad(x), grad(y) on 160 tip elements | max error 5.08e-13 |
| C6 omitted circular segment area | 4.04e-06 µm² total over 80 edges (5.05e-08 per edge) |

C1 dissolves the question rather than answering it: there is no curved isoparametric mapping in this mesh. The quadratic edge is exactly the chord, so the chord tangent is the exact tangent of the mesh, and the affine Jacobian in `element_geom` is exact for the element as built. C5 then verifies that Jacobian and the P2 shape gradients on tip elements to 5e-13. C4's mixed determinant sign is a per-surface orientation convention, not degeneracy; the signed determinant is handled correctly, which C5 proves independently.

### 6.4 Structural finding: a 79.5° corner at x = 2.15

a(x) is C0 but not C1 at x = 2.15: the arc arrives with a vertical tangent (90.0°) and the cap leaves at 10.49°, a **79.51°** discontinuity. At x = 2.05 the kink is only 2.37°. This is frozen geometry and was not modified, but it belongs on the Stage 3 record: any tip-field claim sits on that corner.

### 6.5 The "arc" bucket contains three features

| region | r3 | r4 | ratio | p |
|---|---:|---:|---:|---:|
| arc proper (2.04 < x ≤ 2.14) | 0.014275 (n=90) | 0.010516 (n=124) | 0.737 | 0.91 |
| cap + corner (x > 2.14) | 0.046166 (n=158) | 0.036621 (n=214) | 0.793 | 0.69 |
| wall no-free-charge | 0.010368 (n=724) | 0.007757 (n=1010) | 0.748 | 0.86 |

The previously reported arc figure (0.033922 → 0.026612) mixes a finely meshed quarter arc with the coarser downstream cap and the corner. On the curved arc proper the mismatch is 1.05 % at r4, essentially the same as the flat wall's 0.78 %.
