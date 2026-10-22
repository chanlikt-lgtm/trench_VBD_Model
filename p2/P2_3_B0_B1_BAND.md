# P2.3 — B0 → B1: Band Refinement at Fixed Background

**Record date:** 2026-10-22
**Question:** with geometry, charged region, background target and generation
policy held fixed, what does enabling the frozen-band refinement do to the
bounded neutral area?

## 1. Control design

| Quantity | B0 | B1 |
|---|---|---|
| Size policy | uniform 0.0085 | band 0.0040, threshold 0.04→0.06, else 0.0085 |
| Nodes / elements | 212,779 / 105,774 | 231,568 / 114,721 |
| Near-band elements | 2,282 @ 0.00850 | 10,160 @ 0.00400 |
| Si area | 1.954559924360 | 1.954559924360 |
| Exact clipped charged area | 1.575269192653 | 1.575269192653 |
| Centroid charged area | 1.575255252958 | 1.575498331328 |
| Min angle | 35.95° | 31.09° |

Si and exact charged areas are identical to all printed digits, so the polygonal
geometry and arc discretisation are the same. The band adds 8,947 elements, of
which 7,878 (88 %) lie inside the 0.045 µm band.

## 2. Certification

**B0** — OSQP `solved`, 5,025 iterations, obj `-1.052110e+05`,
pri `4.051e-08`, dua `3.763e-09`, checkpoint reload verified
(`x(142473,)`, `y(303672,)`).

```
active rows        4,435
stationarity       5.824e-11
complementarity    5.016e-17
min slack         -8.147e-15
violations         0
lambda min        +0.000e+00
oxide residual     4.083e-13
overshoot         +4.268e-15
GATE              PASS
```

**B1** — certified in the previous record (stat 8.929e-11).

## 3. Bounded neutral area

| Mesh | Lower | Upper |
|---|---:|---:|
| B0 | 0.021644935001 | 0.021649266844 |
| B1 | 0.021987652318 | 0.021991846424 |

```
B0  -> B1    [+3.383855e-04, +3.469114e-04]   +1.563 % … +1.603 %
L1b -> B1    [+1.419032e-04, +1.497452e-04]   +0.650 % … +0.686 %
L1b -> L1    [+2.252629e-04, +2.330078e-04]   +1.031 % … +1.067 %
```

The band alone moves the area by 1.5× the whole L1b→L1 refinement step, and
2.4× the L1b→B1 difference. At fixed background, boundary resolution from
0.0085 to 0.0040 is the single largest area lever measured in this project.

## 4. Donor-load confound — quantified, and too small

Body-cut audit at x = 0.54 µm:

| Quantity | B0 | B1 |
|---|---:|---:|
| Straddling elements | 180 | 174 |
| Straddle area | 5.598270e-03 | 5.389783e-03 |
| Net mask error | -1.394e-05 | +2.291e-04 |
| Overcharge | +5.203172e-04 | +6.856357e-04 |
| Undercharge | +5.342569e-04 | +4.564970e-04 |
| h near cut (±0.05 µm) | 4,425 @ 0.00850 | 4,421 @ 0.00850 |
| Σf = RHO·A_centroid | 1.966806755e+03 | 1.967110253e+03 |

ΔΣf = +3.035e-01, i.e. **1.54e-4 relative**, against an area change of
**1.56–1.60 % relative** — a factor of ~100.

The mechanism is a near-cancellation: the over- and under-charge terms are each
~5e-4 µm², five times their difference, and the net flips sign between B0 and
B1 (−1.39e-5 → +2.29e-4) even though the element size at the cut is unchanged
at 0.0085 µm. The downstream band changes the *triangulation arrangement*, which
changes which whole elements straddle the body cut. The centroid rule is thus
arrangement-sensitive, not size-sensitive.

**Consequence:** the donor change cannot explain the B0→B1 area difference.
It remains a declared confound at the 0.015 % level and would need the exact
partial-element load to be removed entirely; it does not need to be removed to
interpret this result.

## 5. Front coordinates

ε = 1e-8:

```
x = 2.35   B0  (-0.0649, +0.0685)            single
           B1  (-0.0746, -0.0727) + (-0.0698, +0.0757)

x = 2.41   B0  (-0.2026, +0.1957)            single
           B1  (-0.2017, -0.2005) + (-0.1972, +0.1975)
```

The band-refined front is wider: the upstream edge moves from −0.0649 to
−0.0746 (0.0097 µm) at x = 2.35. B1 also shows the shallow satellite spans,
which merge by ε = 1e-5 and sit at the band edge (0.0040 → 0.0085 transition);
they remain unverified as physical features.

## 6. Findings

1. **Boundary resolution is the dominant area control.** At fixed background,
   refining the boundary from 0.0085 to 0.0040 raises the bounded area by
   1.56–1.60 %, more than the entire L1b→L1 step.
2. **Donor discretisation is a 100× smaller effect** on the same comparison and
   is quantified rather than assumed.
3. **The centroid mask is arrangement-sensitive.** Its net error is a small
   difference of two large opposing terms, so it moves by 2.4e-4 when the mesh
   layout changes far away — an important property for any future mesh study.
4. **No convergence order follows.** B0→B1 is a single jump (factor 2.1 in h).
   It establishes direction and magnitude, not an asymptotic rate.
5. The band-refined front is wider, consistent with the area increase, but
   front coordinates remain layout- and threshold-sensitive.

## 7. Artifacts

- `p2_mesh_B0.py`, `p2_B0_p2mesh.npz`, `p2_B0_mesh.json`
- `p2_B0_osqp_checkpoint.npz`, `p2_B0_osqp.json`
- `p2_B0_scaled.npz`, `p2_B0_scaled.json`
- `p2_B0_post.json`, `p2_body_cut_audit.json`
- commit `e473a2c`
