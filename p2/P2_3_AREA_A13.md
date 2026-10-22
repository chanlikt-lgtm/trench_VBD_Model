# P2.3 — AREA_TOL = 1e-13 Four-Area Audit

**Record date:** 2026-10-22
**Scope:** postprocessing only. No OSQP solve, no polish, certified solutions
untouched. Separate outputs; the `1e-11` enclosures are preserved for
comparison.

## 1. Configuration

```
AREA_TOL = 1e-13     MAX_DEPTH = 24
method   : Bernstein subdivision, unchanged from p2_area_bounds
```

**Depth-limit status: not reached.** Termination needs sub-area < 1e-13; from
each mesh's smallest element area the maximum terminating depth is 13.64 (B0,
B0-X; A_min 1.641e-05) and 12.56 (B1, B1-X; A_min 3.664e-06). Both are below
`MAX_DEPTH = 24`, so no conversion is depth-truncated. Every bound reported here
is therefore area-limited, as intended.

## 2. Eight endpoints (CP-A1, CP-A2)

| Solution | Lower | Upper | Width | Runtime | Mixed |
|---|---:|---:|---:|---:|---:|
| B0 | 0.021645551687359 | 0.021645845650898 | 2.940e-07 | 853 s | 294 |
| B1 | 0.021988339286914 | 0.021988682704330 | 3.434e-07 | 786 s | 629 |
| B0-X | 0.021645551625654 | 0.021645845589193 | 2.940e-07 | 855 s | 294 |
| B1-X | 0.021988322072981 | 0.021988665490453 | 3.434e-07 | 780 s | 629 |

Total runtime 3,275 s. Stored in `p2_area_a13.json` with the `1e-11` values
kept in `p2_area_bounds.json` / the per-mesh post records.

**Nesting check.** Each 1e-13 interval lies strictly inside its 1e-11
counterpart:

```
B0   1e-11 [0.021644935001, 0.021649266844]  4.33e-06
     1e-13 [0.021645551687, 0.021645845651]  2.94e-07   inside
B1   1e-11 [0.021987652318, 0.021991846424]  4.33e-06
     1e-13 [0.021988339287, 0.021988682704]  3.43e-07   inside
```

Width reduction ~14.5× for a 100× tolerance reduction, consistent with the
`perimeter × sqrt(TOL)` scaling.

## 3. Conservative propagation (CP-A3, CP-A4)

Interval subtraction `[L_a − U_b, U_a − L_b]` throughout; no midpoints used as
exact values.

```
ΔA_C  = A_B1  − A_B0    [+3.424936360e-04, +3.431310170e-04]   width 6.374e-07
ΔA_X  = A_B1X − A_B0X   [+3.424764838e-04, +3.431138648e-04]   width 6.374e-07

Δ_src = ΔA_X − ΔA_C     [-6.545331832e-07, +6.202287835e-07]   width 1.275e-06
```

The `Δ_src` width is exactly the sum of the two Δ widths, as conservative
subtraction requires.

Per-solution source corrections, for completeness:

```
δ_B0 = A_B0X − A_B0   [-6.17e-11, -6.17e-11]   (below both widths)
δ_B1 = A_B1X − A_B1   [-1.72e-08, -1.72e-08]   (below both widths)
```

## 4. Zero inclusion (CP-A5)

```
0 ∈ Δ_src :  TRUE
```

The correction interval contains zero. Per the pre-registered decision rule:

> the sign remains unresolved; **do not claim a measured nonzero correction.**

Express as a bound only: any source-induced change to the band increment is
smaller than **±6.5e-7 µm²**, i.e. **±0.19 %** of the increment (~3.43e-4),
under the stated numerical assumptions.

## 5. What changed relative to the 1e-11 result

| Quantity | 1e-11 | 1e-13 |
|---|---|---|
| ΔA_C | [3.383855e-04, 3.469114e-04] | [3.424936e-04, 3.431310e-04] |
| correction bound | ±8.5e-06 (±2.5 %) | ±6.5e-07 (±0.19 %) |
| 0 ∈ Δ_src | yes | **yes** |

The verdict did not change: the source correction is unresolved at both
tolerances. What improved is the *bound* on it — 13× tighter — and the
certified increment itself is now 13× more precise.

## 6. Validity notes

- The tighter intervals are nested inside the looser ones on both meshes, which
  is the expected behaviour for a genuine narrowing of the same enclosure.
- No depth truncation occurred (§1), so the bounds are limited only by
  `AREA_TOL` and by the Bernstein convex-hull test itself.
- The `1e-11` results remain valid as looser enclosures; they are not superseded
  as errors, only as precision.
- Certified solutions, meshes and checkpoints were read-only in this audit.

## 7. Artifacts

- `p2_area_a13.py`, `p2_area_a13.json`, `p2_area_a13_out.txt`
- `_a13_prop.py` (propagation), `_depth_bound.py` (depth arithmetic)
