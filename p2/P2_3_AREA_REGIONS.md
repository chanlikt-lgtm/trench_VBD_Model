# P2.3 — Spatially Decomposed Neutral-Area Bounds

**Record date:** 2026-10-21
**Scope:** mid, L1b, M2, L1 — all F7-certified. Fixed physical regions,
independent of mesh connectivity. No new mesh was built.

## 1. Method

Regions are fixed rectangles in silicon:

```
A    : 2.320 <= x <= 2.370,  |y| <= 0.13
B    : 2.395 <= x <= 2.420,  |y| <= 0.24
REM  : charged Si outside A and B
```

Bernstein bounds are computed by the same de Casteljau subdivision used for
the whole-domain enclosure (`AREA_TOL = 1e-11`, `MAX_DEPTH = 24`,
`EPS = 1e-8`), extended so that each sub-triangle is classified against the
region using its physical vertices:

- fully outside → discarded;
- fully inside → ordinary coefficient test (all ≤ 0 → lower bound; all > 0 → nothing; mixed → subdivide);
- straddling → never counted in a lower bound, subdivided until its area is below `AREA_TOL`, then counted in the upper bound only.

Ownership is exclusive: a sub-triangle belongs to exactly one of A, B, REM,
so no region is double-counted.

## 2. Conservation check

Sum of the three regional bounds against the independent whole-domain
enclosure (same solutions, same `EPS`, same tolerance):

| Mesh | Σ regional lower | Σ regional upper | whole lower | whole upper | Δ lower | Δ upper |
|---|---:|---:|---:|---:|---:|---:|
| mid | 0.021787102964 | 0.021799962373 | 0.021789286073 | 0.021793581617 | -2.18e-06 | +6.38e-06 |
| L1b | 0.021840692831 | 0.021850004899 | 0.021842101270 | 0.021845749070 | -1.41e-06 | +4.26e-06 |
| M2 | 0.021838121024 | 0.021847402993 | 0.021839484478 | 0.021843096605 | -1.36e-06 | +4.31e-06 |
| L1 | 0.022068931948 | 0.022081474801 | 0.022071011922 | 0.022075109085 | -2.08e-06 | +6.37e-06 |

In every case `Σ lower ≤ whole lower` and `Σ upper ≥ whole upper`. The sums
bracket the whole-domain value, with 1.4e-6 … 2.2e-6 slack below and
4.3e-6 … 6.4e-6 slack above, attributable to region-boundary straddling.

## 3. Regional bounds

| Region | mid | L1b | M2 | L1 |
|---|---:|---:|---:|---:|
| A | [0.005219248, 0.005222418] | [0.005184914, 0.005187103] | [0.005194455, 0.005196714] | [0.005267099, 0.005270278] |
| B | [0.009557472, 0.009560167] | [0.009610430, 0.009612156] | [0.009562626, 0.009564329] | [0.009664305, 0.009666824] |
| REM | [0.007010382, 0.007017377] | [0.007045349, 0.007050746] | [0.007081040, 0.007086360] | [0.007137528, 0.007144373] |
| Σ | [0.021787103, 0.021799962] | [0.021840693, 0.021850005] | [0.021838121, 0.021847403] | [0.022068932, 0.022081475] |

## 4. Change intervals

Interval subtraction `[a⁻−b⁺, a⁺−b⁻]`:

**mid → L1b**

| Region | Interval | Excludes zero |
|---|---|---|
| A | [-3.7504e-05, -3.2146e-05] | yes (negative) |
| B | [+5.0263e-05, +5.4684e-05] | yes |
| REM | [+2.7972e-05, +4.0364e-05] | yes |

**L1b → M2 (layout only)**

| Region | Interval | Excludes zero |
|---|---|---|
| A | [+7.3528e-06, +1.1800e-05] | yes |
| B | [-4.9530e-05, -4.6101e-05] | yes (negative) |
| REM | [+3.0294e-05, +4.1011e-05] | yes |
| Σ | [-1.1884e-05, +6.7102e-06] | **no** |

**L1b → L1**

| Region | Interval | Share of Σ |
|---|---|---|
| A | [+7.9996e-05, +8.5364e-05] | 36.0 % |
| B | [+5.2150e-05, +5.6394e-05] | 23.6 % |
| REM | [+8.6781e-05, +9.9024e-05] | 40.4 % |
| Σ | [+2.1893e-04, +2.4078e-04] | — |

**M2 → L1**

| Region | Interval |
|---|---|
| A | [+7.0385e-05, +7.5822e-05] |
| B | [+9.9976e-05, +1.0420e-04] |
| REM | [+5.1168e-05, +6.3333e-05] |
| Σ | [+2.2153e-04, +2.4335e-04] |

## 5. Findings

1. **The L1 increase is distributed, not localised.** Roughly 36 % in window A,
   24 % in window B, 40 % in the silicon remainder. No single front segment
   carries the change, so a front-only refinement strategy is not justified.
2. **Relative growth is highest in A.** Against its own base area: A +1.59 %,
   REM +1.32 %, B +0.56 % at L1b → L1. Window A, the most strongly refined
   region, also grows the most in relative terms, but B is nearly flat.
3. **Refinement redistributes area between regions.** mid → L1b *decreases*
   region A by 3.2–3.8e-5 µm² while B and REM increase. A region-level
   increment can therefore be negative even when the total grows — the
   neutral set moves relative to the fixed analysis boxes.
4. **Layout change is a pure redistribution.** L1b → M2 moves 4.6–5.0e-5 µm²
   out of B into REM (+3.0–4.1e-5) and A (+0.7–1.2e-5), with the total
   interval straddling zero. This is the regional signature of the
   area-invariance already reported: the same total area, placed differently.
5. **Conservation is verified.** The regional sums bracket the independent
   whole-domain enclosure on all four meshes, so the decomposition does not
   lose or duplicate neutral area beyond the stated boundary slack.

## 6. Limitations

- A, B, REM are analysis boxes, not physical features. Their boundaries are
  arbitrary; item 3 shows that region-level increments depend on that choice.
- Boundary straddling contributes 4e-6 … 6e-6 µm² of slack to the upper sums.
  Regional changes of that order are not resolvable.
- No continuum limit is implied. The decomposition localises *where* the
  discrete area changes, not what the limit is.
- The M2 window-B core is 3 % coarser than L1b; the L1b → M2 regional split
  therefore mixes layout and a small resolution difference. The total
  (straddling zero) is insensitive to it.

## 7. Artifacts

- `p2_area_regions.py`, `p2_area_regions.json`, `p2_area_regions_out.txt`
- `_region_deltas.py` — change intervals and shares
