# P2.3 — Front-Conditioning Diagnostic

**Record date:** 2026-10-20
**Scope:** five F7-certified P2 solutions (r3, r4, intermediate, L1b, L1);
cuts x = 2.35 and 2.41 µm. No new mesh was built.

## 1. Interpolation check

P2 mid-side ordering detected geometrically before use:
`((0,1),(1,2),(2,0))` for all five meshes, max |node − edge midpoint| = 0.0.
Gradients below therefore come from the full quadratic interpolant matching
`A.shape_vals`, not from vertex differences.

## 2. Threshold sensitivity of the extracted intervals

`eps = 1e-8 / 1e-6 / 1e-5 V`, main intervals only.

**x = 2.35**

| Mesh | 1e-8 | 1e-6 | 1e-5 |
|---|---|---|---|
| r3 | (-0.0624, 0.0542)+(0.0703, 0.0756) | same | extra (-0.0805,-0.0794) appears |
| r4 | (-0.0697, 0.0714) | same | same |
| mid | (-0.0754, 0.0686)+(0.0787, 0.0799) | same | same |
| L1b | (-0.0774, 0.0719) | same | same |
| L1 | (-0.0718, 0.0684)+(0.0722, 0.0763) | upper 0.0690 | upper 0.0704 |

**x = 2.41**

| Mesh | 1e-8 | 1e-6 | 1e-5 |
|---|---|---|---|
| r3 | (-0.1945, 0.1958) | same | extra (0.2068,0.2083) |
| r4 | (-0.1892, 0.1899)+(0.2037, 0.2050) | same | same |
| mid | 3 intervals: (-0.1977,-0.1914), (-0.1895,0.1903), (0.1922,0.1981) | same | merge to one |
| L1b | (-0.2019, 0.1945) | same | extra (0.1962,0.1975) |
| L1 | (-0.1996,0.1943)+(0.1964,0.1982) | 2nd → 0.2022 | 2nd → 0.2032 |

Threshold sensitivity is confined to the secondary spans and to the L1 upper
edge; the primary lower edges are stable to 1e-5 V.

## 3. Endpoint gradient and conditioning

`|δy| ≈ |δu| / |∂u/∂y|`, evaluated at each endpoint from the P2 gradient.

At x = 2.35, `|∂u/∂y|` spans **3.1e-4 … 1.8e-1 V/µm** (amplification
6 … 3.2e3 µm/V). At x = 2.41 it spans **5.5e-6 … 2.2e-1 V/µm**
(amplification 5 … 1.8e5 µm/V).

Worked example, x = 2.35, lower edge:

```
L1b endpoint  y = -0.0774
du between meshes at that y        -1.80e-04 V
du/dy of L1 at that y              +3.999e-02 V/um
predicted displacement             -0.0045 um
observed L1b -> L1 displacement    +0.0056 um
```

Magnitudes agree to ~20 %; the sign does not. Across the table predicted and
observed shifts are of the same order (1e-3 … 2.5e-2 µm), so the observed
front scatter is consistent with amplification of **~2e-4 V** inter-mesh
potential differences rather than with a large physical difference.

## 4. Where the linearization fails

Several endpoints sit on a plateau where `|∂u/∂y|` falls to 1e-13 … 1e-15
(e.g. mid at y = −0.0754 and +0.0686, x = 2.35; mid at −0.1977, −0.1914,
+0.1903, +0.1922, x = 2.41). At those points the predicted displacement is
unbounded and the estimate is meaningless. These are plateau-adjacent
endpoints, not isolated simple crossings — the same fragments that appear and
vanish with threshold and mesh.

## 5. Gap maxima — are adjacent intervals distinguishable?

Maximum u inside each gap (V):

| Cut | r3 | r4 | mid | L1 | L1b |
|---|---|---|---|---|---|
| 2.35 | 2.9e-4 | single | 7.1e-4 | 4.2e-4 | single |
| 2.41 | single | 7.5e-4 | 1.4e-5, 1.2e-5 | 1.1e-5 | single |

Cross-mesh potential difference near the front is 1.9e-3 – 3.6e-3 V (95th
percentile) and 2e-5 – 4e-4 V at the endpoints themselves.

At **x = 2.35** the gaps (2.9e-4 – 7.1e-4 V) exceed the endpoint-level Δu, so
the two-span structure there is resolvable in principle.

At **x = 2.41** the gaps in mid and L1 are 1.1e-5 – 1.4e-5 V — an order of
magnitude *below* the inter-mesh potential difference. Those secondary spans
are therefore **not** numerically resolvable, and the changing interval count
(1 → 2 → 3 → 1 → 2) should not be read as a topological feature.

## 6. Cross-mesh potential differences

| Pair | max abs Δu (all y) | max in ±0.02 of front | p95 in band |
|---|---:|---:|---:|
| L1b − L1, 2.35 | 5.64e-2 | 4.35e-3 | 1.87e-3 |
| mid − L1b, 2.35 | 5.28e-2 | 5.69e-3 | 3.28e-3 |
| mid − L1, 2.35 | 5.42e-2 | 6.24e-3 | 3.17e-3 |
| L1b − L1, 2.41 | 1.32e-2 | 3.35e-3 | 1.78e-3 |
| mid − L1b, 2.41 | 1.53e-2 | 7.16e-3 | 3.43e-3 |
| mid − L1, 2.41 | 1.08e-2 | 6.41e-3 | 3.56e-3 |

The large differences (5e-2 V) occur in the bulk at y ≈ ±0.33 µm, far from the
front; near the front the solutions agree to ~2e-3 V.

## 7. Conclusions

1. The extracted front is **moderately ill-conditioned**: with
   `|∂u/∂y| ~ 4e-2 V/µm`, a 2e-4 V potential difference displaces the root by
   ~5e-3 µm — the size of the observed L1b/L1 discrepancy.
2. The non-monotone front coordinates are therefore **consistent with
   root-location amplification** of small inter-mesh potential differences.
   This does not prove the mechanism; it does establish that the front is not
   a reliable discriminator at this accuracy.
3. Plateau-adjacent endpoints are degenerate: the linearised sensitivity is
   unbounded and the reported endpoint is threshold-dependent.
4. At x = 2.41 the secondary intervals are separated by less than the
   inter-mesh potential difference, so the interval count there carries no
   information.
5. The area statistic is more robust than the front position — the gaps and
   the aggregation suppress single-endpoint amplification — but it is still
   not converged, and its increments are non-regular (§2 of the L1b record).

## 8. Implication for experiment choice

Because the front moves by ~5e-3 µm in response to ~2e-4 V, a mesh experiment
aimed at the front must either (a) hold the realised core resolution fixed and
change only transition grading/connectivity, or (b) reduce the potential
discretisation error itself (exact partial-element donor load, tighter F7) so
that Δu drops below ~1e-5 V. Another `H_LOCAL` reduction alone would change
both the resolution and the error, and would not separate them.

## 9. Artifacts

- `p2_front_conditioning.py`, `p2_front_conditioning.json` — endpoints,
  gradients, gaps, cross-mesh differences
- `p2_front_amp.json` — predicted vs observed displacement per endpoint
- `p2_front_du_band.json` — Δu in the front-adjacent band
- `p2_front_cond_out.txt` — full console record
