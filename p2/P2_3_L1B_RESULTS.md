# P2.3 — L1b Results (`H_LOCAL = 0.0055 µm`, windows A/B unchanged)

**Record date:** 2026-10-20
**Status:** L1b **F7-certified**; P2.3 remains **OPEN**.
L1 preserved as the certified reference; P1 frozen; Stage 3 suspended.

## 1. F7 full-set gate — PASS

```
active rows        11,130
stationarity       3.481e-10     (threshold 1e-8)
complementarity    6.179e-17     (threshold 1e-10)
min slack         -2.001e-14     (threshold -1e-9)
violations         0
lambda min        +0.000e+00
oxide residual     4.480e-13     (threshold 1e-9)
overshoot         +1.224e-14     (threshold 1e-9)
GATE              PASS
```

Checkpoint: OSQP `solved`, 14,625 iterations, objective `-1.014015e+05`,
primal residual `1.385e-08`, dual residual `1.663e-09`; reload verified
(`x(350190,)`, `y(862794,)`, `idx(350190,)`).

## 2. Bounded neutral area

| Mesh | Lower | Upper | Width |
|---|---:|---:|---:|
| Intermediate | 0.021789286073 | 0.021793581617 | 4.296e-06 |
| **L1b** | **0.021842101270** | **0.021845749070** | **3.648e-06** |
| L1 | 0.022071011922 | 0.022075109085 | 4.097e-06 |

L1b lies **between** intermediate and L1. The area is monotone in realised
local resolution.

### 2.1 Interval-bounded increments

| Comparison | Lower | Upper |
|---|---:|---:|
| intermediate → L1b | 4.851965e-05 | 5.646300e-05 |
| L1b → L1 | 2.252629e-04 | 2.330078e-04 |
| intermediate → L1 | 2.774303e-04 | 2.858230e-04 |

The two partial increments sum to the direct intermediate → L1 interval, as
they must.

### 2.2 The increments do not shrink with resolution

```
mid -> L1b : 5.25e-05   (realised h 0.00780 -> 0.00524, factor 0.672)
L1b -> L1  : 2.29e-04   (realised h 0.00524 -> 0.00399, factor 0.762)
ratio      : [3.990, 4.802]
```

A *smaller* refinement step in `h` produces a **4.0–4.8× larger** area
increment. This is the opposite of a convergent sequence: adding L1b as an
intermediate resolution shows the area is monotone but steepening as the
mesh is refined. No Richardson-style extrapolation is justified.

## 3. Front cuts at the two priority locations

Threshold `eps = 1e-8 V`; intervals in µm.

**x = 2.35**

```
intermediate  (-0.0754, +0.0686)
L1b           (-0.0774, +0.0719)
L1            (-0.0718, +0.0684) + (+0.0722, +0.0763)
```

The lower edge is **non-monotone**: −0.0754 → −0.0774 → −0.0718. L1b moves
0.0020 µm further upstream than the intermediate mesh, then L1 pulls back
0.0056 µm. L1b does not fall between its neighbours.

**x = 2.41**

```
intermediate  (-0.1895, +0.1903) + (+0.1922, +0.1981)
L1b           (-0.2019, +0.1945)            [single interval at 1e-8]
L1            (-0.1996, +0.1943) + (+0.1964, +0.1982)
```

The lower edge moves outward monotonically (−0.1895 → −0.2019 → −0.1996),
with L1b marginally beyond L1. The shallow second interval present in both
neighbours collapses to a single interval in L1b at `1e-8`, reappearing only
at `1e-5` as `(+0.1962, +0.1975)`.

## 4. Context cuts (threshold sweep)

`eps = 1e-8 / 1e-6 / 1e-5 V`

| x | 1e-8 | 1e-6 | 1e-5 |
|---|---|---|---|
| 2.35 | (-0.0774, 0.0719) | same | same |
| 2.38 | (-0.1350, 0.1395) | same | (-0.1363, 0.1395) |
| 2.39 | (-0.1523, 0.1562) | same | (-0.1573,-0.1560),(-0.1523,0.1562),(0.1595,0.1605) |
| 2.40 | (-0.1798, 0.1730) | (-0.1798,0.1730),(0.1757,0.1782) | (-0.1853,-0.1842),(-0.1798,0.1785) |
| 2.41 | (-0.2019, 0.1945) | same | (-0.2019,0.1945),(0.1962,0.1975) |

The extra intervals that appear at `1e-5` are boundary artefacts adjacent to
the main intervals, not separate neutral pockets; they were filtered
(`short` counts 1–3 per cut).

## 5. Interpretation

1. **Area:** L1b confirms monotonicity in resolution but **not** convergence.
   The increment grows by 4.0–4.8× as `h` decreases, so the four heterogeneous
   meshes still do not define a limit.
2. **Front:** the upstream cut at x = 2.35 is non-monotone even with the
   intermediate point inserted. Local element size alone does not control it.
3. **Cause:** unchanged from the L1 assessment — mesh-distribution sensitivity
   is demonstrated; refinement contrast is **not** proven to be the mechanism.
   L1b varies local size and transition placement together, so no unique
   attribution is available from this pair.
4. **Known bias:** the whole-element centroid donor load differs between L1
   and L1b (see `P2_3_L1B_INPUT_AUDIT.md` §4). That difference enters the load
   vector and can move an active-set boundary; it is not corrected here by
   design.

**No continuum limit, no converged front position, and no converged tip field
follow from these results.**

## 6. Method notes (defects found and fixed during this run)

Three defects in the L1b polish path were found and corrected before the
certified result was produced. They are recorded because the same code shaped
the L1 polish.

1. **Active-set seeding.** The first version seeded `W` from primal slack only.
   OSQP returns `y ≤ 0` on active lower bounds, so the dual must be read as
   `λ = −sK·y`. Seeding from `λ > 1e-8` (as in the L1 pipeline) is correct.
2. **Dense independence filter.** `independent()` forms a dense `n×n` Gram
   matrix and runs a pivoted QR — at 11,802 rows this is 1.1 GB and O(n³).
   It accounted for a 50-minute stall. The Bernstein rows are element-local
   and independent by construction; the filter is unnecessary and was removed.
3. **Residual-correction sign.** The correction used `[-stat, -slackW]`; the
   correct residual of `[K u − Cᵀλ; −C u] = [g; cD]` is `[-stat, +slackW]`.
   With the wrong sign the constraint residual **doubled every round**
   (−1.79e-9 → −3.47e-9 → −6.93e-9 → −1.39e-8 → −2.77e-8 → −5.55e-8) while
   stationarity fell, which is what made the polish appear to stall. With the
   sign corrected both quantities converge:

```
refine 0   stat 4.793e-05   slack -1.790e-09
refine 1   stat 3.771e-06   slack -3.875e-10
refine 2   stat 7.321e-08   slack -4.702e-12
refine 3   stat 4.695e-09   slack -3.629e-13
refine 4   stat 3.481e-10   slack -2.001e-14
```

Note: `p2_L1_scaled.py` carries the same `-slack[W]` line. L1's certification
is **unaffected** — the F7 gate recomputes stationarity, feasibility,
complementarity, oxide residual and overshoot independently from the final
iterate, and L1 passed all of them. But that script should carry the corrected
sign if it is reused.

## 7. Artifacts

- `p2_L1b_p2mesh.npz`, `p2_L1b_mesh.json` — mesh and quality record
- `p2_L1b_osqp_checkpoint.npz`, `p2_L1b_osqp.json` — primal/dual checkpoint
- `p2_L1b_scaled.npz`, `p2_L1b_scaled.json` — certified solution and gate
- `p2_L1b_post.json` — area bounds and cuts
- `p2_charge_mask_audit.json`, `p2_size_audit.json` — input audits
