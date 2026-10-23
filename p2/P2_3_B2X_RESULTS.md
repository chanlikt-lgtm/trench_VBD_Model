# P2.3 — B2-X: Second Band Level, Exact Source

**Record date:** 2026-10-23
**Question:** with the band already active, does a further halving of `h_band`
produce a smaller area increment than the first activation did?

## 1. Configuration (frozen, as authorised)

```
band reference   original frozen L1b contour (unchanged)
h_band           0.00283 um      (level 1: 0.0040)
background       0.0085 uniform
grading/tip      identical to B1-X
donor            exact clipped P2 integration
constraints      unchanged Bernstein set
solver           checkpoint-first OSQP + scaled MINRES
area tolerance   1e-13
```

### 1.1 Defect caught by the pre-solve gate

The first B2 build requested `h_band = 0.00283` and realised **0.00340**.
`Mesh.MeshSizeMin` was pinned at `H_TIP = 0.0034`, so gmsh clamped the band
target. At B1 this was invisible because `0.0040 > 0.0034`; it binds only at
level 2. Corrected to `min(H_TIP, H_BAND)` and rebuilt:

```
requested 0.00283  ->  realised A 0.00283  B 0.00283  strip 0.00283
```

Without the realised-size check, B2 would have been a 0.0034 mesh recorded as
0.00283, and the increment test would have compared the wrong `h`.

## 2. Mesh verification

```
nodes 253,636    elements 125,579
min area 1.804e-06   max aspect 2.70   min angle 31.26 deg
Si area 1.954559924360          exact clipped charged 1.575269192653
straddling elements 177
```

Si and exact charged areas are bit-identical to B0 and B1 — geometry and donor
region unchanged; only the band level moved.

Donor conservation:

```
full-element identity max err  6.35e-22
sum f / rho  1.575269192653   =  exact charged area   (match: True)
```

## 3. Certification

```
[B2X-A] OSQP solved, 5,450 iterations, obj -9.485464e+04
        pri 1.228e-08, dua 3.653e-09, checkpoint reload verified
        x(181699,) y(416940,)

active rows      31,723
stationarity     1.0195e-10
complementarity  8.082e-18
min slack       -3.228e-15     violations 0
lambda min       0             oxide 4.509e-13
overshoot        1.679e-15
GATE             PASS
```

## 4. Bounded areas at 1e-13

```
B0-X  [0.021645551625654, 0.021645845589193]   width 2.940e-07
B1-X  [0.021988322072981, 0.021988665490453]   width 3.434e-07
B2-X  [0.022154236625040, 0.022154596319621]   width 3.597e-07   697 s, 887 mixed
```

Depth limit again not binding: `log4(1.804e-06/1e-13) ≈ 13.3 < 24`.

## 5. The increments

```
ΔA_01^X = A_B1X − A_B0X  [+3.424764838e-04, +3.431138648e-04]
ΔA_12^X = A_B2X − A_B1X  [+1.655711346e-04, +1.662742466e-04]

ratio ΔA_12/ΔA_01 = [0.4826, 0.4855]
strictly positive: yes
```

The second increment is **48 % of the first** — it halved, and both intervals
are disjoint from zero and from each other.

## 6. Apparent exponent — observation, not a claim

With `h` ratio `0.00283/0.0040 = 0.7075`, a smooth asymptotic law
`A(h) = A* − C h^p` would give increment ratio `(h1/h0)^p`. The measured ratio
implies

```
p ≈ ln(0.483)/ln(0.7075) ≈ 2.09 … 2.11
```

Consistent with a second-order trend. It is **not** a measured convergence
order: there are only two increments, B0-X carries no active band (so the
triple is not a nested family), and the band geometry, transition grading and
active set all differ between levels. Treat `p ≈ 2` as a hypothesis for a future
designed study, not a result.

## 7. Decision

Per the pre-registered rule, a smaller increment is *encouraging evidence of
diminishing sensitivity to band refinement*. It does not justify a convergence
order, and **no third level is authorised** regardless of this outcome.

## 8. Artifacts

- `p2_mesh_B2.py`, `p2_B2_p2mesh.npz`, `p2_B2_mesh.json`
- `p2_B2X_osqp_checkpoint.npz`, `p2_B2X_osqp.json`
- `p2_B2X_scaled.npz`, `p2_B2X_scaled.json`
- `p2_area_a13_b2x.json`, `_b2x_cmp.py`

## 9. Closure (2026-10-23)

The B0-X → B1-X → B2-X comparison is **accepted and closed** at commit
`2478f82`. The second increment is strictly positive and its supremum
(`1.662742466e-04`) is below the infimum of the first (`3.424764838e-04`), so the
reduction is certified in the interval sense, not merely in midpoints.

**No third band level is authorised, and none is justified by this result.**
The sequence is not a nested family (B0-X has no active band), so its successive
differences must not be used as Richardson data. A future convergence study
requires a genuinely comparable mesh family under one frozen policy.

B2-X full-set F7 diagnostics, archived for audit (`p2_B2X_scaled.json`):

```
active rows      31,723
stationarity     1.0195393954383166e-10
complementarity  8.082037731392984e-18
min slack       -3.2276814617319727e-15
violations       0
lambda min       0.0
oxide residual   4.5094412746102277e-13
overshoot        1.6785376981639928e-15
minres info      0
gate             true
```

Defensible conclusion, no stronger claim:

> Under the frozen L1b-derived boundary-band policy and exact partial-element
donor integration, reducing the band target from 0.0040 to 0.00283 µm produced a
further positive neutral-area increment whose magnitude is less than half the
preceding B0-X → B1-X increment. This establishes diminishing area sensitivity
across the tested sequence, but does not establish continuum convergence or an
asymptotic convergence order.

**Disposition:** P2.3 spatial continuum convergence **OPEN**; P1 **FROZEN**;
Stage 3 **SUSPENDED**. The `15ec811` eight-endpoint audit remains the comparison
baseline.
