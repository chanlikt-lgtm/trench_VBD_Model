# P2.3 — Exact Donor Integration (B0-X / B1-X): Result

**Record date:** 2026-10-22
**Question:** how much of the certified B0→B1 neutral-area increase survives
removal of the centroid whole-element donor-classification error?

## 1. Design

Same meshes as B0 and B1. Only the donor-load assembly changes:

```
f_i^X = rho * int_{T cap {x >= 0.54}} N_i(x,y) dA
```

Everything else — geometry, stiffness, Bernstein constraints, BCs, solver,
acceptance — is untouched. B0/B1 centroid solutions are preserved as the
immutable baseline.

## 2. Assembly validation

Implemented in the parent barycentric frame: clip to `x >= 0.54`, fan
triangulate, 3-point degree-2 exact rule, sub-area from the barycentric
determinant times parent area.

| Test | B0-X | B1-X |
|---|---:|---:|
| Full-element identity `[0,0,0,A/3,A/3,A/3]`, max err | 1.27e-21 | 1.69e-21 |
| `Σf/ρ` (µm²) | 1.575269192653 | 1.575269192653 |
| Exact clipped charged area (µm²) | 1.575269192653 | 1.575269192653 |

`Σf/ρ` equals the exact charged area to all twelve digits **and is identical
across the two meshes**, unlike the centroid totals (1.575255252958 vs
1.575498331328).

**Candidate-set check** (the open item in the session record): the centroid mask
is the constraint selection, not a donor filter.

```
B0:  charged(constraints) 50,612   exact-load candidates 50,702   gained  90   lost 0
B1:  charged(constraints) 59,203   exact-load candidates 59,264   gained  61   lost 0
```

Every triangle with a non-zero clipped charged part contributes; the extra
elements are those whose centroid lies below the cut but which extend above it.
`lost = 0` confirms no contribution is dropped.

## 3. Certification

```
B0-X: stat 5.834e-11  comp 4.166e-17  slack -4.873e-15  viol 0
      oxide 4.888e-13  overshoot 2.983e-15  active 4,435    GATE PASS
B1-X: stat 8.837e-11  comp 9.952e-18  slack -2.853e-15  viol 0
      oxide 4.405e-13  overshoot 1.511e-15  active 16,719   GATE PASS
```

OSQP: B0-X `solved` 5,025 iters, obj `-1.052111e+05`; B1-X `solved` 5,350 iters,
obj `-9.637236e+04`. Both checkpoints saved and reload-verified.

## 4. Bounded areas

| Formulation | B0 | B1 |
|---|---|---|
| Centroid | [0.021644935001, 0.021649266844] | [0.021987652318, 0.021991846424] |
| Exact | [0.021644934927, 0.021649266770] | [0.021987627914, 0.021991822020] |

## 5. The decisive comparison

```
ΔA_centroid  [+3.383854740e-04, +3.469114230e-04]   width 8.53e-06
ΔA_exact     [+3.383611440e-04, +3.468870930e-04]   width 8.53e-06

correction   ΔA_X − ΔA_C = [-8.550e-06, +8.502e-06]     straddles zero
             = [-2.46 %, +2.51 %] of the increment
```

Per-mesh source corrections:

```
δA_B0 = A_B0X − A_B0  [-4.332e-06, +4.332e-06]   straddles zero
δA_B1 = A_B1X − A_B1  [-4.219e-06, +4.170e-06]   straddles zero
```

Both per-mesh intervals straddle zero and are no wider than the enclosure width
itself (4.33e-6).

## 6. Interpretation

Against the pre-registered outcome table:

| Outcome | This result |
|---|---|
| ΔA_X ≈ ΔA_C, strictly positive | **yes** |
| ΔA_X materially different, still positive | no |
| ΔA_X encloses zero or flips sign | no |

**The centroid-source approximation is not material to the B0→B1 increment.**
The band-enabled increase survives exact donor integration essentially intact:
the correction is bounded at ±2.5 % of the increment and is not resolvable at
the current enclosure width.

This also answers the earlier caution that a ~100× ratio between relative
source and area changes does not bound the area response: the obstacle problem
*was* given the chance to amplify the source perturbation, and it did not. The
bound is now empirical, not dimensional.

## 7. Limits

- ±2.5 % is a bound, not a measurement of zero. Tighter enclosures
  (`AREA_TOL = 1e-13`) would be needed to resolve a smaller correction.
- Exact donor integration removes centroid misclassification only. Stiffness
  discretisation, active-set response and mesh arrangement still differ between
  B0-X and B1-X, so this is not a pure resolution effect.
- No continuum limit and no convergence order follow.
- B2 remains unauthorised.

## 8. Artifacts

- `p2_exact_donor.py`, `p2_B0X_*`, `p2_B1X_*`
- `p2_B0X_post.json`, `p2_B1X_post.json`
- `p2_B0X_scaled.json`, `p2_B1X_scaled.json`
