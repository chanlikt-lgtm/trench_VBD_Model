# P2.3-S0 — Freeze Record: Completed Controlled-Refinement Experiment

**Status:** CLOSED and FROZEN as a controlled-refinement experiment.
**Accepted at:** `ebff654` (CP-L2-A), record consolidated after `65746af`.
**Scope of the freeze:** the three-level S0 family, its certified solutions, area
enclosures, increment comparison and arc-region diagnostic. Nothing here is a
continuum result.

## 1. Certified numerical record

| Level | Neutral-area enclosure (µm²) | Width | F7 |
|---|---|---:|---|
| S0-L0 | [0.022074881123311834, 0.022075207649363562] | 3.265e-07 | PASS (`a3a5248`) |
| S0-L1 | [0.022217431654472784, 0.022217790436636530] | 3.588e-07 | PASS (`e1c8a2d`) |
| S0-L2 | [0.022293016138576278, 0.022293346039891330] | 3.299e-07 | PASS (`6e8d9e1`) |

Increments, endpoint subtraction only:

```
ΔA_01^S0  [1.422240051e-04, 1.429093133e-04]   strictly positive
ΔA_12^S0  [7.522570194e-05, 7.591438542e-05]   strictly positive

sup|ΔA_12| = 7.591438542e-05  <  1.422240051e-04 = inf|ΔA_01|
ratio       [0.5264, 0.5338]      (second increment ≈ 53 % of the first)
candidate p [1.81, 1.85]           descriptive only — NOT a certified order
```

The shrinkage is certified by non-overlapping intervals, not by comparing
midpoints.

## 2. Arc-region diagnostic (fixed region, `AREA_TOL = 1e-13`)

```
            arc lower    arc upper      complement
S0-L0       0            1.362955e-07   [0.022074744828, 0.022075207649]
S0-L1       0            1.512802e-07   [0.022217280374, 0.022217790437]
S0-L2       0            1.325195e-07   [0.022292883619, 0.022293346040]

ΔA_arc,01  [-1.362955e-07, +1.512802e-07]   contains zero
ΔA_arc,12  [-1.512802e-07, +1.325195e-07]   contains zero
```

Direct arc contribution ≤ 1.513e-07 µm² = **0.106 % of ΔA_01** and
**0.201 % of ΔA_12**. Arc upper bounds are non-monotone across levels because
they are unresolved straddler ceilings, and the lower bound is zero at every
level because no sub-triangle in the band is entirely neutral. These are
ceilings, not measurements. Regional quantities are **DIAGNOSTIC ONLY**.

## 3. Error-budget ledger (what is closed, what is not)

| Component | Status |
|---|---|
| Algebraic KKT / obstacle certification at each level | **CLOSED** — full-set F7 at L0, L1, L2 |
| Controlled global refinement with one frozen polygon | **CLOSED** — boundary identity 1.110e-16, scaling ratios 0.700–0.721 |
| Diminishing neutral-area sensitivity across tested levels | **CERTIFIED** for these three levels |
| Direct arc-region contribution to the increments | **BOUNDED** ≤ 0.2 % |
| Nonlocal influence of the arc refinement shortfall | **OPEN** — elliptic coupling and active-set response not bounded |
| Fixed-polygon vs curved-domain geometry error | **OPEN** — sagitta 7.53e-06 µm, canonical area −1.631e-05 µm² vs the curved reference |
| Asymptotic regime / convergence order | **NOT ESTABLISHED** — two increments only |
| Finite-element topological nesting | **NOT VERIFIED** |

The consequence is stated once, plainly: a decrease in successive increments on
a fixed polygon is evidence of diminishing sensitivity, not of convergence to a
curved-domain continuum limit.

## 4. Frozen artifacts

```
p2_3_spatial_S0/
  levels/S0_L{0,1,2}_p2mesh.npz          meshes (hashes in CP-M2-A_manifest.json)
  CP-M2-C_L0_polished.npz                certified L0 field
  CP-L1_L1_polished.npz                  certified L1 field
  CP-L2_L2_polished.npz                  certified L2 field
  A_S0-L{0,1,2}.json                     area enclosures
  DA_S0_01.json, DA_S0_12.json           increments
  CP-L2-A6_arc_region.json, _summary.json  arc diagnostic
  P2_3_S0_Results_Tables.xlsx            paper-ready tables Q1–Q6 + notes
```

Commits: `e937926` design · `1c8d4df` repaired pilots · `a5b76e2` manifest ·
`a3a5248` L0 F7 · `76f6f42` L0 area · `e1c8a2d` L1 F7 · `803dcff` L1 area ·
`c07ec2f` L2 pilot · `6e8d9e1` L2 F7 · `5773158` L2 area · `ebff654` arc
diagnostic · `65746af` workbook.

No certified artifact may be modified. A re-solve requires a new authorization.

## 5. What this experiment does and does not conclude

> S0 neutral area increases with refinement, successive increments decrease by a
> rigorously resolved amount, and the direct arc-region contribution is small.
> Continuum convergence remains unestablished.

## 6. Disposition

```
S0-L0/L1/L2 F7 certification            CLOSED
S0 area enclosures and increment test   CLOSED
CP-L2-A6 arc diagnostic                 CLOSED
Diminishing area sensitivity            CERTIFIED for tested levels
Direct arc-area contribution            QUANTITATIVELY BOUNDED
Nonlocal arc-refinement effect          OPEN
Fixed-polygon vs curved-domain error    OPEN
S0-L3                                   NOT AUTHORIZED
Chord-policy modification               NOT AUTHORIZED
Richardson / convergence-order claim    NOT AUTHORIZED
P2.3 continuum convergence              OPEN
P1                                      FROZEN
Stage 3                                 SUSPENDED
```

**Next work, if any, must be designed around the two open error terms** — the
nonlocal influence of arc-region resolution, and the fixed-polygon geometry
error — rather than by extending S0 to a fourth level.
