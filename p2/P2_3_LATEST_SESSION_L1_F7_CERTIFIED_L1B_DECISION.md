# P2.3 — Latest Session: L1 F7 Certification, Spatial Nonconvergence, and L1b Decision

**Record date:** 2026-10-20  
**Status:** L1 **F7-certified algebraically**; **P2.3 still OPEN** because spatial convergence is not established. P1 frozen; Stage 3 suspended.  
**Next approved action:** Build and solve **L1b**, with `H_LOCAL = 0.0055 µm` in the same two windows, unchanged physical model, and mandatory full-set F7 certification before postprocessing conclusions.  
**Evidence basis:** User-reported solver, mesh, and postprocessing outputs in this conversation. The production L1 solve and diagnostics were **not independently rerun** while preparing this document.

---

## 1. Executive outcome

The L1 run has now passed the complete algebraic F7 gate. The successful method was **checkpoint-first OSQP**, then **row-scaled, matrix-free MINRES KKT polishing**, followed by a **residual-correction refinement**. The latter reduced relative stationarity from `3.571e-07` to `1.484e-10`; no large direct SuperLU factorization was required.

However, L1's corrected bounded neutral area increased by **[0.000277430305, 0.000285823012] µm²** relative to the intermediate mesh. That increment is demonstrably larger than the preceding r4-to-intermediate increment **[0.000211061681, 0.000219491826] µm²**, by a factor bounded by **[1.264, 1.354]**. The neutral front also moves nonmonotonically at x = 2.35 µm, and the downstream shallow interval structure changes. Thus:

- **F7 / algebraic KKT:** PASS for L1.
- **L1 spatial-convergence evidence:** NOT PASS / insufficient.
- **Interpretation:** mesh-distribution sensitivity has been demonstrated; **refinement contrast is a hypothesis, not yet a proven cause**.
- **Decision:** Perform one controlled lower-contrast local-mesh experiment (L1b). Do **not** run a global r5 or Stage 3.

## 2. Frozen model and experiment definition

The physical geometry, material assignments, BCs, donor support, and P2 Bernstein obstacle remain unchanged from the matched-model P2.3 baseline. The model uses a conforming silicon/oxide mesh; six-node P2 triangles; the silicon obstacle `u = V_D - phi >= 0` imposed through nonnegative Bernstein coefficients; and a QP solved by OSQP and polished with a rank-aware KKT procedure. A successful polished solve must pass F7 against **all** inequality constraints, not only an active subset.

Relevant baseline mesh: intermediate global target `h = 0.0085 µm` (tip size reported as `0.0034 µm`). L1 introduced two locally finer regions against approximately the same global background:

- **Window A (upstream):** `2.32 <= x <= 2.37 µm`, `|y| <= 0.13 µm`.
- **Window B (downstream):** `2.395 <= x <= 2.420 µm`, `|y| <= 0.24 µm`.
- **L1 local target:** `H_LOCAL = 0.0040 µm`.
- **Approved L1b local target:** `H_LOCAL = 0.0055 µm` with the **same windows**, graded transitions, and unchanged physics.

## 3. L1 mesh — built and geometrically healthy

| Metric | L1 reported |
|---|---:|
| Mesh nodes | 451,894 |
| Mesh elements | 225,245 |
| Elements above intermediate | 23,627 |
| Si elements | 171,407 |
| SiO2 elements | 31,728 |
| Gate elements | 6,448 |
| Shield elements | 15,662 |
| Minimum triangle area | `2.639e-06` (mesh area units) |
| Maximum aspect ratio | 2.55 |
| Minimum interior angle | 36.75° |
| Si area | 1.954544961 µm² |
| Charged Si area | 1.575290811 µm² |

**Scope of health claim:** The reported geometric measures are favorable. They do not, by themselves, prove smoothness of the *actual* element-size field at the boundaries of the refinement windows or equivalence of all source/BC tags. L1b must retain explicit mesh/assembly checks.

## 4. L1 QP solution and resource history

| Quantity | L1 reported |
|---|---:|
| P2 DOFs | 451,448 |
| Free DOFs | 375,812 |
| Bernstein constraints | 939,738 |
| Charged Si elements | 156,623 |
| Body elements | 14,784 |
| OSQP status | `solved` |
| OSQP iterations | 15,350 |
| OSQP objective | `-1.024806e+05` |
| OSQP primal residual | `1.386e-08` |
| OSQP dual residual | `1.638e-09` |

### 4.1 Failed direct polish

The initial rank-reduced KKT method attempted to factorize the free 375,812-DOF operator with SuperLU. The factorization consumed approximately **81 CPU minutes** and **5.9 GB** without reaching the active-set stage. It was stopped as a **resource failure**, not a failed numerical solution. No valid L1 spatial conclusion was drawn from that run.

### 4.2 Checkpoint recovery and API issue

An initial save-first rerun completed OSQP but failed while writing solver status because the installed OSQP build used `res.info.prim_res` rather than `res.info.pri_res`. A version-tolerant field lookup (`prim_res` with `pri_res` fallback) fixed the script. A subsequent complete solve wrote and reloaded:

- `p2_L1_osqp_checkpoint.npz` — reported ~4.9 MB; contains `x(375812,)`, `y(939738,)`, `idx(375812,)`.
- `p2_L1_osqp.json` — status and metadata.

A separate wait-loop process confirmed that the saved checkpoint was readable. **Important:** These filenames refer to outputs on the user's compute environment; they are not claimed as attachments to this Markdown file.

### 4.3 Iterative polish and successful refinement

The raw OSQP iterate did not meet F7 stationarity: approximately `6.42e-05`. An initial matrix-free MINRES polish reduced this to approximately `4.05e-07`, still above the `1e-08` stationarity threshold. The later scaled/refined KKT solve reported:

```text
scaled + refined MINRES:
  refine 0: stat_rel = 3.571e-07
  refine 1: stat_rel = 1.484e-10
```

The refinement improvement factor was approximately **2,406×** for these two reported stationarity values. The crucial improvement was the **residual-correction round**, not diagonal scaling by itself. The procedure avoided the unsuccessful large SuperLU factorization.

## 5. L1 final algebraic F7 gate — PASS

```text
stationarity      1.484e-10
complementarity   5.778e-14
minimum slack    -4.222e-11
full violations   0
minimum lambda    0
oxide residual    4.766e-13
overshoot         2.127e-11
F7 GATE           PASS
```

The minimum slack is slightly negative in floating-point arithmetic and is reported transparently; it is **within the applied F7 tolerance**, rather than an assertion of exact mathematical feasibility. The zero-violation result refers to the gate's defined tolerance. L1 is certified **as a discrete algebraic obstacle solution**, not as a converged continuum field or front.

**Pipeline change to retain:** Save the OSQP primal/dual checkpoint *before* attempting any polish; prefer matrix-free row-scaled iterative polishing with explicit residual correction; then recompute F7 on the original unscaled full system, including all Bernstein constraints, oxide residual, and physical overshoot.

## 6. F2 correction history — retractions and current area bounds

A prior 'zero-gap exact Bernstein area' report was invalidated. Its recursion evaluated subtriangle Bernstein coefficients using **parent-frame nodes**, a coordinate-frame bug. Accordingly, the earlier apparent exact areas

```text
r3   0.021318701188
r4   0.021851690577
mid  0.022032403652
```

are **RETRACTED** and **MUST NOT** be used in any comparison, graph, fit, or future certification. Also retracted: the claim that P2/P1 disagreement vanished to +0.03% at r4. The genuine approximate **−1.3% P2/P1 offset** remains.

The corrected Bernstein subdivision implementation produced **nonzero enclosure widths**. The linear manufactured test `u = L1`, `eps = 0.5`, on a triangle of area 0.5 gave **[0.375, 0.375 + 1.907e-06]**, containing the exact area 0.375 rather than falsely collapsing to zero gap.

### 6.1 Current accepted numerical area intervals

Units below are µm²; values are **bounds**, not exact areas.

| Mesh | Lower | Upper | Width |
|---|---:|---:|---:|
| r3 | 0.020926837673 | 0.020930251232 | 0.000003413559 |
| r4 | 0.021574089791 | 0.021578224392 | 0.000004134601 |
| Intermediate | 0.021789286073 | 0.021793581617 | 0.000004295544 |
| **L1** | **0.022071011922** | **0.022075109085** | **0.000004097163** |

The corrected recursion used reported subdivision setting `AREA_TOL = 1e-11 µm²`; this must **not** be confused with a guaranteed total area enclosure width of `1e-11 µm²`. The global widths above are approximately `3.4e-06`–`4.3e-06 µm²`. The previously proposed `<= 1e-06 µm²` final-width target is **not yet attained**. A tighter `AREA_TOL = 1e-13` evaluation remains a separate final-certification task, not a prerequisite for L1b.

P1 exact clipping comparator remains:

- r3 `0.021196729 µm²`.
- r4 `0.021845237 µm²`.
- Corrected P2/P1 range roughly `0.9873–0.9874` at r3 and `0.9876–0.9878` at r4, about **−1.2% to −1.3%**.

The earlier equal-area centroid quadrature values also remain **superseded as area estimators**; they are not rigorous bounds.

## 7. Interval propagation: the L1 area increment is larger

For area enclosures `A_i in [L_i,U_i]`, propagate a signed difference without replacing it by midpoint subtraction:

`A_b - A_a in [L_b - U_a, U_b - L_a]`.

| Comparison | Lower increase (µm²) | Upper increase (µm²) |
|---|---:|---:|
| r3 → r4 | 0.000643838559 | 0.000651386719 |
| r4 → intermediate | 0.000211061681 | 0.000219491826 |
| **Intermediate → L1** | **0.000277430305** | **0.000285823012** |

The first two steps give a tightly bounded reduction ratio (later increment / prior increment) approximately **[0.324, 0.341]**, consistent with a roughly threefold decrease. The next, heterogeneous local-refinement step gives:

```text
(min L1 increment) / (max preceding increment) = 1.264
(max L1 increment) / (min preceding increment) = 1.354
```

Thus the **L1 area increment increased by a factor between 1.26 and 1.35**, even after rigorous propagation of the reported area intervals. This is a real change in the measured discrete quantity at the current area-bound resolution, **not** something that could be explained by its ~4e-06 µm² enclosure width alone.

**Interpretive limitation:** r3, r4, intermediate, and L1 are **not a uniform one-parameter family of nested globally refined meshes**. L1 changes the *spatial distribution* of resolution. The data challenge spatial convergence and reveal mesh-design sensitivity; they do not prove divergence of the underlying FEM sequence or establish the mechanism.

## 8. F3 front-cut comparisons: intermediate → L1

All entries in µm; reported threshold `eps = 1e-8 V`. Parenthesized pairs are sublevel intervals along the vertical cut, not uncertainty intervals.

| x (µm) | Intermediate intervals | L1 intervals | Key change |
|---|---|---|---|
| 2.35 | `(-0.0754, +0.0686)` | `(-0.0718, +0.0684)` and `(+0.0722, +0.0763)` | Lower edge moves inward by 0.0036 µm; small extra interval appears |
| 2.38 | `(-0.1314, +0.1331)` | `(-0.1397, +0.1398)` | Cut widens on both sides (~0.0083 / 0.0067 µm) |
| 2.39 | `(-0.1564, +0.1593)` | `(-0.1541, +0.1595)` | Smaller asymmetric displacement |
| 2.40 | `(-0.1816, +0.1746)` | `(-0.1773, +0.1783)` | Both edges move inward/outward differently |
| 2.41 | Three reported intervals | `(-0.1996, +0.1943)` and `(+0.1964, +0.1982)` | Still a shallow multiple-interval configuration |

### 8.1 Upstream front nonmonotonicity

At `x = 2.35 µm`, the lower endpoint across four meshes is:

```text
r3            -0.0624 µm
r4            -0.0697 µm
intermediate  -0.0754 µm
L1            -0.0718 µm
```

The sequence drifts outward through intermediate, then **reverses by +0.0036 µm** for L1. Therefore the front cannot be called stabilized or monotone with respect to this mixed refinement sequence.

### 8.2 Downstream shallow intervals

At `x = 2.41 µm`, L1 shows two intervals at `eps = 1e-8 V`, which **merge at `eps = 1e-5 V`**. This threshold sensitivity is not evidence of a robust isolated neutral pocket or a topological bifurcation in the continuum solution. Threshold sweeps and careful interval stitching/deduplication remain required.

### 8.3 What is *not* yet known

No credible limiting front position, rigorous continuum neutral area, tip-field convergence, or isolated-pocket topology follows from these L1 cuts. Algebraic KKT accuracy does not supply a FEM discretization error estimate.

## 9. Diagnosing refinement contrast: hypothesis, not proof

One candidate mechanism is the approximately `0.0040 µm` local target embedded within a `0.0085 µm` background, with resolution transition near the neutral boundary. Such a change can shift the discrete active-set geometry and front representation. But the existing measurements do **not** isolate this cause from local resolution changes, transition placement, or intrinsic sensitivity of the discrete sharp-depletion obstacle.

**L1b tests a specific hypothesis:** A more moderate `H_LOCAL = 0.0055 µm` will reduce sensitivity to the local/global size contrast. It is a diagnostic comparison, not an optimization to force a smaller area increment.

If L1b stabilizes the front, that *supports* mesh-distribution involvement but still does **not prove** the contrast mechanism, because the central-window resolution has also changed. A future controlled test should vary grading/transition placement while holding the central resolution fixed.

## 10. Approved L1b execution plan

### 10.1 Build and input freeze

1. Create a **new** conforming L1b mesh with exactly the same geometry, material regions, interfaces, electrode surfaces, donor-support exclusions, BCs, and two refinement windows.
2. Set `H_LOCAL = 0.0055 µm`; retain the `0.0085 µm` nominal background. Use smooth grading at the window boundaries; record **actual triangle sizes and transition gradients**, not only Gmsh target values.
3. Check physical material/charged/body areas, element counts, positive orientation/Jacobians, element angles/aspects, and interface conformity; inspect any geometry changes around the tip/cap.
4. Verify assembly consistency, especially charged-element tags, body exclusion, fixed DOFs, and P2 Bernstein coefficient mapping.

### 10.2 Solve and certify

5. Run the existing P2 Bernstein-obstacle QP with OSQP.
6. **Immediately** save primal and dual iterates, mapping, solver status and residuals, objective, model/mesh fingerprints. Reload the checkpoint as a durability check.
7. Compute raw F7 diagnostics. If necessary, run row-scaled matrix-free MINRES polishing with the proven residual-correction round and controlled iteration/memory limits.
8. Apply the **original full-set F7 gate** in unscaled coordinates: stationarity (target `<=1e-8`), feasibility, multiplier sign, complementarity, oxide residual, and physical overshoot. **No certified spatial comparison unless the full gate passes.**

### 10.3 Postprocess only after F7

9. Compute **Bernstein-subdivision neutral-area lower and upper bounds** using the fixed child-coordinate transformation; report total enclosure width and its propagation in comparisons with intermediate and L1.
10. Extract complete **quadratic-aware analytic front-cut intervals**, including slivers, at x = `2.35, 2.38, 2.39, 2.40, 2.41 µm` and thresholds `eps = 1e-8, 1e-6, 1e-5 V`.
11. Report neutral-region element edge-connectivity and its limitations; interface normal-flux mismatch (wall/arc/cap); symmetry measures; and **fixed-physical-location** electric-field probes.
12. Archive solver logs, OSQP checkpoint, polished solution, F7 gate report, mesh-quality and input/tag checks, area bounds, front intervals, and a single side-by-side comparison table.

### 10.4 Required cross-mesh comparison

| Measurement | Intermediate | L1 (certified) | L1b (pending) |
|---|---|---|---|
| `H_LOCAL` | n/a (global baseline) | 0.0040 µm | 0.0055 µm |
| F7 | PASS | **PASS** | Required |
| Neutral-area bracket (µm²) | `[0.021789286073, 0.021793581617]` | `[0.022071011922, 0.022075109085]` | To compute |
| Increment vs intermediate (µm²) | — | `[0.000277430305, 0.000285823012]` | To bound |
| x=2.35 lower edge at 1e-8 V (µm) | -0.0754 | -0.0718 | To extract |
| x=2.41 topology at 1e-8 V | Three intervals | Two intervals | To extract |
| x=2.41 behavior at 1e-5 V | Threshold-sensitive | L1 intervals merge | To test |
| Interface flux, symmetry, fixed-probe E | Baseline diagnostics | Compare retained records | To report |

## 11. L1b interpretation and decision tree

| Observation | Interpretation justified | Follow-up |
|---|---|---|
| Area movement decreases and front moves toward intermediate | Mesh-distribution sensitivity supported | Second independently graded mesh before claiming convergence |
| L1b nearly reproduces L1 | Sensitivity persists under lower contrast | Examine active-set/front and transition location |
| Area comparison improves but front still oscillates | Global area hides local geometric instability | Prioritize cut/root and active-set diagnostics |
| Area increases again, or cut structure changes again | Existing spatial convergence claims remain unsupported | Change experimental design; no Stage 3 |
| F7 fails | L1b spatial metrics cannot be accepted as certified | Repair solver/polish before interpreting metrics |

Do not label a front reversal as a definitive 'mesh-contrast artifact' based solely on L1b; the experiment varies contrast and absolute local resolution together.

## 12. Remaining audit and project gates

- **F7:** CLOSED for r3, r4, intermediate, and now **L1**. Re-run on every newly generated mesh.
- **F9:** CLOSED for the previously tested fully neutral *element* edge-connectivity. Repeat the diagnostic on L1b; it does not certify topology of the full continuous P2 sublevel set.
- **F2:** Corrected area-bound method operational; erroneous exact numbers explicitly retracted. Final proposed global width `<=1e-6 µm²` is still open; tighter `AREA_TOL = 1e-13` evaluation planned **after** controlled local experiments, prior to final 'certified area' claims.
- **F3:** Analytic local front diagnostics available; *spatial front convergence remains OPEN*.
- **F4:** Make maps quadratic-aware (full P2 interpolation, not vertex-only linearized contouring).
- **F5:** Manufactured dielectric-interface patch test must exercise production `interface_flux`.
- **F11:** Verify streamline electric-field vectors, interpolation, and orientation.
- **F6/F8 / legacy:** Quarantine obsolete `p2_build.py` assembly and `conforming_fem.solve_complementarity` paths.
- **P1:** Frozen; separate independent archive/provenance audit still pending.
- **P2.3:** OPEN; no established continuum neutral-area limit, resolved front geometry, or tip-field convergence.
- **Global r5:** NOT AUTHORIZED by current decision.
- **Stage 3 tip-radius sweep:** SUSPENDED pending P2.3 convergence.

## 13. Explicit final decision

> **GO — L1b at `H_LOCAL = 0.0055 µm` in unchanged windows A and B, with unchanged physics.** Preserve the algebraically certified L1 result and its corrected bounded-area/front-cut diagnostics as a reference. Apply checkpoint-first OSQP and the scaled matrix-free iterative KKT polish; require a full-set F7 PASS before neutral-area or front-cut conclusions. Do not run global r5 or Stage 3.

The immediate objective is to distinguish **mesh-distribution dependence** from a robust continuum trend. It is **not** to fit a convergence order or report a final neutral area from the current four heterogeneous meshes.

---

### Record-integrity notes

1. This session document supersedes any previous session record stating the retracted zero-gap Bernstein values as valid or calling F2 'exact'.
2. `P2_3_CONVERGENCE.md` §10 was reported by the user to contain the corrected interval/retraction; this document does not claim to have edited that file.
3. Results for L1b are **pending**; none have been fabricated or extrapolated.
4. All numeric mesh, QP, KKT, area, and front-cut data above are transcribed from user-reported runs. Difference intervals and comparison ratios were arithmetically recomputed from the supplied endpoints.
