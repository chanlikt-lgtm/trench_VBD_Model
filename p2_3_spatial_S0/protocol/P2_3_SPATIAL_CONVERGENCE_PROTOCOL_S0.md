# P2.3-S0 — Controlled Spatial-Convergence Study: Design Protocol and Checkpoint Gates

**Document status:** PROPOSED / DESIGN-ONLY — awaiting review and approval  
**Scope:** Plan a new P2 conforming FEM spatial-convergence experiment; **no mesh generation or production solve is authorized by this document**.  
**Baseline provenance:** Closed B0-X/B1-X/B2-X band-sensitivity experiment (`15ec811` → `2478f82` → `4e5e6d4`).  
**Disposition to preserve:** P2.3 continuum convergence **OPEN**; P1 **FROZEN**; Stage 3 **SUSPENDED**; any third level of the **old** band experiment **NOT AUTHORIZED**.  
**Prepared:** 2026-10-23  
**Proposed record name:** `P2_3_SPATIAL_CONVERGENCE_PROTOCOL_S0_CHECKPOINTS.md`

> **Evidence provenance.** Values identified below as *archived observations* are drawn from results supplied in the project conversation; the production code, commits and raw files have **not** been inspected in this document. Numerical targets and checkpoint thresholds labeled *proposed* are design criteria, not already verified facts.

---

## 1. Purpose and decision boundary

The closed B0-X/B1-X/B2-X study established **diminishing bounded neutral-area sensitivity along that tested sequence**. It did **not** produce a controlled asymptotic refinement family: B0-X has no active boundary-size reduction, grading differs between levels, and the active set can vary. An apparent exponent near 2.1 from that sequence **must not** be used as a FEM convergence order.

The next study should isolate *spatial discretization* as far as practicable by freezing the physical problem, exact source integration, interface/geometry policy, P2 obstacle formulation, neutral-area estimator, and a single mesh-size-field family. The first milestone is **S0: a reviewed, reproducible protocol and resource model**, not another large solve.

**Design gate:** approve S0 design only. Further, separate authorizations are required for (i) pilot mesh generation, (ii) production mesh-family generation, and (iii) solution and area audits. Approval of a protocol is **not** implicit approval of these actions.

### Required S0 deliverables

1. An unambiguous master geometry and mesh-size-field definition, with its version/hash.
2. A choice between fixed-geometry (discrete-domain) convergence and geometry-refined physical-domain convergence, with a budget for geometry error.
3. A proposed level table plus a resource-feasibility forecast with stated uncertainties.
4. Frozen assembly, obstacle, solver-polish, F7, area-interval, and diagnostic policies.
5. A checkpointed execution plan including archive/rollback policy, quantitative pass/hold/fail conditions, and an explicit stopping rule.
6. A signed decision note specifying exactly which later actions, if any, become authorized.

---

## 2. Established baseline: do not re-interpret as Richardson data

### 2.1 Certified exact-source neutral-area enclosures (µm²)

| Solution | Lower bound | Upper bound | Width |
|---|---:|---:|---:|
| B0-X | 0.021645551625654 | 0.021645845589193 | 2.940e-07 |
| B1-X | 0.021988322072981 | 0.021988665490453 | 3.434e-07 |
| B2-X | 0.022154236625040 | 0.022154596319621 | 3.597e-07 |

Certified increments, by conservative interval subtraction:

\[
\Delta A_{01}^{X}\in[3.424764838,3.431138648]\times10^{-4}\;\mu\mathrm m^2,\qquad
\Delta A_{12}^{X}\in[1.655711346,1.662742466]\times10^{-4}\;\mu\mathrm m^2.
\]

Thus

\[
\sup(\Delta A_{12}^{X})=1.662742466\times10^{-4}
<3.424764838\times10^{-4}=\inf(\Delta A_{01}^{X}).
\]

The ratio enclosure reported for the closed sequence is **[0.4826, 0.4855]**. This is a certified reduction in successive area-increment magnitude **for that sequence only**.

### 2.2 B2-X archived F7 evidence

| Diagnostic | Archived result |
|---|---:|
| Active constraint rows | 31,723 |
| Stationarity | 1.0195393954383166e-10 |
| Complementarity | 8.082037731392984e-18 |
| Minimum Bernstein slack | -3.2276814617319727e-15 |
| Violations | 0 |
| Minimum multiplier | 0.0 |
| Oxide residual | 4.5094412746102277e-13 |
| Overshoot | 1.6785376981639928e-15 |
| MINRES info | 0 |
| F7 gate | PASS |

Other archived provenance: exact donor integral \(\sum_i f_i/\rho=1.575269192653\,\mu\mathrm m^2\) on the tested meshes; charged-support audit confirms nonzero clipped portions contribute even when a triangle centroid lies on the body side. The tightly bounded *difference of band increments caused by donor integration* on B0/B1 was within approximately **±0.19%** of the first increment, with unresolved sign. This is **not** an independently estimated continuum error.

### 2.3 Unresolved scientific question

Does the corrected P2 sharp-depletion solution's neutral area (and, separately, front geometry) stabilize under a **controlled refinement family** for the same physical model? A positive answer needs several resolved differences and a defensible geometry/source/numerical error budget. F7 certification at any individual mesh does not answer it.

---

## 3. Freeze the mathematical model and discrete invariants

All entries in this section require a concrete file path, git commit, config digest, and audit output at **CP-S0.2**; listing them here alone does not establish their provenance.

| Item | Frozen requirement | Evidence to archive |
|---|---|---|
| Physics | 2D electrostatics, same silicon/SiO₂ regions, electrodes, drain bias, symmetry conditions, donor density and depletion convention | Geometry/BC/material manifest |
| Drain voltage | 34.274208536 V (archived model parameter) | Exact configuration value |
| Donor density | \(N_D=6.9\times10^{16}\,\mathrm{cm}^{-3}\) | Physical-units manifest |
| Source scaling | \(\rho=qN_D/\epsilon_0\approx1248.5638\,\mathrm{V/\mu m^2}\) (archived) | Unit test; constant version |
| Body cutoff | \(x=0.54\,\mu\mathrm m\), same physical clipping predicate and boundary convention | Clip oracle fixtures |
| FE space | Six-node conforming P2 triangles, quadratic basis and the existing Dirichlet elimination | Local reference-element tests |
| Stiffness | Same material tags, dielectric coefficients, quadrature and numerical assembly conventions | Assembly checksums / patch tests |
| Donor load | **Exact partial-element P2 integration** over \(T\cap\{x\ge0.54\}\), including all partially charged triangles | Local identity and global charge checks |
| Obstacle | \(u=V_D-\phi\ge0\) in silicon via the same Bernstein inequalities | Constraint builder/hash |
| Solver | Checkpoint-first OSQP; existing rank-reduced scaled MINRES/KKT polish; corrected residual sign | Solver options manifest |
| F7 | Full-set nonnegative Bernstein slacks and active multipliers, zero inactive multipliers, complementarity, stationarity, oxide check, overshoot check | Full diagnostic JSON |
| Area | Corrected, coordinate-frame-safe Bernstein branch-and-bound estimator, outward intervals | Area algorithm/tests/version |
| Front | Full quadratic interpolation, correct barycentric arguments, fixed x-stations and epsilon sensitivity protocol | Diagnostics config / tests |

**No silent changes:** a change in geometry, quadrature, solver tolerances, constraint normalization, area definition, or front extraction invalidates the simple comparison unless explicitly evaluated as an independent control.

### 3.1 Geometry-control fork (must be resolved in S0)

**Option G-fixed (recommended first):** All levels solve on an identical *polygonal* domain with identical physical/material interfaces and body-cut representation; only the FE interior triangulation is refined. This gives a clean *discrete-domain* mesh study, **not yet** an error estimate relative to the intended curved physical domain. Boundary vertices/segments may be subdivided **only along the same polygonal segments**, without moving the polygonal domain.

**Option G-physical:** Refine a curved/geometry-aware boundary approximation in a controlled sequence, and independently bound or measure geometry approximation error. This targets the physical domain more directly but introduces another varying error source.

**Checkpoint rule:** Record the chosen target. Never mix G-fixed and G-physical meshes in a single Richardson sequence. If G-fixed is chosen, explicitly label all area-limit claims as limits for the frozen polygonal geometry until a separate geometry study closes the physical-domain gap.

### 3.2 Body-source representation under refinement

The physical clipped charged region changes only if the geometry changes. On G-fixed, the clipped charged area must be invariant to numerical tolerance across levels. The discrete *load vector* will differ with mesh even though its integrated total matches the same area. **Charge conservation does not establish load-vector identity or eliminate spatial source discretization.** Retain the exact cut-element scheme on every level.

---

## 4. Freeze a single mesh-size-field family

Define a continuous/nonnegative master target size function \(h_0(x,y)\) over the full frozen domain, including the tip, frozen neutral-boundary band, inter-window strip, grading transition, background, oxide, and interfaces. Then prescribe

\[
 h_k(x,y)=r^k h_0(x,y),\qquad r=2^{-1/2}\approx0.70710678.
\]

**Important:** The global scaling rule is the preferred design for a controlled family; if a fixed-size tip cap, boundary cap, or unscaled background is retained, the resulting study is no longer globally similar and must be explicitly called a *localized-refinement family*. Do not quietly claim a global Richardson ratio. An alternative fixed-band-only study is possible, but would answer a narrower question and must carry its own error budget for the unrefined remainder.

### 4.1 Frozen band geometry

- Reuse the archived L1b-derived neutral boundary at \(\epsilon=10^{-8}\,\mathrm V\), if and only if its exact coordinates and provenance can be recovered and hashed.
- Archived extraction: a connected component of length approximately 1.2492 µm, with x-range about 2.333–2.420 µm and y-range approximately ±0.575 µm, resampled to 625 points at ~0.002 µm spacing.
- Band core: distance ≤0.04 µm from the fixed curve; outward transition width approximately 0.02 µm. **Keep the physical widths constant across levels.**
- Ensure that the inter-window x strip [2.35, 2.40] µm is explicitly covered in all mesh diagnostics, not inferred from A and B coverage.
- Do **not** move the band to follow the computed free boundary on finer meshes: that would make the size field solution-dependent and break the frozen family.
- The archived B1 background target was 0.0085 µm and tip target 0.0034 µm, but these are **historical references only** until the new master field is approved and globally scaled.

### 4.2 Proposed levels (not solve authorizations)

For a master field whose band target is 0.00400 µm:

| Level | Global factor \(r^k\) | Band target h (µm) | Role | Authorization |
|---|---:|---:|---|---|
| S0-L0 | 1.0000 | 0.00400 | Reference level | Design candidate only |
| S0-L1 | 0.7071 | 0.0028284 | First refinement | Design candidate only |
| S0-L2 | 0.5000 | 0.00200 | Second refinement | Conditional candidate |
| S0-L3 | 0.3536 | 0.0014142 | Further asymptotic check | **Not authorized** |

These numbers resemble old B1-X/B2-X band targets, but **the new levels are not the archived meshes**: the full size field, interface treatment, and possibly geometry policy must be rebuilt consistently. No old area result may be inserted as an S0-Lk data point without a demonstrable exact policy match.

### 4.3 Realized-mesh quality tests

At each *authorized future* level, compute distributions of realized triangle sizes \(h_T\), not just target values, in each of:

1. Band core (distance ≤0.04 µm to frozen curve).
2. Outward transition region (0.04–0.06 µm).
3. Legacy windows A and B and their overlapping/inter-window strip.
4. Tip and immediate high-gradient region.
5. Body-cut neighborhood \(x\approx0.54\,\mu\mathrm m\).
6. Oxide/interface neighborhoods and global remainder.

Collect element counts, min/median/p90/p95/max \(h_T\), aspect/quality measures, and size ratio \(h_{k+1}/h_k\) by matching **physical regions** (not necessarily matching triangles). A proposed feasibility gate is that the regional median ratio be within **±10% of** \(r\) in regions that are intended to scale, with no persistent unscaled bottleneck. This is a *design threshold* to revisit based on Gmsh behavior, not a mathematical certification threshold.

For strict nested mesh construction, record parent-child refinement mapping and constraint/boundary conformity; for independently remeshed levels, document mesh alignment diagnostics and acknowledge the loss of exact nesting. Controlled refinement can exist without topological nesting, but the absence of nesting increases the burden of demonstrating a stable asymptotic regime.

---

## 5. Resource budget before mesh-generation approval

**Do not extrapolate directly from B1-X→B2-X element counts**: those runs refined only selected regions, while the proposed \(h_k=r^k h_0\) scales the **whole domain**. For 2D meshes, a first-order estimate is \(N_T\propto h^{-2}\), so a global step at \(r=1/\sqrt{2}\) can approximately **double** elements; two steps can approximately **quadruple** them. Solver memory and wall time need not scale linearly with element count. Active constraints may grow differently than DOFs.

### 5.1 Required pilot estimation

At CP-S0.3 produce for **each candidate level**:

- Expected triangles, P2 total/free DOFs, silicon Bernstein rows, interface nodes, nonzeros in K and constraint matrices.
- Estimated OSQP memory, factorization/preconditioner memory, scaled-MINRES memory, peak RAM, disk for checkpoints and diagnostic arrays.
- Estimated solve and area-bound wall time as *ranges*, separately; include a factor for active-set repair and mixed-triangle area recursion.
- Hardware resources and safety margin, including workspace/disk checkpoint limits.
- A conservative abort threshold (e.g., projected peak RAM >70–80% of physically available RAM) and evidence explaining the chosen percentage.

**Known reference (not an S0 budget):** B2-X checkpoint arrays: primal `x(181699,)` and dual `y(416940,)`, OSQP 5,450 iterations; B2-X area run 697 s and 887 mixed elements. Those figures are insufficient to predict a globally scaled L2/L3 solve without mesh-generation pilots.

### 5.2 Resource gates

- **PASS:** feasible mesh/solver/area estimates with headroom and recoverable checkpoints; explicit authorizations can be considered.
- **HOLD:** uncertainty is too large; conduct read-only profiling or a lightweight mesh pilot (only after separate permission).
- **FAIL:** credible projection exceeds compute/RAM/disk limits; redesign before production.

---

## 6. Checkpoints and hard gates

**Rules common to every checkpoint:** write `PASS`, `HOLD`, or `FAIL` in an immutable checkpoint record, with timestamp, commit, command line, input hashes, output paths, summary metrics, reviewer, and next authorized action. Do not advance on `HOLD` or `FAIL`. Use distinct folders for candidate and accepted outputs; keep the closed B0/B1/B2 artifacts unchanged.

### CP-S0.0 — Archive lock and scope declaration (design phase)

**Evidence:** Verify the provenance list `15ec811`, `2478f82`, `4e5e6d4`; reference `P2_3_B2X_RESULTS.md` §9 and the archived F7 JSON; record current branch/commit and confirm no current solves are to be restarted.

**PASS iff:** accepted numerical intervals and conclusions are reproduced **verbatim in the plan**, and the old three-mesh chain is marked CLOSED.  
**HOLD if:** files/hashes are unavailable or a provenance inconsistency appears.  
**On PASS:** proceed with a *design-only* configuration review.

### CP-S0.1 — Scientific target and geometry policy

**Evidence:** State the target scalar \(A_h\) definition (neutral-area set, boundary inclusion, epsilon and estimator policy), fixed physical PDE, G-fixed vs G-physical geometry path, and separate front diagnostics. State what continuum claim could eventually be supported.

**PASS iff:** the target domain is unambiguous, the geometry-error policy is explicit, and area-vs-front conclusions are separated.  
**HOLD if:** polygonal geometry and physical curved geometry are being conflated.  
**On PASS:** draft frozen numerical formulation.

### CP-S0.2 — Numerical formulation and source-integration freeze

**Evidence:** Hash P2 basis/stiffness code, Dirichlet/Neumann/interface treatment, constraint builder, exact donor clipper/quadrature, OSQP/MINRES configuration, area and front routines. Include unit tests:

- Fully charged element: local load equals \(\rho A[0,0,0,1/3,1/3,1/3]\) to floating-point tolerance.
- Fully body-side element: exactly zero donor load.
- Cut element: degree-2 quadrature in **parent barycentric coordinates** and positive physical subareas.
- Global \(\sum_i f_i/\rho\) matches exact clipped charged area within a documented roundoff budget.
- All partly charged triangles contribute, including centroid-body-side straddlers; constraint counts never used as a donor mask.
- Corrected Bernstein area estimator passes parent/child frame-consistency fixtures; front extractor uses full P2 shapes and correct barycentric order.

**PASS iff:** all tests pass and references/config hashes are frozen.  
**HOLD if:** any mismatch or undocumented code/config change.  
**On PASS:** specify master mesh field.

### CP-S0.3 — Master field, family construction and resource model

**Evidence:** Versioned master curve and distance-field coordinates, h0 functional prescription over **all** regions, scale ratio, domain/geometry handling, nested-vs-remeshed policy, candidate level matrix, resource forecast, and mesh diagnostic script specification.

**PASS iff:** one deterministic mesh-family policy is defined, the intended local/global scope is explicit, geometry handling is consistent, and resource projections are credible with headroom.  
**HOLD if:** the tip/background/interface are accidentally unscaled or geometry changes silently.  
**On PASS:** assemble an S0 package for review. **Do not generate a mesh unless separately approved.**

### CP-S0.4 — Protocol review and design freeze

**Evidence:** A reviewer-signed version of this protocol; level table; config/geometry/size-field hashes; expected storage/runtime; risks; and explicit authorization text.

**PASS iff:** the S0 protocol is accepted in writing.  
**HOLD if:** decisions are unresolved or reproducibility is incomplete.  
**On PASS:** S0 design is CLOSED. Future pilot meshes require a **new** authorization.

### CP-M1 — Pilot mesh / geometry invariance (future, not authorized)

**Entry condition:** Separate approval to generate at least the pilot meshes.  
**Evidence:** same intended polygonal/curved geometry policy, identical physical cut-off, silicon and clipped donor areas, mesh quality, region-wise h ratios, verified frozen-band coverage, and any nested-parent maps.

**PASS iff:** mesh realizations match their intended size field and geometry policy; cut-source conservation and interface conformity checks pass.  
**HOLD if:** background/tip grading or the inter-window strip fails regional checks; geometry areas drift outside the predeclared budget.  
**FAIL if:** degenerate/inverted triangles, materially wrong tagging, or inconsistent charged support.  
**On PASS:** production solves remain **not authorized** until a distinct compute gate passes.

### CP-M2 — Solve authorization and reproducible checkpoints (future)

**Entry condition:** Separate solve approval after CP-M1 and updated measured resource estimates.  
**Evidence:** OSQP checkpoint, reload verification of primal/dual arrays and dimensions, unaltered assembly/BC/constraints, solver log, hardware peak-RAM/wall-clock records.

**PASS iff:** OSQP returns a converged status with usable residuals and verified checkpoint; no resource thresholds exceeded.  
**HOLD if:** solver stalls, checkpoint cannot reload, or resource estimate is invalidated.  
**On PASS:** perform scaled-MINRES polish; OSQP status alone is not F7.

### CP-M3 — Full-set F7 algebraic certification (future)

**Evidence:** rank-reduced active-set / MINRES details (including corrected residual sign and `info=0`); full-set stationarity, complementarity, min slack, violations, lambda min, inactive lambda, oxide residual, overshoot, active set and repairs; archived JSON.

**PASS iff:** all existing F7 gates pass under **unchanged documented tolerances**; small admissible numerical negative slacks are treated under those tolerances, not silently clipped.  
**HOLD if:** any gate fails or a tolerance was relaxed post hoc; investigate, repair and rerun before comparing area.  
**On PASS:** authorize area postprocessing only within the approved study budget.

### CP-M4 — Bernstein-area enclosure and estimator QA (future)

**Evidence:** interval \([L_k,U_k]\), width, area tolerance, mixed triangle count, subdivisions, max observed depth, minimum relevant parent triangle area, recursion stop reasons, runtime, and nesting against looser tolerance if available.

**PASS iff:** lower ≤ upper; depth is demonstrably **not** a truncation mechanism; certified bounds pass internal estimator invariants; widths are small enough for the specific **next** difference comparison. `AREA_TOL=1e-13` is an initial parameter, **not** a guarantee of a 1e-13-wide area interval.

**HOLD if:** area intervals too broad to distinguish changes, even if estimator terminated normally. In that case tighten **area postprocessing** first; no new solve is justified by interval width alone.  
**On PASS:** compute comparison intervals with outward rounding.

### CP-M5 — Controlled differences and convergence eligibility (future)

For \(A_k\in[L_k,U_k]\), calculate

\[
D_k=A_{k+1}-A_k\in[L_{k+1}-U_k,\;U_{k+1}-L_k].
\]

For **signed** differences, first check whether the intervals exclude zero and whether successive differences have compatible sign. If a nonzero ratio is meaningful, use conservative interval division to bound \(D_{k+1}/D_k\); for magnitude shrinkage, compare \(\sup|D_{k+1}|\) with \(\inf|D_k|\). If either difference contains zero, do **not** infer its log ratio or an empirical order.

Only if the study has a consistent refinement ratio, coherent geometry, robustly resolved differences, credible dominance of spatial rather than geometry/source/postprocessing error, and multiple levels in an apparent asymptotic regime, consider

\[
p_k=\frac{\ln|D_{k+1}/D_k|}{\ln r}.
\]

**PASS for a *sensitivity statement*** if a stated interval inequality is certified.  
**PASS for *Richardson eligibility*** only if the stronger family-consistency and asymptotic-evidence requirements are independently documented.  
**HOLD** if ratios depend on intervals containing zero, spatial changes are dominated by area uncertainty, or geometry errors are uncontrolled.  
**Never claim continuum convergence solely from decreasing two increments.**

### CP-M6 — Front and regional response (future)

**Evidence:** full quadratic neutral sets at frozen x-stations; epsilon sweep; satellite/shallow-gap qualification; regional neutral areas including the A/B/inter-window/remainder decomposition; spatial maps and sensitivity of \(\partial_yu\) near crossing points.

**PASS iff:** front data are presented explicitly as separately conditioned diagnostics, with uncertainty qualification.  
**HOLD if:** vertex-only interpolations, ambiguous multi-span crossings, or near-zero gradients are interpreted as certified front movement.  
**On PASS:** write scientific disposition without overstating front convergence.

### CP-M7 — Stop/go and archive (future)

**Evidence:** every checkpoint result; level-by-level provenance table; solver and area uncertainty; resource use; interval conclusions; any independent geometry study needed; limitations and project status.

**PASS iff:** a reviewer accepts the **exact** level of claim justified by the evidence.  
**Stop** if a full-set F7 gate cannot be reached within budget, the resource ceiling is breached, mesh-policy consistency fails, or the approved maximum level is reached. A stalled/asymptotically inconclusive study is a legitimate reportable outcome, not permission to refine indefinitely.

---

## 7. Conservative interval arithmetic and stopping criteria

### 7.1 Basic rule

Store all area endpoints at full output precision, with conservative outward rounding. For intervals \(A_i=[L_i,U_i]\) and \(A_j=[L_j,U_j]\), define:

\[
A_j-A_i=[L_j-U_i,\; U_j-L_i].
\]

The width of this difference is \((U_j-L_j)+(U_i-L_i)\). No midpoint subtraction is a substitute for that enclosure.

### 7.2 Quantitative *proposed* postprocessing gate

For consecutive area intervals, require

\[
\mathrm{wid}(D_k)\le\eta\,|\mathrm{mid}(D_k)|,
\]

with **provisional \(\eta=0.1\)** only when \(D_k\) excludes zero. For precise rate analysis target a significantly stricter \(\eta\) after checking available compute cost. This is a **planning heuristic**, not proof that spatial discretization error dominates. Preserve exact interval endpoints and also report \(\mathrm{wid}(D_k)\) separately.

If widths overwhelm differences, first tighten the Bernstein estimator or audit its recursion policy. If the area estimator cannot resolve meaningful changes within the allowed budget, stop and report that limitation.

### 7.3 Required stop conditions

- **Immediate STOP/HOLD:** inconsistent physics, geometry policy, donor support, mesh size field, area estimator, constraint builder, or boundary tagging.
- **Immediate STOP/HOLD:** F7 failure or unresolved active-set/KKT issues.
- **STOP:** peak memory/disk/runtime exceeds preapproved ceilings or the next step is not authorized.
- **STOP AND REPORT:** successive spatial differences do not shrink or their signs/ratios oscillate without a convincing explanation; do not infer an order.
- **STOP AND REPORT:** intervals cannot resolve the desired comparison at feasible postprocessing cost.
- **STOP:** approved family reaches its maximum level, regardless of apparent order; extension requires a new review.
- **Scientific restraint:** F7 PASS is necessary for trusted discrete solutions, but neither sufficient for spatial convergence nor for a converged free-boundary location.

---

## 8. Suggested directory layout and reproducibility schema

The following is a **proposal**, not a claim that files already exist:

```text
p2_3_spatial_S0/
  protocol/
    P2_3_SPATIAL_CONVERGENCE_PROTOCOL_S0_CHECKPOINTS.md
    DECISIONS.md
    model_and_geometry_manifest.json
    size_field_manifest.json
    solver_and_area_manifest.json
    resource_forecast.csv
  checkpoints/
    CP-S0.0.json
    CP-S0.1.json
    CP-S0.2.json
    CP-S0.3.json
    CP-S0.4.json
    # CP-M1 ... CP-M7 only after corresponding approvals
  levels/             # future, not generated by S0
    L0/
    L1/
    L2/
  results/            # future
    areas_intervals.csv
    interval_differences.csv
    region_front_diagnostics.csv
    final_disposition.md
```

Minimal checkpoint record fields:

```json
{
  "checkpoint": "CP-S0.2",
  "status": "HOLD",
  "phase": "design_only",
  "utc_timestamp": null,
  "git_commit": null,
  "source_hashes": {},
  "config_hashes": {},
  "commands": [],
  "tests_and_metrics": {},
  "output_artifacts": [],
  "deviations": [],
  "reviewer": null,
  "next_authorized_action": "none_until_review"
}
```

The sample status `HOLD` is intentionally conservative; the checklist has not been executed merely because this document exists.

---

## 9. Design-review checklist (S0 only)

- [ ] **CP-S0.0** Closed B0-X/B1-X/B2-X provenance verified and immutable.
- [ ] **CP-S0.1** Geometry target and geometry-error strategy selected.
- [ ] **CP-S0.1** Neutral area defined separately from epsilon-sensitive front geometry.
- [ ] **CP-S0.2** P2/stiffness/source/constraints/solver/estimator code and tolerances hashed.
- [ ] **CP-S0.2** Exact partial-element donor quadrature tests and charge-support audit specified.
- [ ] **CP-S0.3** Frozen L1b band curve coordinates/provenance recoverable.
- [ ] **CP-S0.3** Whole-domain h0 function and r declared; global vs localized study labeled accurately.
- [ ] **CP-S0.3** Physical tip, body-cut, interfaces, strip and oxide-region grading covered.
- [ ] **CP-S0.3** Candidate level mesh/DOF/RAM/disk/runtime estimates documented with uncertainty.
- [ ] **CP-S0.4** Stop/go criteria and budgets approved in writing.
- [ ] **CP-S0.4** Future mesh and solve permissions explicitly separated.
- [ ] **CP-S0.4** Final design-frozen commit/hash recorded.

**S0 exit condition:** All required design checkpoints PASS and an approved signed/committed protocol exists. **S0 exit does not mean any production mesh or solve has started.**

---

## 10. Exact decision requested at next review

**Requested approval now:** **GO for S0 design/review only**, including offline manifesting, defining h0, choosing G-fixed vs G-physical, and preparing a resource estimate.  
**Not requested and not granted:** producing a new mesh, running OSQP/MINRES, changing the existing solver, rerunning the closed B2-X study, starting the old third band level, or resuming P1/Stage 3.

**Proposed progression:**

```text
CP-S0.0 archive lock
  -> CP-S0.1 scientific target & geometry fork
  -> CP-S0.2 numerical source / obstacle / estimator freeze
  -> CP-S0.3 master field + levels + resource feasibility
  -> CP-S0.4 protocol acceptance
  -> [SEPARATE AUTHORIZATION]
  -> CP-M1 pilot meshes and geometric audit
  -> [SEPARATE SOLVE AUTHORIZATION]
  -> CP-M2 OSQP checkpoints
  -> CP-M3 full F7
  -> CP-M4 certified area enclosures
  -> CP-M5 interval comparison and Richardson eligibility check
  -> CP-M6 front/regional diagnostics
  -> CP-M7 accepted scientific disposition and archive
```

**Closing scientific boundary:** A controlled mesh family and shrinking, resolved neutral-area increments can support an empirical convergence study. An asymptotic order requires stronger evidence; a physical continuum limit also requires adequate control of geometric approximation and any remaining model or estimator biases. No new result should be characterized as a Richardson estimate until those conditions have been separately established.

---

### Archive reminder

`15ec811` = exact-source four-area `1e-13` baseline audit; `2478f82` = B2-X enclosure and strict increment comparison; `4e5e6d4` = closing acceptance and B2-X F7 archival. Preserve all certified centroid and exact-source artifacts unchanged. The current activity described in this document is **protocol preparation only**, with **no solve running or queued by this document**.
