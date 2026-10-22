# P2.3 Latest Session — Exact Donor Integration Assembly Validation (B0-X / B1-X)

**Session scope:** Exact partial-element donor-load implementation and validation; B0-X solver underway, B1-X queued.  
**Reference repository state:** B0/B1 certified baseline and associated records at commit `ba879b7` (as previously reported); no new commit for the X variant has been supplied in this session.  
**Project disposition:** P2.3 **OPEN**; P1 **FROZEN**; Stage 3 **SUSPENDED**; B2 **ON HOLD**.

> **Evidence boundary:** This document records results reported by the user and their technical interpretation. It does not claim independent execution or verification of the production scripts or their numerical archives. Results for B0-X and B1-X F7 or neutral areas are **not yet available**.

## 1. Purpose and experimental design

The centroid-source B0/B1 comparison demonstrated a large, strictly positive change in the Bernstein-bounded P2 neutral area when the frozen neutral-boundary refinement band was enabled. However, globally regenerated triangulations changed whole-element charge classification at the body cut `x = 0.54 µm`, despite identical continuum charged geometry.

The X variants replace centroid whole-element donor loading with exact partial-element P2 load integration, **without changing either mesh**:

| Variant | Mesh | Donor integration | Status |
|---|---|---|---|
| B0 | B0 original | Centroid whole-element | F7 certified; area bounded |
| B1 | B1 original | Centroid whole-element | F7 certified; area bounded |
| B0-X | **Same B0 mesh** | Clipped partial-element quadratic integration | Assembly tested; OSQP running at last report |
| B1-X | **Same B1 mesh** | Clipped partial-element quadratic integration | Assembly tested; queued at last report |

All other discretization and physical choices are to remain unchanged: geometry, material tags, permittivities, Dirichlet/Neumann conditions, stiffness matrix, obstacle/Bernstein constraints, body treatment, and solver acceptance criteria.

## 2. Implemented exact donor-load assembly

For each silicon triangle `T`, clip to the charged domain `T ∩ {x ≥ 0.54 µm}`. Fan-triangulate the clipped polygon, evaluate the P2 basis in the **parent barycentric coordinate frame**, and integrate using a three-point degree-two-exact rule. Subtriangle areas are computed from the barycentric determinant multiplied by parent area.

The desired six-component element load is

\[
  f_i^{X}=\rho\int_{T\cap\{x\ge0.54\}}N_i(x,y)\,\mathrm dA.
\]

For a fully charged P2 element of area `A` with standard vertex/midside ordering, the constant-source identity is

\[
  \mathbf f_T = \rho[0,\,0,\,0,\,A/3,\,A/3,\,A/3]^\mathsf T.
\]

Thus interior fully charged elements were already integrated exactly under the old method. The X change affects triangles **straddling the charge boundary**, with the caveat that fully uncharged elements should remain zero.

## 3. Exact-integration self-test results

Reported numerical tests:

| Test | B0-X assembly | B1-X assembly |
|---|---:|---:|
| Maximum full-element P2 load-identity error | `1.27e-21` | `1.69e-21` |
| Integrated load divided by rho, `Σf/ρ` (µm²) | `1.575269192653` | `1.575269192653` |
| Independently clipped exact charged area (µm²) | `1.575269192653` | `1.575269192653` |
| Reported source-integral conservation | **PASS** | **PASS** |

The full-element identity agrees at floating-point roundoff. The global load sum matches the exact clipped charged area to all twelve reported decimal digits on **both** meshes.

**What is established:** The user-reported checks provide strong algebraic evidence for the correct barycentric quadrature mapping, subtriangle-area handling, and global charge conservation. Unlike the original centroid mask, the X variants have equal *total integrated charge* for B0 and B1.

**What is not established by this check alone:** Equality of total charge does **not** imply equality of the spatially distributed assembled load vectors on the two finite-element spaces. The meshes have different nodes, basis supports, and connectivity. B0-X/B1-X still differ in the discretized source distribution and stiffness/active-set response. The X experiment removes the centroid whole-element misclassification error; it does not erase all mesh-dependent source-discretization effects.

### Additional assembly verification before final interpretation

The original assembly counts remain `charged=50612` and `body=12206` for B0-X at the reported stage. Confirm those counts are diagnostic/constraint-selection counts, **not** a filter that drops exact clipped donor contributions for triangles whose centroids fall below the cut. Every straddling triangle with nonzero clipped charged area must contribute appropriately. The global charge identity strongly suggests correct behavior, but the actual assembly selection should be checked explicitly.

## 4. Certified centroid-source reference

B0/B1 centroid-source configurations are previously reported as fully F7-certified, including B0's stationarity `5.824e-11` and B1's stationarity `8.929e-11`. Their Bernstein area intervals are:

| Mesh/formulation | Neutral area enclosure (µm²) |
|---|---|
| B0 centroid | `[0.021644935001, 0.021649266844]` |
| B1 centroid | `[0.021987652318, 0.021991846424]` |

The rigorously propagated difference for B0 → B1 under the centroid-source formulation is

\[
  \Delta A_C=A_{B1}-A_{B0}
  \in[3.383855,\,3.469114]\times10^{-4}\ \mu\mathrm m^2,
\]

corresponding to a positive relative change of approximately **1.56–1.60%** as previously reported.

The previously audited source perturbation was:

- B0 centroid charged area: `1.575255252958 µm²`.
- B1 centroid charged area: `1.575498331328 µm²`.
- Exact clipped area on both meshes: `1.575269192653 µm²`.
- Change in centroid-classified charged area: approximately `+2.43078e-4 µm²`.
- Change in integrated source load: `ΔΣf = +0.303498` in the assembled source units, approximately `1.54e-4` relative.

**Interpretation guardrail:** A ~100-fold ratio between *relative* neutral-area and integrated-source changes is **not** a mathematical bound on source-induced neutral-area sensitivity. A sharp-depletion obstacle solution and its neutral sublevel area can respond nonlinearly to small spatial source changes. The exact-source control is needed to measure that response.

## 5. Current solver sequence for exact-source variants

At last report:

```text
[B0X-A] dofs=212441 free=142473 constraints=303672
         charged=50612 body=12206
B0-X OSQP: iterating (pid 34536 at report time)
B1-X OSQP: queued behind B0-X
```

No B0-X or B1-X F7 validation, final neutral-area bound, or front-cut result has yet been reported.

Required steps **for each** X variant:

1. Complete OSQP and **save the primal/dual checkpoint immediately** before polishing.
2. Verify checkpoint reload and preserve mapping, mesh provenance, model configuration, and solver diagnostics.
3. Polish with scaled matrix-free MINRES and residual correction as in the certified B0/B1 workflow.
4. Apply the **full-set F7 gate**, not merely the active-set Krylov residual: stationarity, feasibility/slack, violations, nonnegative multipliers, complementarity, oxide residual, and overshoot.
5. Only after F7 PASS, evaluate the same Bernstein-bounded neutral area and canonical quadratic front-cut diagnostics.
6. Keep original centroid-source solutions and their audit records immutable. Store X artifacts with explicit `B0-X` and `B1-X` labels.

## 6. Decisive bounded comparisons

After X certification, calculate the intervals conservatively:

\[
  \Delta A_X=A_{B1-X}-A_{B0-X},
\]

\[
  \Delta_{\mathrm{correction}}=\Delta A_X-\Delta A_C.
\]

Also report within-mesh source corrections separately:

\[
  \delta A_{B0}=A_{B0-X}-A_{B0},
  \qquad
  \delta A_{B1}=A_{B1-X}-A_{B1}.
\]

For intervals `[L,U]`, subtract conservatively via `[L_a-U_b,\ U_a-L_b]`; do not treat midpoint differences as exact or independently certified.

| Observed exact-source outcome | Supported interpretation |
|---|---|
| `ΔA_X` close to centroid `ΔA_C`, strictly positive | Centroid source quadrature has a relatively small impact on the *observed increment* for these meshes. |
| `ΔA_X` materially different but still strictly positive | Source quadrature contributes to magnitude; a positive band-enabled response survives. |
| `ΔA_X` encloses zero or changes sign | Centroid-source approximation was a major confound for the original band comparison. |

Regardless of outcome, do **not** claim a continuum-converged neutral area or a causal, pure boundary-resolution effect solely from this pair; spatial mesh arrangement and active-set effects remain.

## 7. Decision and project status

**GO:** Finish B0-X, then B1-X, under the unchanged mesh and physical configuration. Complete full F7 certification and compare Bernstein-bounded areas using the formulas above.

**HOLD:** Do not start B2 or silently replace the frozen centroid-source baseline with exact-source assembly.

| Workstream | Current disposition |
|---|---|
| B0 centroid | F7 certified; bounded area established |
| B1 centroid | F7 certified; bounded area established |
| Exact P2 donor-load assembly | Self-tests and global charge conservation PASS, as reported |
| B0-X | OSQP running at last report; F7 pending |
| B1-X | Queued; F7 pending |
| X-variant neutral-area comparison | Pending |
| B2 | On hold |
| P2.3 | **OPEN** |
| P1 | **FROZEN** |
| Stage 3 | **SUSPENDED** |

**Main unresolved quantitative question:** How much of the certified B0 → B1 neutral-area increase survives removal of centroid-based donor-source integration error?
