# P2.3 Latest Session — B0 OSQP Checkpoint and Body-Cut Charge Audit

**Session scope:** Latest B0 results through the body-cut audit and the decision to finish B0 scaled polishing before further numerical experiments.  
**Project:** Conforming Gmsh P2 Bernstein-obstacle electrostatic model of a shield-gate trench MOSFET.  
**Numerical disposition:** B0 OSQP solved and checkpoint verified; scaled MINRES polishing running at the time of this report; **B0 not yet F7-certified**. B1 is F7-certified. B2 is on hold.  
**Physical/model disposition:** P2.3 **OPEN**; P1 frozen; Stage 3 suspended.  
**Source of record:** User-reported run outputs and audit measurements in the latest conversation. This record does not claim independent execution or verification of production code.

## 1. Executive summary

The B0 uniform-background control was constructed to compare with B1 under the same intended CAD geometry, charge-classification rule and mesh-generation policy. The two meshes have the same reported total silicon and exact geometrically clipped charged area. Unlike B1, B0 does not apply the frozen neutral-front refinement band; B1 has a much finer band. The design change nevertheless modifies triangulation near the body cut at \(x=0.54\,\mu\mathrm m\), even though the realized local element sizes there are the same.

The body-cut audit found that the centroid whole-element donor rule overcharges some cut elements and undercharges others. These two errors nearly cancel, making the *net* donor-load error sensitive to connectivity and triangle arrangement. From B0 to B1, the integrated donor load rises by \(0.303498\) in the assembly's reported units (about \(1.54\times10^{-4}\) relative), despite identical exact charged geometry. Consequently, B0→B1 is an interpretable test of **enabling the refinement band within the existing centroid-classification formulation**, but it is *not* a pure boundary-resolution comparison holding the assembled donor source fixed.

B0 OSQP has solved successfully and its primal/dual checkpoint was saved and reload-verified. The complete F7 gate, bounded neutral area and B0→B1 area difference are **not yet available** in this session.

## 2. B0 and B1 construction comparison

### 2.1 Geometry and mesh

| Metric | B0 | B1 | Interpretation |
|---|---:|---:|---|
| Nodes | 212,779 | 231,568 | Different refinement level |
| Triangles | 105,774 | 114,721 | B1 has 8,947 more |
| Silicon area (µm²) | 1.954559924360 | 1.954559924360 | Identical to printed precision |
| Exact clipped charged area (µm²) | 1.575269192653 | 1.575269192653 | Identical to printed precision |
| Centroid-classified charged area (µm²) | 1.575255252958 | 1.575498331328 | Different discrete source support |
| Near-band elements | 2,282 | 10,160 | B1 has 7,878 additional elements in this band |
| Near-band median \(h\) (µm) | 0.00850 | 0.00400 | Controlled local-size contrast |
| Region A / B / inter-window strip median \(h\) (µm) | 0.00850 / 0.00835 / 0.00850 | Approximately 0.00400 throughout | B1 refines the full targeted boundary band |
| Minimum angle | 35.95° | 31.09° | Both reported geometrically healthy |

The 7,878 additional near-band elements represent approximately **88%** of the total 8,947 additional triangles. The nominal geometric and physical domains are unchanged. Matching areas provide strong evidence of identical *area measures* of the polygonal silicon and clipped donor regions, although area equality by itself is not a proof of vertexwise polygon identity.

### 2.2 Why the B0 control is useful

Previous comparisons between L1b and B1 simultaneously changed near-front resolution and the realized far-field/background grading. By contrast, B0 and B1 share the uniform-background policy and differ intentionally by activating the fixed neutral-boundary size field. The resulting comparison isolates the effect of the **mesh-generation intervention**, including global remeshing consequences. It does not hold each element's charge mask or assembled load invariant.

## 3. B0 OSQP solve and durable checkpoint

The user-reported OSQP output was:

```text
[B0-A] solved, 5,025 iterations, obj -1.052110e+05
       pri 4.051e-08, dua 3.763e-09
       checkpoint reload verified: x(142473,) y(303672,)
```

| Field | Result |
|---|---:|
| OSQP status | `solved` |
| Iterations | 5,025 |
| Objective | \(-1.052110\times10^5\) |
| OSQP primal residual | \(4.051\times10^{-8}\) |
| OSQP dual residual | \(3.763\times10^{-9}\) |
| Primal checkpoint length | 142,473 |
| Dual checkpoint length | 303,672 |
| Checkpoint reload | Verified |
| Scaled MINRES polishing | Running at the time of the report (reported PID 34084) |
| Full-set F7 | **Pending; no certification implied by OSQP status** |

The B0 mesh/assembly had previously been reported as 212,441 DOFs, 142,473 free DOFs and 303,672 Bernstein constraints. **The primal checkpoint must remain preserved independently of polishing.** Neither an F7 residual nor a B0 neutral-area estimate is inferred from these OSQP figures.

## 4. Body-cut audit at \(x=0.54\,\mu\mathrm m\)

### 4.1 Reported numerical measurements

| Diagnostic | B0 | B1 |
|---|---:|---:|
| Cut-straddling Si elements | 180 | 174 |
| Total area of straddling elements (µm²) | \(5.598270\times10^{-3}\) | \(5.389783\times10^{-3}\) |
| Exact charged Si area (µm²) | 1.575269192653 | 1.575269192653 |
| Centroid-classified charged area (µm²) | 1.575255252958 | 1.575498331328 |
| Net centroid-mask area error (µm²) | \(-1.394\times10^{-5}\) | \(+2.291\times10^{-4}\) |
| Overcharged area (µm²) | \(+5.203172\times10^{-4}\) | \(+6.856357\times10^{-4}\) |
| Undercharged area, magnitude (µm²) | \(5.342569\times10^{-4}\) | \(4.564970\times10^{-4}\) |
| Elements within \(\lvert x-0.54\rvert\le0.05\,\mu\mathrm m\) | 4,425 | 4,421 |
| Local median element size (µm) | 0.00850 | 0.00850 |
| Sum of discrete donor load \(\sum_i f_i\) | 1966.806755 | 1967.110253 |

The sign convention is **overcharge minus undercharge**, so the approximate area error is

\[
E_{\rm mask}=A_{\rm over}-A_{\rm under}.
\]

At B0, the overcharged and undercharged contributions are nearly equal; the negative net error is only about \(1.39\times10^{-5}\,\mu\mathrm m^2\). At B1, a redistribution of straddling triangles yields a positive net error of about \(2.29\times10^{-4}\,\mu\mathrm m^2\).

### 4.2 Integrated source difference

The change in centroid-classified charged area is

\[
\Delta A_{\rm charge,centroid}
=1.575498331328-1.575255252958
=0.000243078370\,\mu\mathrm m^2.
\]

The reported load sums give

\[
\Delta\sum_i f_i
=1967.110253-1966.806755
=+0.303498,
\]

or a relative change of approximately

\[
\frac{0.303498}{1966.806755}=1.54\times10^{-4}
=0.0154\%.
\]

The same donor coefficient and P2 load assembly are reported for both meshes. The audit attributes the load difference to mesh-dependent whole-element centroid classification, *not* to a change in continuum donor support or the intended rule.

### 4.3 Mechanism and its evidence boundary

The cut neighborhood contains roughly the same number of elements (4,425 versus 4,421) with the same median size (0.00850 µm), but their connectivity/placement differs. The refinement field introduced near \(x\approx2.4\,\mu\mathrm m\) leads to a different global triangulation, including the set of triangles that cross the body cut. The measured overcharge and undercharge magnitudes are around \(5\times10^{-4}\,\mu\mathrm m^2\) each, while their net may be an order of magnitude smaller. Modest variations in those contributions can therefore alter the sign and size of the net error.

**Supported conclusion:** The error of whole-element centroid donor classification is sensitive to triangle arrangement, even at approximately fixed local \(h\). The data support global remeshing as a plausible mechanism. Demonstrating the exact propagation of mesher topology from downstream band controls to the upstream body cut would require vertex/connectivity correspondence, but is not needed for the present source-error audit.

**Important non-equivalence:** The \(0.0154\%\) donor-load change is *not* a \(0.0154\%\) neutral-area error. Obstacle active-set changes and field conditioning can alter the response; no linear area/load relation has been established.

## 5. B1 reference already F7-certified

B1 was previously certified at repository commit **`0f865ee`**. The full-set F7 values were:

| F7 diagnostic | B1 result |
|---|---:|
| Independent active rows | 16,719 |
| Final stationarity | \(8.929\times10^{-11}\) |
| Complementarity | \(9.417\times10^{-18}\) |
| Minimum slack | \(-2.420\times10^{-15}\) |
| Full-set violations | 0 |
| Minimum multiplier | 0 |
| Oxide residual | \(4.405\times10^{-13}\) |
| Potential overshoot | \(1.314\times10^{-15}\) |
| Gate | **PASS** |

The successful B1 gate included one active-set repair: two multiplier-negative rows (approximately \(-1.04\times10^{-6}\)) were dropped and the re-solve passed.

B1's certified whole-domain Bernstein area enclosure is

\[
A_{\mathrm{B1}}\in
[0.021987652318,\;0.021991846424]\,\mu\mathrm m^2.
\]

No B0 area interval was available at the time of this session, so **no B0→B1 neutral-area difference can yet be calculated**.

## 6. Interpretation of the pending B0→B1 comparison

Once B0 passes full-set F7, define

\[
\Delta A_{01}=A_{\mathrm{B1}}-A_{\mathrm{B0}}.
\]

If \(A_{\mathrm{B0}}\in[L_0,U_0]\) and \(A_{\mathrm{B1}}\in[L_1,U_1]\), interval propagation requires

\[
\Delta A_{01}\in[L_1-U_0,\;U_1-L_0].
\]

The difference will measure the combined response to activating the neutral-boundary refinement field under the existing centroid-charge implementation, including (i) improved boundary resolution, (ii) changed triangulation elsewhere, (iii) the measured change in discrete donor load, and (iv) any active-set changes. It **cannot** yet be identified as a pure local-resolution effect with fixed assembled source.

The next controlled scientific question is whether the resulting neutral-area change remains after replacing the centroid source approximation with **exact partial-element donor-load integration**. Such a control must integrate every P2 basis function over the geometrically clipped charged portion of each cut triangle; substituting the exact scalar charged area alone is insufficient.

Do **not** modify the source formulation within the currently running B0/B1 comparison. Any exact-integration experiment is separate and must be clearly labeled as a modified discrete source model.

## 7. Pending validation gates and recommended sequence

1. **Finish B0 scaled MINRES** and residual correction, preserving the original reload-verified OSQP checkpoint.
2. **Apply complete F7 to B0**: original-coordinate stationarity, all Bernstein slacks, all inequality multipliers, complementarity, oxide equilibrium, overshoot and active-set consistency. Do not accept Krylov convergence alone.
3. **Only after F7 PASS**, calculate B0's bounded whole-domain neutral area using the corrected canonical quadratic extractor and Bernstein subdivision.
4. **Propagate B0/B1 intervals** rather than using area midpoints as exact values. Include regional A/B/remainder decomposition and the established area-conservation enclosure checks.
5. **Inspect front cuts** at \(x=2.35,2.40,2.41\,\mu\mathrm m\), including threshold-dependent satellites and gap barriers; do not interpret adjacent shallow fragments as independently certified physical pockets.
6. **Review B0→B1 attribution** with the measured body-cut source discrepancy explicitly included. Then decide whether to authorize B2 or prioritize an exact partial-element charge-load control.

### Decision gate after B0

- **Large, non-overlapping positive area increment:** strong response to enabling the refinement-band *mesh intervention*; prioritize exact donor integration before assigning the magnitude solely to boundary resolution.
- **Area intervals overlap zero increment:** inspect redistribution by region and the donor-source contribution before proceeding to B2.
- **Negative increment:** the previous monotone-area pattern is not reproduced under this controlled band policy; investigate mesh-layout and source effects before another refinement level.

**B2 remains unauthorized** pending B0 F7 certification and the resulting bounded B0/B1 comparison.

## 8. Project status at close of this session

| Workstream | Status |
|---|---|
| B0 construction / geometry comparison | Reported PASS |
| B0 OSQP | Solved, checkpoint saved and reload-verified |
| B0 scaled MINRES | Running at the last user update; completion not established here |
| B0 full-set F7 | **PENDING** |
| B0 neutral-area enclosure | **PENDING** |
| B0/B1 body-cut audit | Complete; discrete donor-load difference quantified |
| B1 full-set F7 | Certified PASS |
| B1 neutral area | Certified interval reported above |
| B2 | On hold / not authorized |
| P2.3 | **OPEN** — spatial convergence not established |
| P1 | Frozen |
| Stage 3 tip-radius study | Suspended |

### Final disposition

**Proceed without changing B0's configuration.** Finish the running polish, certify against full-set F7, and obtain an interval-bounded B0→B1 neutral-area difference. Preserve the centroid-mask triangulation dependence as a declared source discretization confound. Do not start B2 and do not claim continuum area/front convergence from the available runs.

---

*Record assembled from the user's reported B0/B1 mesh and solver outputs and the existing certified B1 reference. No new simulation was executed to create this document.*
