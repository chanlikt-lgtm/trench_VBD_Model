# P2.3-G1 — Arc Resolution and Geometry Sensitivity Controls

**Status:** DESIGN ONLY. No mesh generation, no solve, no area evaluation, no
chord-policy implementation is performed or authorised here.
**Reference study:** P2.3-S0, CLOSED and FROZEN at `ec4bc80`.
**Reference level:** S0-L1 (`e1c8a2d` solve, `803dcff` area) —
`A_S0-L1 ∈ [0.022217431654472784, 0.022217790436636530] µm²`, F7 PASS.
**Disposition:** S0 CLOSED · P2.3 continuum convergence OPEN · P1 FROZEN ·
Stage 3 SUSPENDED · S0-L3 NOT AUTHORISED.

## 1. Question

S0 established shrinking increments on one frozen polygon. Two error terms were
left open: the **nonlocal** effect of the arc-adjacent refinement shortfall, and
the **fixed-polygon vs curved-domain** geometry error. G1 separates them.

Two distinct perturbations are required, because a single finer-chord mesh changes
both at once and therefore measures only their sum:

| Control | Polygon | Arc-adjacent resolution | Measures |
|---|---|---|---|
| **G1-C0** | frozen S0 polygon G₀ | archived (L1) | reference |
| **G1-C1** | **G₀ unchanged** | **refined** | δA_mesh = A(G₀,M₁) − A(G₀,M₀) |
| **G1-C2** | **G₁ (finer chords)** | matched to C1 | δA_geometry = A(G₁,M₁′) − A(G₀,M₁) |

## 2. G1-C1 — same polygon, finer arc resolution

**Geometry.** The S0 polygon is untouched: the four arcs remain 64 straight chords
with 65 vertices each. Refinement is obtained by **collinear subdivision of those
existing straight edges**, which cannot move the boundary — verified in S0 by the
1.110e-16 boundary drift between L0/L1/L2.

**Size field.** Everything at S0-L1 values except a local band:

```
arc band      |distance to the two arc polylines| <= 0.02 µm  → h = 0.0020
grading       0.02 – 0.04 µm                                 → linear to L1 value
elsewhere     S0-L1 size field, unchanged
```

`h = 0.0020` is the L2 band target, so C1 asks for L2-grade resolution **only**
in the arc band while holding the rest of the mesh at L1.

**What it isolates.** If the arc shortfall acts only through direct local area,
δA_mesh should be at or below the measured arc-region ceiling (≤1.5e-7 µm²).
If it acts nonlocally through the elliptic coupling and active set, δA_mesh will
be larger — and that is the quantity S0 could not bound.

**Expected cost.** The band contains ~500–900 elements at L1 sizes; refining to
0.0020 adds roughly 1–3 k elements, so the mesh stays near L1: ~462 k nodes,
~230 k elements, ~462 k compacted DOFs, ~712 k Bernstein rows.

## 3. G1-C2 — finer-chord polygon

**Geometry.** Increase the arc discretisation from 64 to **128 chords per 90°**,
choosing the new vertices so that the 64-chord vertices are a **subset** — i.e. add
one midpoint per existing chord. The two polygons are then nested approximations
of the same circle and the change is a monotone geometry refinement rather than
an unrelated re-mesh.

```
64 chords   chord 2.454e-03 µm   sagitta 7.53e-06 µm
128 chords  chord 1.227e-03 µm   sagitta 1.88e-06 µm  (4× smaller)
```

**Size field.** Identical to C1, so the mesh-resolution perturbation is nominally
the same and the geometry change is isolated as far as practical.

**Why this is not a pure geometry error.** Changing the chord count also changes
boundary vertex count (65 → 129 per arc), local triangulation and possibly the
active set. δA_geometry therefore carries a residual mesh-layout component that
must be reported, not assumed away. A scale prediction, not a bound: the Si-area
deficit of the 64-chord polygon is −1.631e-05 µm²; recovering three-quarters of it
at 128 chords suggests a domain-area change of order +1.2e-05 µm², which is ~8 % of
ΔA_01. If the neutral area responds proportionally, geometry uncertainty is of the
same order as the quantity being used to argue convergence — which is precisely
why this control matters.

## 4. Pre-registered quantities

```
δA_mesh      = A(G₀, M₁)  − A(G₀, M₀)
δA_geometry  = A(G₁, M₁′) − A(G₀, M₁)
```

Both by endpoint subtraction, both at `AREA_TOL = 1e-13`, both only after a full-set
F7 PASS on the new solutions. `M₁′` must be matched to `M₁` (same size field, same
background, same grading); the residual triangulation near the arc is reported.

Secondary diagnostics, per control: arc-region h distributions, minimum angle,
aspect tail, zero-area count, per-face winding, P2 midside deviation, exact donor
conservation, and the arc-region area enclosure at the same fixed region used in
CP-L2-A6.

## 5. Checkpoint gates

| Gate | Required evidence | PASS condition |
|---|---|---|
| **G1-0** design freeze | reference level, both geometry definitions, size fields, diagnostic list, resource plan | one deterministic specification, no ambiguity |
| **G1-1** same-geometry control | polygon identity: Si area to ≤1e-12, boundary drift ≤1e-12; arc h median at least 25 % finer than C0; quality | no boundary movement; band genuinely refined |
| **G1-2** finer-geometry control | new polygon hash, chord count, sagitta, Si/charged area drift vs G₀; quality; triangulation-difference report | geometry change quantified and nested |
| **G1-3** resources | nodes, elements, DOFs, rows, assembly RSS, solve forecast vs 70 % ceiling | headroom with margin |
| **G1-4** solve authorisation | separate decision after G1-0…G1-3 | explicit GO |
| **G1-5** numerical comparison | F7-certified solutions; δA_mesh and δA_geometry as intervals | intervals resolve or are reported unresolved |

Enclosure widths of order 3e-07 µm² are already achievable (S0), so differences of
1e-6 and above are resolvable; a δ near 1e-7 may not be.

## 6. Resource estimate (planning only)

| Item | G1-C1 | G1-C2 |
|---|---|---|
| Elements (est.) | ~230 k | ~232 k |
| Compacted DOFs (est.) | ~462 k | ~465 k |
| Bernstein rows (est.) | ~712 k | ~716 k |
| Assembly peak RSS (est.) | ~350 MB | ~360 MB |
| Solve peak RSS (est.) | ~1.0–1.2 GB | ~1.0–1.2 GB |
| Wall clock (est.) | ~20 min | ~20 min |

Anchors: L1 measured 998 MB peak and 1,191.8 s; the 70 %-of-physical ceiling is
17,127 MB. Forecasts are not upper bounds.

## 7. Decisions required before any mesh is generated

1. **C1 band width and target.** Proposed 0.02 µm / h = 0.0020 — this is a
   deliberate choice, and a smaller band or a coarser target weakens the test.
2. **C2 chord count.** Proposed 128 (nested with the 64-chord set). A higher count
   reduces the geometry perturbation further but grows the boundary-node
   difference against C1.
3. **Matching criterion for M₁′.** Proposed: same size field and background, with
   the boundary-vertex difference (65 vs 129 per arc) explicitly reported.
4. **Acceptance thresholds** in G1-1/G1-2 as written, or amended.

## 8. Authorisation boundary

**GO:** protocol design, control specification, geometry policy planning, resource
estimation — this document.

**HOLD:** mesh generation, OSQP/MINRES, F7, area evaluation, chord-policy
implementation, and any modification of frozen S0 artefacts.

S0 remains closed at `ec4bc80`. G1 must use separate directories, separate hashes
and separate provenance; nothing in G1 may overwrite or reinterpret an S0 record.
The purpose is to decide whether arc-related effects are small enough to justify a
future continuum study — **not** to retrofit a convergence order onto S0.
