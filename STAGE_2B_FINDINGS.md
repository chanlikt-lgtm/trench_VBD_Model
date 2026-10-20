# Stage 2B.2 / 2B.3 — Depletion-front convergence on a conforming trench mesh

**Status:** 2A passed · 2B.1 passed · 2B.2 open · 2B.3 open · Stage 3 suspended
**Date:** 2026-10-19
**Solver version:** P1 (frozen reference, `conforming_gmsh/`)

---

## 1. Scope and fixed inputs

Geometry, materials, boundary values and the depletion formulation were held
fixed across every mesh in this study:

| Quantity | Value |
|---|---|
| Tip radius r_tip | 0.10 µm (circular arc, unchanged) |
| Mesa half-width cap | A_MAX = P/2 − t_ox = 0.575 µm |
| Oxide thickness | 0.10 µm, **vertical offset** (not normal offset — see §6) |
| Doping N_D | 6.9e16 cm⁻³ |
| ρ = qN_D/ε_0 | 1248.5638 V/µm² |
| ε_Si / ε_SiO2 | 11.7 / 3.9 |
| V_D | 34.274208536 V |
| V_source, V_body, V_G, V_SG | 0 V |

Discretisation: unstructured conforming triangles, P1 elements,
`div(ε ∇φ) = −ρ`, obstacle posed on `u = V_D − φ`:

    u ≥ 0,  r = (Ku − f_u) ≥ 0,  u·r = 0   (silicon, free DOFs only)
    r = 0                                   (oxide, unconstrained)
    f_u = V_D·K·1 − b_phi

---

## 2. Conclusion 1 — Validated

| Item | Evidence |
|---|---|
| Conforming geometry | 7 faces (Si 1, SiO2 2, Gate 2, Shield 2), area sum 3.267 µm² vs analytic 3.267, rel. err 1.4e-16 |
| Material topology | 0 triangles spanning two materials, measured at element centroids |
| Electrode isolation | Gate ∩ Shield = 0 nodes; neither intersects silicon |
| Mesh quality | r5: 298,355 tri, min angle 36.00°, max aspect 2.51, 0 inverted |
| Charge assembly | integrated 2440.3906763934638 vs analytic 2440.3906763934638, rel. err **0.0** |
| Linear FEM operator | R_free = 1.63e-15 (u-form), finite |
| Silicon-only complementarity | primal min u_Si = 0.0; dual min r_Si = −9.74e-13; comp. normalised 2.30e-12; oxide equilibrium 7.08e-12; interface reaction 1.87e-12 |
| Orphan check | 96 unused nodes removed; 3 lay on Dirichlet planes, inspected, not connected to any element |

The linear operator test was run **before** the obstacle was activated, so
assembly correctness is established independently of the complementarity.

---

## 3. Conclusion 2 — Converging

**φ(1.30, 0)** — increments −0.013548 then −0.009504, ratio 1.43 against an
h-ratio of 1.43: clean first-order. Richardson limit ≈ **19.8085 V**.

**|E|(1.30, 0)** — 0.14455 → 0.14286 → 0.14317 MV/cm, a 1.2 % band. The
gradient-recovery spread halves each refinement (±0.00713 → ±0.00362 →
±0.00189), confirming the local gradient converges even where the point value
wobbles.

**Front symmetry** — x = 2.40 asymmetry 0.0010–0.0026 µm across r1–r5; x = 2.35
collapsed 0.0455 → 0.0253 → 0.0009 and stayed ≈0.002. With local refinement the
2.40 asymmetry reached **0.0000**.

**Drain connectivity** — a single component at every resolution, verified
through real elements rather than a distance test: 27 → 42 → 61 → 84 bridging
triangles, 13 → 21 → 31 → 42 active nodes on the drain face, gap 22.1 → 6.4 nm.
Zero active interior nodes on the Dirichlet plane.

---

## 4. Conclusion 3 — Not converged

Global refinement r1 → r5 (16,964 → 298,355 triangles):

| | r1 | r2 | r3 | r4 | r5 |
|---|---|---|---|---|---|
| Active nodes | 34 | 70 | 146 | 285 | 584 |
| Lumped area (µm²) | 0.02118 | 0.02033 | 0.02071 | 0.02133 | 0.02145 |
| P1 zero-set area (µm²) | 0.02060 | 0.01990 | 0.02039 | 0.02108 | 0.02129 |
| bbox y | ±0.177 | ±0.190 | ±0.193 | ±0.200 | ±0.210 |

The P1 area is **non-monotonic** (dip at r2) with a 5.7 % spread and increments
that grow rather than shrink. Both the lumped and the P1 estimator show the same
shape and differ by only 1.2–2.8 %, so the fluctuation is not an estimator
artefact — it is in the discrete neutral set itself.

Local refinement r5 → fA → fB (h_loc 0.0035 → 0.00175, distant mesh unchanged):

| | r5 | fA | fB |
|---|---|---|---|
| Triangles | 298,355 | 290,465 | 328,698 |
| Active nodes | 584 | 2,058 | 8,243 |
| P1 area (µm²) | 0.02129 | 0.02148 | 0.02163 |
| max \|y\| extent | 0.21258 | 0.21466 | 0.21826 |

Area increments 0.00019 then 0.00015 — shrinking, consistent with local control.
But the **outer extent keeps growing** (0.2126 → 0.2147 → 0.2183) at roughly
0.002 µm per halving of the local size. So the area is being driven by local
resolution, while the outermost boundary has not stopped moving.

---

## 5. Conclusion 4 — What the reconstruction established

Two reconstructions were used, because the obstacle makes `u ≥ 0` with a region
of exact zeros, so conventional sign-change marching has nothing to bracket:

1. **Exact discrete zero set** — triangles whose three nodes have u ≤ 1e-9, and
   the boundary of that connected set traced as a polyline. Mesh-aligned.
2. **Positive level sets u = δ**, δ → 0+, giving sub-element positions.

Exact zero-set intersections at fixed cuts (negative / positive y):

```
cut      r5                    fA                    fB
x=2.35   -0.0756/+0.0649      -0.0715/+0.0712      -0.0712/+0.0711
x=2.38   -0.1277/+0.1302      -0.1334/+0.1340      -0.1350/+0.1358
x=2.39   -0.1514/+0.1533      -0.1541/+0.1575      -0.1552/+0.1554
x=2.40   -0.1781/+0.1697      -0.1758/+0.1740      -0.1738/+0.1750
x=2.41   -0.1931/+0.1991      -0.1952/+0.1939      -0.1973/+0.1946
```

**Findings:**

- At **x = 2.35** the zero boundary is stable to 0.0004 µm, while the old
  window metric was still drifting by 0.003 µm and sat ~0.04 µm away. There the
  apparent drift was **measurement error**, not boundary motion.
- At **x = 2.40** the zero boundary *contracts* toward symmetry
  (−0.1781 → −0.1738) while the window metric *expanded* (−0.2092 → −0.2174).
  Opposite trends: the window was contaminated.
- At **x = 2.38 and 2.41** the zero boundary keeps moving outward monotonically
  (2.41: −0.1931 → −0.1952 → −0.1973, ~0.002 µm per halving), and the max
  extent grows. **This survives the corrected measurement**, so it is genuine
  discretisation drift, not extraction error.
- **Level-set sensitivity** confirms the interior field converges: at δ = 1e-3
  and 1e-2 all three meshes agree to four decimals; δ = 1e-12 and 1e-6 coincide
  with the zero-set trace.

**Interpretation.** The P1 obstacle zero set is a union of mesh-aligned
elements. Its boundary rides the last element where u is exactly zero, and as h
decreases that band resolves more of the near-zero region, so the traced
boundary creeps outward toward the true level set at a rate set by the local
gradient. This is a property of the discretisation, not of the measurement.

**Therefore** the P1 front position and neutral area are not certified
converged. The provisional `A ≈ 0.0211 ± 0.0006 µm²` is a descriptive spread of
five discrete solutions, **not** a certified 3 % uncertainty.

---

## 6. Limitations and open items

1. **Oxide offset convention.** The 0.10 µm oxide is a *vertical* offset of the
   profile, not a constant normal thickness. These differ where the profile is
   steep. Not changed during this study and not resolved.
2. **Oxide peak field.** 3.9168 → 4.3464 → 4.7180 MV/cm across r1–r3, location
   pinned to the rounded tip, flipping between symmetric corners. This shows the
   *sampled maximum* is not converged. It does **not** establish a mathematical
   singularity at finite radius, and no extrapolation is warranted.
3. **Neutral-area definition.** Two estimators are reported (lumped nodal area
   and exact P1 zero-set clipping). They agree to ~1 % at the finest grids; the
   reported convergence failure is not an estimator artefact.
4. **Non-nested meshes.** Independently regenerated Gmsh meshes are progressively
   refined but not nested. Acceptable for this study, and recorded because it
   rules out exact order-of-accuracy fitting.
5. **Boundary conditions vs tags.** Physical-group tags and assembled Dirichlet
   rows were audited separately; the tags alone were not treated as evidence.

---

## 7. Disposition

- **P1 is frozen** as the validated reference. No further P1 refinement.
- Stage 2B.2/2B.3 remain **open**: regular field converged, free boundary and
  neutral area not.
- **P2 is justified on evidence**, as a separately versioned solver, with a
  correct higher-order obstacle constraint — non-negative P2 nodal values do not
  imply `u_h ≥ 0` inside a quadratic triangle.
- **Stage 3 (tip-radius sensitivity) stays suspended.**
