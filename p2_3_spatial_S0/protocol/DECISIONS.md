# P2.3-S0 — Design Decisions (checkpoints CP-S0.0 … CP-S0.3)

**Status:** DESIGN ONLY. No mesh generated, no solve run, no production action.
**Protocol:** `P2_3_SPATIAL_CONVERGENCE_PROTOCOL_S0.md` (filed verbatim,
SHA256 `FCCD16A8F19F9D0DC6E57B446CA8A026294C127895A505DAF31BD4066D12D97F`).
**Prepared:** 2026-10-23.

## CP-S0.0 — Archive lock: PASS

```
15ec811 -> 15ec8119f3476f160d15cb57c64cd905c1ab59ca
2478f82 -> 2478f821e8bfb2b7f7840c15070576038d85c7ef
4e5e6d4 -> 4e5e6d4a233a56982b4cfa04d7589f0e52e2ce12
```

All three milestones resolve in the repository. Archived values reproduced
verbatim in §2 of the protocol: the eight 1e-13 endpoints, the two increments,
the ratio enclosure [0.4826, 0.4855], and the B2-X F7 record. The three-mesh
chain is marked CLOSED. No solve is running or queued.

## CP-S0.1 — Scientific target and geometry policy: PASS (decision made)

**Target scalar.** `A_h` = Bernstein-bounded area of `{u ≤ 1e-8}` inside charged
silicon, reported as an interval `[L_h, U_h]`. Front coordinates are a
*separate* diagnostic series and are never folded into the area claim.

**Geometry fork: G-fixed chosen.** All levels share one frozen polygonal domain,
one material partition, one body cut at x = 0.54 µm, and identical interface
polylines. Only the interior triangulation is refined; boundary vertices may be
subdivided *along the same segments* without moving the polygonal boundary.

**Consequence, recorded now:** every limit obtained is a limit for the frozen
polygonal geometry. It is **not** a physical-domain limit until a separate
geometry study bounds the arc/chord approximation error. Measured evidence that
this matters: B0/B1 Si area `1.954559924360` at h≈0.0085 and B2 the same, but
the L1b-vs-B1 comparison earlier showed a `+1.50e-05` µm² drift when band size
changed near the arc — i.e. the discrete area of a *fixed* polygon is not
perfectly mesh-invariant because the arc is chorded. G-fixed keeps the polygon
fixed but does not remove chord error at the boundary; that is the geometry
budget's job.

**Geometry error budget (declared):** to be established by an arc-refinement
study; until then, area differences below ~1e-4 µm² cannot be attributed to the
continuum operator with confidence, which is why the study is designed around
increments ≫ that scale.

## CP-S0.2 — Numerical freeze and fixtures: PASS

Code hashes (SHA256, first 16):

```
p2_audit.py            208419F13A8DC5FE
p2_exact_donor.py      6202E14C9FBE608A
p2_area_bounds.py      CE325BEE9CCECD74
p2_exact_metrics.py    B05A27A552FCE6BB
p2_nested_consistency.py EB3100696F8B97F0
p2_topology_clean.py   BBBE488CEA1A7DAE
p2_mesh_B2.py          B0746D15BE064CF5
p2_band_contour.npz    D2A00BECA0D38F0183669790CBABC8864D6D467C92C0B609AC85BE0ED62753B3
```

Fixture results (`p2_s0_tests.json`):

| Fixture | Result |
|---|---|
| Full element load identity `ρA[0,0,0,1/3,1/3,1/3]` | max err 1.39e-17 |
| Fully body-side element | exactly 0 |
| Cut element: quadrature sum vs clipped area | err 8.33e-17 |
| Cut element: negative vertex components | present (−0.0144) — expected, `N=L(2L−1)` is sign-changing |
| Area, u ≡ 0 | [0.5, 0.5] on A = 0.5 |
| Area, u ≡ +1 | [0, 0] |
| Area, linear manufactured, eps = 0.5 | [0.375, 0.375 + 1.91e-06] |
| Front extractor `2x−3y`: value / gradient | 2.2e-16 / 8.9e-16 |
| Global `Σf/ρ` vs exact charged (B2) | 1.5752691926529643 vs 1.5752691926529394 |

**Two fixture defects were found and corrected during this checkpoint** — the
first run FAILED:

1. the synthetic triangle had all mid-side dofs aliased to node 0, so the
   interpolant was not quadratic;
2. `nodes[tris[0]][:, :3]` sliced columns, handing six points to the clipper.

Both were harness errors, not production defects; the corrected suite passes.
Recording them matters: a checkpoint that has never failed has not been tested.

## CP-S0.3 — Master field, levels, resources: PASS (estimates, not measurements)

**Master field.** `h_k(x,y) = r^k h_0(x,y)`, `r = 2^(−1/2)`, with `h_0` frozen
as: band core 0.0040 µm inside 0.04 µm of the frozen L1b contour, transition to
0.0085 µm by 0.06 µm, background 0.0085 µm, tip 0.0034 µm scaled by the same
`r^k` (no unscaled cap — the study is therefore **global**, not localized).

**Frozen band.** One component, length 1.2492 µm, x ∈ [2.333, 2.420],
y ∈ [−0.575, 0.575], 625 resampled points, SHA256 as above. The band is never
re-derived from a later solution.

| Level | r^k | band h | background h | role |
|---|---:|---:|---:|---|
| S0-L0 | 1.0000 | 0.00400 | 0.00850 | reference |
| S0-L1 | 0.7071 | 0.00283 | 0.00601 | first refinement |
| S0-L2 | 0.5000 | 0.00200 | 0.00425 | conditional |
| S0-L3 | 0.3536 | 0.00141 | 0.00301 | not authorized |

**The archived B1/B2 meshes are not S0 levels.** B1/B2 scale only the band; S0
scales the whole field. No archived area may be inserted as an S0 data point.

**Resource model** (`N_T ∝ h^-2`, anchored on B1 = 114,721 elements at L0-like
sizes; uncertainty ±30 % until a pilot mesh exists):

| Level | est. elements | est. DOFs | est. OSQP RAM | est. solve | est. area @1e-13 |
|---|---:|---:|---:|---:|---:|
| L0 | ~115 k | ~230 k | ~0.7 GB | ~10 min | ~12 min |
| L1 | ~230 k | ~460 k | ~1.6 GB | ~25 min | ~18 min |
| L2 | ~460 k | ~920 k | ~3.5 GB | ~60 min | ~28 min |

Anchors: B2-X 125,579 elements / 252,290 DOFs → 680 MB peak, 5,450 OSQP
iterations, 697 s area at 887 mixed elements. Abort threshold proposed at 70 %
of physical RAM. These are projections, not measurements, and CP-S0.3 requires
a pilot before any production authorization.

## CP-S0.4 — Status: package assembled, review requested

Deliverables: this file, the filed protocol, `manifests.json`, the four
checkpoint records, and `p2_s0_tests.json`. **No mesh generation, no solve, and
no modification of any certified artifact.**

Requested: review and, if accepted, authorization limited to a pilot mesh at
L0 and L1 only. Production solves remain separately gated.
