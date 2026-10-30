# P2.3-G1 — G1-0 Design Freeze (paired boundary)

**Status:** G1-0 PASS — design amended and frozen. **Mesh-only pilots authorised.**
**Solves, F7, area evaluation:** HOLD (separate G1-4 / G1-5 decisions).
**S0:** CLOSED at `ec4bc80`, read-only.
**Protocol:** `G1_PROTOCOL.md`, amended §4A (paired boundary subdivision).

## 1. Reference

```
level            S0-L1
mesh SHA256      7FA1C69DA097612960A2668335A014894A0DA4D9FF5C96D623E10FCAE0A3D1BB
solution SHA256  47B1E2B44F47AFE8B499EAE36B35C0E9E41AD5D0B7696956AA722658F46C89E9
area             [0.022217431654472784, 0.022217790436636530] µm²
F7               PASS, stationarity 1.652e-10  (commit e1c8a2d)
```

## 2. Frozen geometry definitions

Arc geometry: radius `R = |A_BOT − A_ARC_END| = 0.1 µm`, centres
`(TIP_START, ±A_ARC_END)`, four arcs (upper and lower trench profiles).

**G0 — the S0 polygon (unchanged).** 64 straight chords per 90° arc, 65 vertices at
angles `k·(π/2)/64`, `k = 0..64`.

**C1 — same domain, doubled vertices.** For each chord `P_i → P_{i+1}` insert
`Q_i = (P_i + P_{i+1})/2` **on the straight segment**. Result: 128 collinear
segments, 129 vertices per arc; the polygon is exactly G0.

**C2 — refined polygon.** Same 129 vertex slots, but the inserted point is
`Q_i'` **on the analytic circle at the half-angle**. Result: 128 chords of the
inscribed 128-gon.

**Correspondence.** Vertex index `j = 0..128` per arc is shared by both controls:

```
j even  → original G0 vertex;  displacement Q' − Q = 0 exactly
j odd   → inserted midpoint;   |Q_i' − Q_i| = R[1 − cos(π/256)] = 7.530e-06 µm
```

That displacement is the **old chord sagitta** — a bilateral, checkable prediction.
C2's own sagitta is `R[1 − cos(π/512)] = 1.882e-06 µm`.

## 3. Frozen size field

Defined in **common physical coordinates relative to the analytic circle**, not
recomputed from either polygon:

```
distance to analytic arc <= 0.02 µm        h = 0.0020
0.02 – 0.04 µm                             linear ramp to the L1 value
elsewhere                                  S0-L1 field: band 0.002828,
                                           background 0.006010, tip 0.002404
```

Both controls use this identical field. The chord count is a geometry property,
not a size-field property, and is the only intended difference between C1 and C2.

## 4. Frozen invariants

Unchanged from S0 and not to be altered in G1: conforming straight-sided P2;
exact clipped partial-element donor integration at `x = 0.54 µm` (no centroid
fallback); silicon-only Bernstein positivity obstacle; material tags, body cut,
electrode BCs; checkpoint-first OSQP and scaled MINRES with the corrected
residual sign; full-set F7 thresholds; Bernstein area enclosure at
`AREA_TOL = 1e-13` for later authorised postprocessing.

## 5. Acceptance criteria (pre-registered)

| Check | Threshold |
|---|---|
| C1 vs G0 boundary drift | ≤ 1e-12 µm, exact point-to-segment |
| C1 vs G0 Si area, exact charged area | summation roundoff |
| C1 arc-band realised h ratio (C1/C0), p50 **and** p90 | ≤ 0.75 |
| C1 P2 midside deviation | 0.00e+00 |
| Degenerate / zero-area elements, both controls | 0 |
| C1 vs C2 vertex count and index correspondence | exact, 129 per arc |
| C2 vertex displacement, even / odd indices | 0 / 7.530e-06 µm |
| C2 polygon area vs G0 | reported, not thresholded |
| Donor conservation, both controls | Σf/ρ − exact ≤ 1e-11 |

## 6. Implementation constraint (recorded)

Preserving G0 while doubling the vertex count **cannot** be done by raising the
transfinite node count on the OCC circle — that places the extra nodes on the arc
and silently turns C1 into C2. C1 requires a **polyline-arc geometry builder** in
which the arc is explicit straight segments through the canonical vertices. The
frozen `mesh_stage2b.py` is not modified; G1 uses its own geometry module.

## 7. Artifacts and authorisation

```
p2_3_G1/
  meshes/G1_C1_p2mesh.npz      G1_C2_p2mesh.npz
  G1_C1_mesh.json              G1_C2_mesh.json
  G1-0_freeze.json             (this record, machine-readable)
```

**Authorised:** C1 and C2 pilot mesh generation; geometry, quality, size,
donor-load and provenance diagnostics; assembly-only resource measurement.
**Not authorised:** OSQP, MINRES, F7, neutral-area evaluation; any change to S0
archives; any chord-count or band change outside this document.

Failure of any threshold returns the work to HOLD with evidence preserved — no
automatic policy change.
