# P2.3-S0 — CP-L2 Mesh-Only Pilot Report

**Date:** 2026-10-26
**Scope:** mesh generation and validation only. No OSQP, no MINRES, no F7, no
area postprocessing. Nothing archived has been modified.

## CP-L2.0 — Provenance lock: PASS

```
frozen band SHA256   D2A00BEC…62753B3   matches              yes
design / geometry    e937926 / 1c8d4df  (referenced, unchanged)
manifest             a5b76e2            (referenced, unchanged)
generator            p2_mesh_S0.py, gmsh 4.15.2
L0 / L1 solutions    untouched
```

## CP-L2.1 — Frozen geometry identity: PASS

```
L2 boundary nodes vs L0 canonical polygon   1,928 nodes
  max exact point-to-segment distance        1.110e-16 µm
  p99                                        1.110e-16 µm
Si area (L0 = L1 = L2)                       1.9545436137453363
exact clipped charged area (all levels)      1.5752528820380198
P2 midside deviation from edge midpoint      0.00e+00
```

The polygon is identical to roundoff; no level-dependent arc re-approximation
survives. 913,595 compacted DOFs / 455,772 elements / 915,727 nodes.

## CP-L2.2 — Realized global refinement: PASS except the arc band

Ratios L2/L1 against the target `r = 0.70710678`:

| Region | p50 ratio | p90 ratio | Verdict |
|---|---:|---:|---|
| **arc ≤ 0.01 µm** | **0.7843** | 0.7096 | **+10.9 % off target** |
| band ≤ 0.04 µm | 0.7072 | 0.7002 | pass |
| transition 0.04–0.06 | 0.7002 | 0.7207 | pass |
| tip | 0.7071 | 0.7051 | pass |
| body cut | 0.7071 | 0.7071 | pass |
| background | 0.7071 | 0.7071 | pass |

Realised L2 medians: band 0.002000, tip 0.002000, body cut 0.004250,
background 0.004250 — all exactly on target. **Arc-adjacent elements sit at
0.003606 µm**, i.e. 1.8× coarser than the band target.

**This is the predicted fixed-chord floor.** The canonical arcs carry 65 nodes /
64 chords, chord length

```
l_chord = 2R sin(pi/256) = 0.002454 µm
```

which is longer than the L2 tip target (0.0020 µm). The boundary discretisation
therefore caps refinement near the arc, and the L2/L1 ratio there degrades from
0.7071 to 0.7843 — a ~11 % shortfall exactly where the prediction put it.

**Consequence to weigh before any L2 solve:** the neutral boundary passes near
the trench tip, so an under-resolved arc could bias the L2 area and corrupt the
second increment. This is a design decision, not a numerical error. Per the
authorization, no chord-count change has been made.

## CP-L2.3 — Quality and orientation: PASS

| Metric | L1 | L2 |
|---|---:|---:|
| Minimum angle | 28.92° | **21.96°** |
| Aspect max | 3.564 | 3.564 |
| Aspect p99 | 1.729 | 1.618 |
| Elements aspect > 3 | 60 | 11 |
| Elements aspect > 5 | 0 | 0 |
| Zero-area elements | 0 | 0 |

Orientation is consistent **per face**: Si all one winding (284,230 negative),
SiO₂ the opposite (84,396 positive), Gate and Shield each exactly half/half
because each is two surfaces with opposed normals. No within-face inconsistency.

Minimum angle degrades from 28.9° to 22.0° — still valid, no inverted or
degenerate cells, but the trend is downward and the worst elements should be
located before a solve.

## CP-L2.4 — Donor-load integrity: PASS

```
sum f / rho     1.575252882035949
exact clipped   1.5752528820380198
difference      -2.07e-12   (roundoff)
```

## CP-L2.5 — Resources: PASS for planning

```
nodes 915,727   elements 455,772
DOFs  913,595 total · 638,124 free · 275,471 fixed
charged elements 235,587   Bernstein rows 1,413,522

K   10,490,957 nnz  123.5 MB
Kff  7,306,208 nnz   86.0 MB
C_F  8,463,138 nnz   99.3 MB

assembly 3.4 s · donor 28.4 s · Bernstein 1.5 s · total 35.3 s
assembly-only peak RSS 997.5 MB
physical RAM 24,467 MB · available at start 5,268 MB
```

**Forecast:** L1 solver peak was 998 MB at 459,288 compacted DOFs. L2 is
1.99× that in DOFs and 1.99× in rows, so a linear projection gives
**~2.0 GB** peak for OSQP plus polish. Against the 70 % ceiling of 17,127 MB
this is comfortable; against the 5,268 MB actually available it leaves
~3.2 GB margin. As always the projection is not an upper bound — sparse
fill-in and active-set growth are not linear in DOF count.

## CP-L2.6 — Reproducibility: PASS

```
L2 mesh SHA256   7E6CF12400B0A1308DC23BC7DEC89DA96EF30766784C90D0CB433166523C9136
band SHA256      D2A00BECA0D38F0183669790CBABC8864D6D467C92C0B609AC85BE0ED62753B3
generator        p2_mesh_S0.py (level 3 build call), gmsh 4.15.2
diagnostics      CP-L2_geometry.json, CP-L2_resources.json
```

## Summary table

| Checkpoint | Result | Note |
|---|---|---|
| CP-L2.0 provenance | PASS | band hash matched |
| CP-L2.1 geometry identity | PASS | 1.110e-16 |
| CP-L2.2 realized refinement | **PASS with arc deviation** | +10.9 % p50 near arc |
| CP-L2.3 quality / orientation | PASS | min angle 21.96° |
| CP-L2.4 donor integrity | PASS | −2.07e-12 |
| CP-L2.5 resources | PASS (planning) | ~2.0 GB forecast |
| CP-L2.6 reproducibility | PASS | hashes archived |

**No L2 solve, no F7, no area work has been performed.** The open question for
review is whether the arc-floor deviation is acceptable for the second
increment, or whether the chord policy needs a separate review first.
