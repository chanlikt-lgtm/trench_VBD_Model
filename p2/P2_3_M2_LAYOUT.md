# P2.3 — M2 Mesh-Layout Experiment (fixed core, wider transition)

**Record date:** 2026-10-21
**Question:** does changing mesh layout at approximately fixed local resolution
produce front and area shifts comparable to the L1b↔L1 differences?

## 1. Design

| Parameter | L1b (reference) | M2 (test) |
|---|---|---|
| Window A realised core | 0.00524 µm | 0.00526 µm (+0.4 %) |
| Window B realised core | 0.00527 µm | 0.00543 µm (+3.0 %) |
| Transition `Thickness` | 0.012 | **0.030** |
| Connectivity | original | **regenerated** |
| Physics / charge rule / BCs | frozen | identical |
| Solver | checkpoint OSQP + scaled MINRES | identical |
| Acceptance | full-set F7 | full-set F7 |

Mesh: 427,047 nodes, 212,838 elements, min angle 37.11°, aspect 2.52,
charged Si area 1.575326220 µm² (L1b 1.575280763).

## 2. F7 gate — PASS

```
active rows        11,169
stationarity       3.430e-10
complementarity    5.035e-17
min slack         -2.267e-14
violations         0
lambda min        +0.000e+00
oxide residual     4.367e-13
overshoot         +1.148e-14
GATE              PASS
```

OSQP `solved`, 14,775 iterations, obj `-9.932317e+04`, pri `1.406e-08`,
dua `1.627e-09`; checkpoint reload verified.

## 3. Bounded area — layout-insensitive

| Mesh | Lower | Upper |
|---|---:|---:|
| L1b | 0.021842101270 | 0.021845749070 |
| **M2** | **0.021839484478** | **0.021843096605** |
| L1 | 0.022071011922 | 0.022075109085 |

Interval-bounded differences:

```
L1b -> M2   [-6.264592e-06, +9.953350e-07]     straddles zero
M2  -> L1   [+2.279153e-04, +2.356246e-04]
L1b -> L1   [+2.252629e-04, +2.330078e-04]
```

Changing layout at fixed core moves the area by at most **6.3e-6 µm²**,
whereas the refinement step moves it by **2.3e-4 µm²** — a factor of ~37.
The neutral area is controlled by local resolution, not by layout.

## 4. Front coordinates — layout-sensitive

Lower edges at ε = 1e-8:

```
x = 2.35   L1b -0.0774   M2 -0.0731   L1 -0.0719
           |M2 - L1b| = 0.0043      |L1 - L1b| = 0.0055   (0.78x)

x = 2.41   L1b -0.2014   M2 -0.1922   L1 -0.1996
           |M2 - L1b| = 0.0092      |L1 - L1b| = 0.0018   (5.1x)
```

At x = 2.35 the layout change reproduces 78 % of the refinement-induced
shift. At x = 2.41 it produces **five times** the refinement-induced shift.

### 4.1 The B-core drift does not explain it

M2's window B core is 3 % coarser. The observed resolution sensitivity
between L1b and L1 is 0.0018 µm for a 32 % change in h, i.e. ≈ 0.0056 µm per
100 % of h. A 3 % change would give ≈ 1.7e-4 µm — two orders below the
0.0092 µm measured. The x = 2.41 shift is therefore layout-driven, not
resolution-driven.

## 5. Conditioning

| Cut | endpoint | u (V) | ∂u/∂y (V/µm) | δy per 1e-8 V |
|---|---|---:|---:|---:|
| 2.35 | -0.0731 | 2.5e-08 | -4.09e-04 | 2.4e-05 |
| 2.35 | +0.0708 | 2.8e-07 | +2.83e-03 | 3.5e-06 |
| 2.41 | -0.1922 | 1.3e-08 | -1.79e-04 | 5.6e-05 |
| 2.41 | +0.1948 | 1.4e-08 | +3.46e-04 | 2.9e-05 |

Signs are physical (negative at lower boundaries, positive at upper). Slopes
are shallower than L1b's (|∂u/∂y| 1.8e-4 … 2.8e-3 versus 4.8e-4 … 3.3e-2),
so M2's front coordinates are *more* amplified than L1b's — consistent with
layout, not resolution, controlling the extracted front.

## 6. Conclusion

1. **Area:** layout-invariant to within 6e-6 µm². The monotone area trend is a
   resolution effect and remains non-convergent (increments still growing,
   see the L1b record).
2. **Front:** layout is a first-order control. At fixed core, changing the
   transition width and connectivity moves the x = 2.35 edge by 0.0043 µm and
   the x = 2.41 edge by 0.0092 µm.
3. Therefore the nonmonotone front at x = 2.35 reported across r3→L1 is
   substantially a mesh-layout/conditioning artefact, not evidence of a
   physical reversal. The area sequence does not share that artefact.
4. P2.3 remains open: no converged area, no converged front. But the two
   statistics now have identified, different sensitivities — area to
   resolution, front to layout — which was not established before.

## 7. Artifacts

- `p2_M2_p2mesh.npz`, `p2_M2_mesh.json`
- `p2_M2_osqp_checkpoint.npz`, `p2_M2_osqp.json`
- `p2_M2_scaled.npz`, `p2_M2_scaled.json`
- `p2_M2_post.json`
- `p2_mesh_M2.py`, `p2_M2_osqp.py`, `p2_M2_polish.py`, `p2_M2_post.py`
