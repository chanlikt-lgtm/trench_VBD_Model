# P2.3 — First physical trench obstacle solve (P1-matched model)

**Date:** 19 October 2026
**Status:** solved; neutral set nonempty; KKT stationarity at 1.9e-5 (solver tolerance, not machine precision); P1 frozen; Stage 3 suspended.

## 1. Model

- Corrected P2 operator (QW sums to 0.5).
- P1-matched body: silicon elements with centroid x ≤ 0.54 carry no donor charge and all their nodes are pinned at u_D = V_D (equipotential).
- Gate and shield pinned at V_D; source pinned at V_D; entire drain plane pinned at 0.
- Bernstein b ≥ 0 on charged silicon elements only (50,037 elements).
- u = V_D − φ.

## 2. Problem size and solver

| Quantity | Value |
|---|---:|
| DOFs | 152,489 |
| Free DOFs | 123,752 |
| Bernstein constraint rows | 300,222 |
| Charged Si elements | 50,037 |
| Body elements (uncharged, pinned) | 5,573 |
| Total donor charge | 1,966.639102 |
| OSQP status / iterations | solved / 9,825 |
| Objective | −6.4635886179e+04 |

## 3. KKT diagnostics

| Residual | Value |
|---|---:|
| min slack | −3.459e-10 |
| Active constraints | 2,229 / 300,222 |
| λ min / max | −2.804e-17 / 3.074e-04 |
| Stationarity (relative) | 1.933e-05 |
| Complementarity (max) | 1.383e-12 |

Primal feasibility and complementarity are at solver precision. Stationarity is 1.9e-5 relative — acceptable for a 124k-variable QP at OSQP's tolerance, but **not** at the 1e-14 level of the manufactured test. A KKT polish on the identified active set is the next refinement.

## 4. Solution

| Quantity | Value |
|---|---:|
| φ range | [0.000000, 34.274209] |
| Overshoot | +1.080e-11 V |
| min u | −1.080e-11 V |
| min u free Si location | (2.34372, 0.00776) |

The obstacle enforces φ ≤ V_D to machine precision.

## 5. Neutral set

| Definition | Area (µm²) | Elements |
|---|---:|---:|
| Bernstein-active / any node \|u\| < 1e-8 | 0.035064617 | 480 |
| All six nodes \|u\| < 1e-8 | **0.020917675** | 291 |

P1 reference: active area 0.020331 µm² (70 active nodes), charged Si 1.585442 µm².

The fully-neutral P2 area (0.020918) is **2.9 % above** P1's 0.020331. The any-node figure (0.035065) is inflated by fixed drain nodes (u = 0) in adjacent elements and should not be compared directly with P1's active-node definition.

Clamped region: Si nodes with |u| < 1e-6 span x 2.3320–2.4200, y ±0.5750. The upstream edge of the neutral region is x ≈ 2.332.

## 6. Next

1. KKT polish on the 2,229-row active set (equality-constrained solve) to bring stationarity to solver precision.
2. Neutral-set map and front location against P1's archived active nodes.
3. Electric-field comparison at fixed probes.
