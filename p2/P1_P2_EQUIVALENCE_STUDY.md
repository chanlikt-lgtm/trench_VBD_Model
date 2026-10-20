# P1/P2 controlled physical-equivalence study

**Date:** 19 October 2026
**Status:** charge-support difference ruled out as the cause of the activation discrepancy; element formulation isolated; both discretizations **not converged** in the tip-to-drain region; P2.3 remains blocked; P1 frozen; Stage 3 suspended.

## 1. Reproduction fixes (all in the reproduction harness, none in frozen code)

The archived P1 mesh carries **138 isolated free nodes** (zero diagonal, singleton graph components). `K_FF` is exactly singular until they are compacted away; no degenerate elements exist (min |det| = 1.46e-05). After compaction: 38,260 nodes, 75,970 elements, 7,450 fixed nodes.

Two convention errors in my own harness were found and corrected:

- P1's audited obstacle variable is `u = V_D - phi` with `u_D = V_D` at grounded electrodes and `0` at the drain, and right-hand side `-f - K u_D` (`bc_audit_u.json`). `conforming_fem.boundary_sets` returns potential-convention values, so using it directly solves for `phi`, not `u`.
- `K @ val` must use a vector non-zero **only on Dirichlet nodes**; filling free entries with `V_D` corrupts the reduced right-hand side.

## 2. P1 reproduction (archived r3 mesh)

| Quantity | Value |
|---|---:|
| charged Si area | 1.585441860 µm² |
| donor charge | 1979.525384 |
| max(phi) − V_D over Si | **+0.232581090 V** |
| min u over free Si | −0.232581090 V at (2.35115, −0.00261) |
| archived overshoot | 0.212549973 V |

Same sign and scale; 9.4 % larger. The archived mesh's geometric Si area is 1.964871 µm² versus 1.954558 µm² in the frozen record, so the archived field was produced on a different mesh. The reproduction is formulation-faithful, not bit-identical.

## 3. P2 variant matching P1's physics

Body silicon elements (centroid x ≤ 0.54): donor charge removed, **every node of those elements pinned at 0** (equipotential, as P1 does for all METAL-element nodes), silicon permittivity retained. Drain plane pinned across all materials, matching P1.

| Quantity | P1 | P2 variant |
|---|---:|---:|
| charged Si area | 1.585441860 | 1.575120976 |
| donor charge | 1979.525384 | 1966.639102 |
| max(phi) − V_D over Si | +0.232581090 | **0.000000000** |
| min u over free Si | −0.232581090 | **+0.044884180** |
| location of min u | (2.35115, −0.00261) | (0.54033, 0.32448) |

Charged areas differ by 0.65 % from mesh geometry; on a common mesh the charges are identical (see §4).

## 4. Same-mesh test — element order isolated

P1 linear elements assembled on the **P2 mesh's own vertex triangles**, with identical material tags, body treatment, charge support (1966.639102), and Dirichlet sets:

| | P1 linear, P2 mesh | P2 quadratic, P2 mesh |
|---|---:|---:|
| donor charge | 1966.639102 | 1966.639102 |
| max(phi) − V_D | **+0.220543445** | **0.000000000** |
| min u free Si | −0.220543445 | +0.044884180 |

**Charge support is therefore ruled out as the cause.** With everything else identical, the linear formulation overshoots by 0.22 V and the quadratic does not.

## 5. Field comparison at identical coordinates (84 probes, x 2.25–2.415)

| | P2 variant vs P1 (own meshes) | P2 quadratic vs P1 linear (same mesh) |
|---|---:|---:|
| max abs dphi | 3.555039 V | 3.527556 V |
| rms dphi | 2.084089 V | 2.060222 V |
| mean dphi | −1.820257 V | −1.799232 V |
| near-drain (x > 2.38) max abs | 0.846 V (from same-mesh run) | 0.840587 V |

The field difference is an order of magnitude larger than the 0.22 V overshoot that decides activation. The two discretizations are **not converged** in the tip-to-drain region; the 79.51° corner at x = 2.15 is the likely amplifier.

## 6. Conclusion

1. The P1/P2 activation discrepancy is **not** caused by the body charge-support difference.
2. It is caused by the **element formulation** (linear + vertex load versus quadratic + midside load), which is the only remaining difference on the common mesh.
3. Neither activation state can be trusted yet: the two fields differ by ~2 V rms where the decision is made, far exceeding the 0.22 V margin/overshoot.
4. P1's overshoot is robust at ~0.21–0.23 V across three meshes; P2's margin is positive at both resolutions (+0.0797, +0.0666 V) and shrinking.

## 7. Next steps

1. Extend the same-mesh test to r4 (and r5): P1 linear on the P2 vertex triangles vs P2 quadratic, at fixed probes.
2. Require the probe rms difference to **decrease with h** before interpreting either activation state.
3. Quantify the x = 2.15 corner's contribution (probes inside/outside a corner neighborhood) and consider excluding it from the convergence statistic.
4. No P2 trench obstacle solve and no further refinement of the obstacle until the unconstrained fields converge.
