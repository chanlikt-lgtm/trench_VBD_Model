#!/usr/bin/env python3
"""Re-validation on the corrected quadrature.

For each resolution:
  - assemble with the corrected QW (sum 0.5)
  - solve the linear P2 problem with the same Dirichlet data as the superseded runs
  - report R_free, potential range, overshoot, min u, Bernstein census
  - rerun the frozen interface audit (edge-matched quadrature, unchanged)
  - save the corrected solution as <tag>_fixed.npz (superseded files untouched)
"""
from __future__ import annotations
import json
import numpy as np
import scipy.sparse.linalg as spla

import p2_audit as A
import p2_interface as PI
from p2_bernstein import B

M = A.M
V_D = A.V_DRAIN


def rabs(recs):
    if not recs:
        return None
    w = np.array([q["ds"] for q in recs])
    ds = np.array([q["d_si"] for q in recs])
    do = np.array([q["d_ox"] for q in recs])
    return float(np.sum(w*np.abs(ds-do))/np.sum(w*(np.abs(ds)+np.abs(do))/2))


def run(fn, tag):
    z = np.load(fn)
    nodes, tris, mat = z["nodes"], z["tris"], z["mat"]
    x = nodes[:, 0]
    K, f = A.assemble(nodes, tris, mat)

    si_el = mat == M.SI
    si_nodes = np.unique(tris[si_el])
    gate = np.unique(tris[mat == M.GATE])
    shield = np.unique(tris[mat == M.SHIELD])
    src = si_nodes[x[si_nodes] <= 1e-12]
    drn = si_nodes[x[si_nodes] >= M.BOTTOM_DEPTH - 1e-9]
    body = si_nodes[x[si_nodes] <= A.BODY_BOTTOM + 1e-12]

    fixed = np.zeros(len(nodes), bool); val = np.zeros(len(nodes))
    for s, v in ((src, A.V_SRC), (body, A.V_SRC), (drn, V_D),
                 (gate, A.V_G), (shield, A.V_SG)):
        fixed[s] = True; val[s] = v
    idx = np.flatnonzero(~fixed)

    phi = val.copy()
    phi[idx] = spla.spsolve(K[idx][:, idx].tocsc(), (f - K @ val)[idx])
    r = K @ phi - f; r[fixed] = 0.0
    R = float(np.linalg.norm(r[idx]) /
              (np.linalg.norm(K[idx][:, idx] @ phi[idx]) + np.linalg.norm(f[idx])))

    u = V_D - phi
    b = u[tris[si_el]] @ B.T
    free_si = si_nodes[~fixed[si_nodes]]
    imin = free_si[np.argmin(u[free_si])]

    recs = PI.interface_flux(nodes, tris, mat, phi)
    wall = [q for q in recs if q["region"] == "wall"]
    arc = [q for q in recs if q["region"] == "arc"]

    out = dict(R_free=R, phi_min=float(phi.min()), phi_max=float(phi.max()),
               overshoot=float((phi[si_nodes]-V_D).max()),
               min_u=float(u[free_si].min()),
               min_u_loc=[float(nodes[imin, 0]), float(nodes[imin, 1])],
               min_b=float(b.min()), negatives=int((b < -1e-12).sum()),
               n_si_elements=int(si_el.sum()),
               elements_affected=int((b < -1e-12).any(axis=1).sum()),
               R_abs_all=rabs(recs), R_abs_wall=rabs(wall), R_abs_arc=rabs(arc),
               interface_samples=len(recs))
    np.savez_compressed(f"{tag}_fixed.npz", nodes=nodes, tris=tris, mat=mat,
                        phi=phi, fixed=fixed)
    print(f"\n[{tag}] corrected")
    for k, v in out.items():
        print(f"  {k:20s} {v}")
    return out


if __name__ == "__main__":
    out = {}
    for tag, fn in (("r3", "p2_r3_solution.npz"), ("r4", "p2_r4_solution.npz")):
        out[tag] = run(fn, tag)
    json.dump(out, open("p2_revalidation.json", "w"), indent=2)
    print("\nwrote p2_revalidation.json")
