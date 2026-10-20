#!/usr/bin/env python3
"""P2.1 second resolution — build, solve, and audit with the FROZEN diagnostic.

Identical geometry builder, identical material tags, identical electrode data,
identical assembly, and the identical interface routine used on r3
(p2_interface.interface_flux). The only thing that changes is h.

D_ref is read from p2_nearzero.json and NOT recomputed, so the near-zero
classification is bit-identical between resolutions.
"""
from __future__ import annotations
import json
import numpy as np
import scipy.sparse.linalg as spla

import p2_build as B
import p2_audit as A
import p2_interface as PI

M = A.M
EPS_SI, EPS_OX = A.EPS_SI, A.EPS_OX

H, H_TIP = 0.0100, 0.0040
PREFIX = "p2_r4"


def stats(recs, dref):
    if not recs:
        return {}
    w = np.array([q["ds"] for q in recs])
    ds = np.array([q["d_si"] for q in recs])
    do = np.array([q["d_ox"] for q in recs])
    dscale = (np.abs(ds) + np.abs(do)) / 2
    den = np.sum(w * (np.abs(ds) + np.abs(do)) / 2)
    rel = np.abs(ds - do) / np.maximum(2 * dscale, 1e-30)
    nz = dscale < 1e-3 * dref
    keep = ~nz
    return dict(n=int(len(recs)),
                R_abs=float(np.sum(w * np.abs(ds - do)) / den),
                R_signed=float(abs(np.sum(w * (ds - do))) / den),
                rel_median=float(np.median(rel[keep])) if keep.any() else None,
                rel_mean=float(np.mean(rel[keep])) if keep.any() else None,
                rel_p95=float(np.percentile(rel[keep], 95)) if keep.any() else None,
                rel_max=float(np.max(rel[keep])) if keep.any() else None,
                near_zero=int(nz.sum()),
                same_sign=int((np.sign(ds) == np.sign(do)).sum()),
                opp_sign=int((np.sign(ds) != np.sign(do)).sum()),
                mean_d_si=float(ds.mean()), mean_d_ox=float(do.mean()))


def main():
    with open("p2_nearzero.json") as fh:
        dref = json.load(fh)["reference"]["d_ref"]
    print(f"[frozen] D_ref = {dref:.6f}  h={H} h_tip={H_TIP}")

    nodes, tris, mat, counts = B.build_p2(H, H_TIP, PREFIX)
    print(f"[mesh] elements={len(tris)} nodes_raw={len(nodes)} counts={counts}")
    for code, name in ((M.SI, "Si"), (M.OX, "SiO2"), (M.GATE, "gate"), (M.SHIELD, "shield")):
        print(f"       {name:6s} {int((mat == code).sum())}")

    used = np.unique(tris.ravel())
    o2n = np.full(len(nodes), -1, dtype=int)
    o2n[used] = np.arange(len(used))
    nodes_c, tris_c = nodes[used], o2n[tris]

    K, f = A.assemble(nodes_c, tris_c, mat)
    x, y = nodes_c[:, 0], nodes_c[:, 1]
    n = len(nodes_c)

    nodes_of = {}
    for code in (M.SI, M.OX, M.GATE, M.SHIELD):
        sel = mat == code
        nodes_of[code] = set(np.unique(tris_c[sel]).tolist()) if sel.any() else set()
    si_n, ox_n = nodes_of[M.SI], nodes_of[M.OX]
    src = {i for i in si_n if x[i] <= 1e-12}
    drn = {i for i in si_n if x[i] >= M.BOTTOM_DEPTH - 1e-9}
    body = {i for i in si_n if x[i] <= A.BODY_BOTTOM + 1e-12}

    fixed = np.zeros(n, bool); val = np.zeros(n)
    for s, v in ((src, A.V_SRC), (body, A.V_SRC), (drn, B.V_DRAIN),
                 (nodes_of[M.GATE], A.V_G), (nodes_of[M.SHIELD], A.V_SG)):
        for i in s:
            fixed[i] = True; val[i] = v
    idx = np.flatnonzero(~fixed)

    rhs = f - K @ val
    phi = val.copy()
    phi[idx] = spla.spsolve(K[idx][:, idx].tocsc(), rhs[idx])
    r = K @ phi - f; r[fixed] = 0.0
    R = float(np.linalg.norm(r[idx]) /
              (np.linalg.norm(K[idx][:, idx] @ phi[idx]) + np.linalg.norm(f[idx])))
    print(f"[solve] dofs={n} free={len(idx)} R_free={R:.3e} finite={np.all(np.isfinite(phi))} "
          f"phi=[{phi.min():.6f},{phi.max():.6f}]")

    np.savez_compressed(f"{PREFIX}_solution.npz", nodes=nodes_c, tris=tris_c,
                        mat=mat, phi=phi, fixed=fixed)

    recs = PI.interface_flux(nodes_c, tris_c, mat, phi)
    wall = [q for q in recs if q["region"] == "wall"]
    arc = [q for q in recs if q["region"] == "arc"]
    out = dict(h=H, h_tip=H_TIP, elements=int(len(tris)), dofs=int(n),
               R_free=R, phi_min=float(phi.min()), phi_max=float(phi.max()),
               d_ref=dref, all=stats(recs, dref), wall=stats(wall, dref),
               arc=stats(arc, dref))
    print(json.dumps({k: v for k, v in out.items() if k != "all"}, indent=2))
    print("all:", json.dumps(out["all"], indent=2))
    with open(f"{PREFIX}_audit.json", "w") as fh:
        json.dump(out, fh, indent=2)
    print(f"wrote {PREFIX}_solution.npz and {PREFIX}_audit.json")


if __name__ == "__main__":
    main()
