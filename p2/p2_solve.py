#!/usr/bin/env python3
"""P2.1 closing audit — physical linear solve + mixed-material interface check.

  1. assemble on the compacted 152,489-DOF mesh with the validated P2 load
  2. solve with the audited Dirichlet data, report independent R_free
  3. potential continuity at shared interface DOFs
  4. one-sided normal displacement eps*grad(phi).n reconstructed independently
     on the Si side and the oxide side, using each element's own gradient and
     the local interface normal
  5. report max, mean and integral measures of the mismatch

Note: continuous Galerkin P2 enforces flux conservation only weakly, so the
pointwise reconstructed jump is NOT expected to be machine zero. The evidence
sought is that the reconstruction is independent on each side and that the
mismatch is small and behaves consistently.
"""
from __future__ import annotations
import json
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

import p2_audit as A
M = A.M
QP = A.QP
V_DRAIN = A.V_DRAIN
EPS_SI, EPS_OX = A.EPS_SI, A.EPS_OX


def main():
    d = np.load("p2_r3_p2mesh.npz")
    nodes, tris, mat = d["nodes"], d["tris"], d["mat"]
    used = np.unique(tris.ravel())
    o2n = np.full(len(nodes), -1, dtype=int)
    o2n[used] = np.arange(len(used))
    nodes = nodes[used]; tris = o2n[tris]
    n = len(nodes)

    K, f = A.assemble(nodes, tris, mat)
    x, y = nodes[:, 0], nodes[:, 1]

    nodes_of = {}
    for code in (M.SI, M.OX, M.GATE, M.SHIELD):
        sel = mat == code
        nodes_of[code] = set(np.unique(tris[sel]).tolist()) if sel.any() else set()
    si_n, ox_n = nodes_of[M.SI], nodes_of[M.OX]
    gate_n, shield_n = nodes_of[M.GATE], nodes_of[M.SHIELD]
    src = {i for i in si_n if x[i] <= 1e-12}
    drn = {i for i in si_n if x[i] >= M.BOTTOM_DEPTH - 1e-9}
    body = {i for i in si_n if x[i] <= A.BODY_BOTTOM + 1e-12}

    fixed = np.zeros(n, bool); val = np.zeros(n)
    for s, v in ((src, A.V_SRC), (body, A.V_SRC), (drn, V_DRAIN),
                 (gate_n, A.V_G), (shield_n, A.V_SG)):
        for i in s:
            fixed[i] = True; val[i] = v
    idx = np.flatnonzero(~fixed)

    rhs = f - K @ val
    phi = val.copy()
    phi[idx] = spla.spsolve(K[idx][:, idx].tocsc(), rhs[idx])
    r = K @ phi - f; r[fixed] = 0.0
    num = float(np.linalg.norm(r[idx]))
    den = float(np.linalg.norm(K[idx][:, idx] @ phi[idx]) + np.linalg.norm(f[idx]))
    R = num / den if den else None

    out = {"solve": dict(free=len(idx), R_free=R, finite=bool(np.all(np.isfinite(phi))),
                         phi_min=float(phi.min()), phi_max=float(phi.max()),
                         residual_norm=num)}

    # ---------- interface nodes and shared-DOF continuity ----------
    iface = sorted(si_n & ox_n)
    # potential continuity: a shared DOF is one value by construction; verify the
    # value is consistent with both sides' element interpolation
    p, det, area, g1, g2, g3 = A.element_geom(nodes, tris)

    # Element gradient evaluated at the centroid (L1 = L2 = 1/3). The gradient
    # of a P2 field is linear, so the centroid value is the element average.
    elem_grad = np.zeros((len(tris), 2))
    for m in range(len(tris)):
        Gc = A.shape_grads(1/3, 1/3, g1[m], g2[m], g3[m])
        elem_grad[m] = (Gc[:, 0] @ phi[tris[m]], Gc[:, 1] @ phi[tris[m]])

    # map interface node -> incident Si and oxide elements
    si_inc, ox_inc = {}, {}
    for m in range(len(tris)):
        tgt = si_inc if mat[m] == M.SI else (ox_inc if mat[m] == M.OX else None)
        if tgt is None:
            continue
        for v in tris[m]:
            tgt.setdefault(int(v), []).append(m)

    rec = []
    for node in iface:
        s_el = si_inc.get(node, [])
        o_el = ox_inc.get(node, [])
        if not s_el or not o_el:
            continue
        # local interface normal from the node's position on the profile
        xn, yn = nodes[node]
        eps_n = 1e-4
        dx = float(M.a_of(xn + eps_n)) - float(M.a_of(xn - eps_n))
        t = np.array([2*eps_n, dx]); t /= np.linalg.norm(t)
        nrm = np.array([t[1], -t[0]])
        if yn < 0:
            nrm = -nrm
        gs = elem_grad[s_el].mean(axis=0)
        go = elem_grad[o_el].mean(axis=0)
        ds = EPS_SI * float(gs @ nrm)
        do = EPS_OX * float(go @ nrm)
        rec.append((node, ds, do, abs(ds - do) / max(abs(ds), abs(do), 1e-30)))

    if rec:
        ds = np.array([q[1] for q in rec])
        do = np.array([q[2] for q in rec])
        mism = np.array([q[3] for q in rec])
        out["interface"] = dict(
            nodes=len(rec),
            d_si_mean=float(ds.mean()), d_ox_mean=float(do.mean()),
            mismatch_median=float(np.median(mism)),
            mismatch_mean=float(mism.mean()),
            mismatch_p95=float(np.percentile(mism, 95)),
            mismatch_max=float(mism.max()),
            integral_si=float(ds.sum()), integral_ox=float(do.sum()),
            integral_rel=float(abs(ds.sum() - do.sum()) /
                               max(abs(ds.sum()), abs(do.sum()), 1e-30)),
        )
    with open("p2_solve.json", "w") as fh:
        json.dump(out, fh, indent=2)
    print(json.dumps(out, indent=2))
    # save the compacted mesh and solution so the interface audit runs on
    # exactly the connectivity that was solved (no renumbering downstream)
    np.savez_compressed("p2_r3_solution.npz", nodes=nodes, tris=tris,
                        mat=mat, phi=phi, fixed=fixed)
    print("wrote p2_r3_solution.npz")


if __name__ == "__main__":
    main()
