#!/usr/bin/env python3
"""Controlled P1/P2 physical-equivalence study.

A. Reproduce the archived P1 unconstrained solution independently (read-only on
   the frozen baseline), using P1's own element_material / assemble /
   boundary_sets, and compare the overshoot with the archived 0.21254997 V.

B. Build a labelled P2 variant that matches P1's body treatment exactly:
     - body Si elements (centroid x <= 0.54) carry NO donor charge
     - every node of those elements is pinned at 0  (equipotential, as P1 does
       for all nodes of METAL elements)
     - silicon permittivity retained (P1 gives METAL eps_SI)
   Everything else stays as in the validated P2 solver.

C. Compare the two unconstrained potentials at identical physical coordinates
   near the drain, and compare the deficit u = V_D - phi (minimum and location).

D. Confirm elementwise charge integrals and boundary-node counts.
"""
from __future__ import annotations
import json
import sys
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

sys.path.insert(0, "..")
import conforming_fem as CF
import p2_audit as A

M = A.M
V_D = A.V_DRAIN
BLOCKS = {"Si": CF.SI, "SiO2": CF.OX, "Gate": CF.METAL, "Shield": CF.METAL}


# ---------------- P1 ----------------
def p1_load():
    d = np.load("../r3_mesh.npz")
    nodes = d["nodes"]
    tris, mats = [], []
    for k, v in BLOCKS.items():
        arr = np.asarray(d[k], dtype=int)
        tris.append(arr)
        mats.append(np.full(len(arr), v))
    tris = np.vstack(tris)
    # the archive carries unused nodes (138 isolated free DOFs); compact them
    used = np.unique(tris.ravel())
    o2n = np.full(len(nodes), -1, dtype=int)
    o2n[used] = np.arange(len(used))
    return nodes[used], o2n[tris], np.concatenate(mats)


def p1_area(nodes, tris):
    p = nodes[tris]
    det = ((p[:, 1, 0]-p[:, 0, 0])*(p[:, 2, 1]-p[:, 0, 1]) -
           (p[:, 2, 0]-p[:, 0, 0])*(p[:, 1, 1]-p[:, 0, 1]))
    return 0.5*np.abs(det)


def p1_unconstrained():
    nodes, tris, _ = p1_load()
    mat = CF.element_material(nodes, tris, 0.10)
    K, f = CF.assemble(nodes, tris, mat)
    fixed, _ = CF.boundary_sets(nodes, tris, mat, 0.10)
    # audited P1 u-formulation: u = V_D - phi, so u_D = V_D at grounded
    # electrodes and 0 at the drain, with rhs = -f - K u_D (bc_audit_u.json).
    drn_nodes = nodes[:, 0] >= CF.BOTTOM_DEPTH - 1e-12
    # Dirichlet values on Dirichlet nodes ONLY; free entries must stay zero
    # or K @ val double-counts them in the reduced right-hand side.
    val = np.zeros(len(nodes))
    val[fixed] = V_D
    val[drn_nodes] = 0.0
    free = ~fixed
    idx = np.flatnonzero(free)
    u = val.copy()
    u[idx] = spla.spsolve(K[idx][:, idx].tocsc(), (-f - K @ val)[idx])
    phi = V_D - u
    area = p1_area(nodes, tris)
    si = mat == CF.SI
    si_nodes = np.unique(tris[si])
    free_si = si_nodes[free[si_nodes]]
    print(f"[P1] nodes={len(nodes)} elements={len(tris)} Si={int(si.sum())} "
          f"charged Si area={area[si].sum():.9f} charge={CF.RHO_SI*area[si].sum():.6f}")
    print(f"[P1] fixed nodes={int(fixed.sum())}  max(phi)-V_D over Si = "
          f"{(phi[si_nodes]-V_D).max():+.9f}")
    print(f"[P1] min u over free Si = {u[free_si].min():+.9f} V at "
          f"(x={nodes[free_si[np.argmin(u[free_si])],0]:.5f}, "
          f"y={nodes[free_si[np.argmin(u[free_si])],1]:.5f})")
    print(f"[P1] archived overshoot = 0.21254997312667828")
    np.savez_compressed("p1_r3_uncon.npz", nodes=nodes, tris=tris, mat=mat,
                        phi=phi, u=u, fixed=fixed)
    return nodes, tris, mat, phi, u, fixed, area


# ---------------- P2 variant ----------------
def p2_area(nodes, tris):
    return np.abs(A.element_geom(nodes, tris)[2])


def assemble_p2(nodes, tris, mat, charge):
    n = len(nodes)
    _, det, area, g1, g2, g3 = A.element_geom(nodes, tris)
    eps = np.where(mat == M.OX, A.EPS_OX, A.EPS_SI)
    rows = np.zeros((len(tris), 36), dtype=np.int64)
    cols = np.zeros((len(tris), 36), dtype=np.int64)
    vals = np.zeros((len(tris), 36))
    f = np.zeros(n)
    for m in range(len(tris)):
        ke = np.zeros((6, 6))
        for q in range(3):
            L1, L2 = A.QP[q]
            G = A.shape_grads(L1, L2, g1[m], g2[m], g3[m])
            ke += A.QW[q] * abs(det[m]) * eps[m] * (G @ G.T)
        ids = tris[m]
        rows[m] = np.repeat(ids, 6)
        cols[m] = np.tile(ids, 6)
        vals[m] = ke.ravel()
        if charge[m]:
            f[ids[3:]] += A.RHO * abs(area[m]) / 3.0
    K = sp.coo_matrix((vals.ravel(), (rows.ravel(), cols.ravel())),
                      shape=(n, n)).tocsr()
    return K, f


def p2_variant():
    z = np.load("p2_r3_solution.npz")
    nodes, tris, mat = z["nodes"], z["tris"], z["mat"]
    x = nodes[:, 0]
    cxs = nodes[tris][:, :3, 0].mean(axis=1)
    si = mat == M.SI
    body = si & (cxs <= A.BODY_BOTTOM + 1e-12)
    charge = si & ~body

    K, f = assemble_p2(nodes, tris, mat, charge)
    area = p2_area(nodes, tris)

    fixed = np.zeros(len(nodes), bool)
    val = np.zeros(len(nodes))
    body_nodes = np.unique(tris[body])
    fixed[body_nodes] = True; val[body_nodes] = 0.0
    si_nodes = np.unique(tris[si])
    src = si_nodes[x[si_nodes] <= 1e-12]
    # P1 pins the ENTIRE drain plane (Si and oxide), not just silicon nodes
    drn = np.flatnonzero(x >= M.BOTTOM_DEPTH - 1e-9)
    fixed[src] = True; val[src] = 0.0
    fixed[drn] = True; val[drn] = V_D
    gate_nodes = np.unique(tris[mat == M.GATE])
    shield_nodes = np.unique(tris[mat == M.SHIELD])
    fixed[gate_nodes] = True; val[gate_nodes] = 0.0
    fixed[shield_nodes] = True; val[shield_nodes] = 0.0

    idx = np.flatnonzero(~fixed)
    u = val.copy()
    u[idx] = spla.spsolve(K[idx][:, idx].tocsc(), (f - K @ val)[idx])
    # this solve uses the potential convention (phi_D = 0 grounded,
    # V_D at the drain), so u IS phi; do not flip it again
    phi = u
    free_si = si_nodes[~fixed[si_nodes]]
    print(f"\n[P2v] elements={len(tris)} body elements={int(body.sum())} "
          f"charged Si area={area[charge].sum():.9f} charge={A.RHO*area[charge].sum():.6f}")
    print(f"[P2v] fixed nodes={int(fixed.sum())}  max(phi)-V_D over Si = "
          f"{(phi[si_nodes]-V_D).max():+.9f}")
    print(f"[P2v] min u over free Si = {u[free_si].min():+.9f} V at "
          f"(x={nodes[free_si[np.argmin(u[free_si])],0]:.5f}, "
          f"y={nodes[free_si[np.argmin(u[free_si])],1]:.5f})")
    np.savez_compressed("p2_r3_variant.npz", nodes=nodes, tris=tris, mat=mat,
                        phi=phi, u=u, fixed=fixed, body=body)
    return nodes, tris, mat, phi, u, fixed, area, body, charge


# ---------------- probes ----------------
def locate(nodes, tris, pts):
    Aa = nodes[tris[:, 0]]; Bb = nodes[tris[:, 1]]; Cc = nodes[tris[:, 2]]
    det = ((Bb[:, 0]-Aa[:, 0])*(Cc[:, 1]-Aa[:, 1]) -
           (Cc[:, 0]-Aa[:, 0])*(Bb[:, 1]-Aa[:, 1]))
    el = np.full(len(pts), -1, int)
    bar = np.zeros((len(pts), 2))
    for k, p in enumerate(pts):
        L1 = ((Bb[:, 1]-Cc[:, 1])*(p[0]-Cc[:, 0]) +
              (Cc[:, 0]-Bb[:, 0])*(p[1]-Cc[:, 1]))/det
        L2 = ((Cc[:, 1]-Aa[:, 1])*(p[0]-Cc[:, 0]) +
              (Aa[:, 0]-Cc[:, 0])*(p[1]-Cc[:, 1]))/det
        L3 = 1.0 - L1 - L2
        m = (L1 >= -1e-9) & (L2 >= -1e-9) & (L3 >= -1e-9)
        if m.any():
            j = int(np.flatnonzero(m)[0])
            el[k] = j; bar[k] = (L1[j], L2[j])
    return el, bar


def eval_p1(nodes, tris, phi, pts):
    el, bar = locate(nodes, tris, pts)
    out = np.full(len(pts), np.nan)
    for k in range(len(pts)):
        if el[k] < 0:
            continue
        L1, L2 = bar[k]
        ids = tris[el[k]]
        out[k] = L1*phi[ids[0]] + L2*phi[ids[1]] + (1-L1-L2)*phi[ids[2]]
    return out, el


def eval_p2(nodes, tris, phi, pts):
    el, bar = locate(nodes, tris, pts)
    out = np.full(len(pts), np.nan)
    for k in range(len(pts)):
        if el[k] < 0:
            continue
        L1, L2 = bar[k]
        out[k] = float(np.dot(phi[tris[el[k]]], A.shape_vals(L1, L2)))
    return out, el


def main():
    n1, t1, m1, f1, u1, fx1, a1 = p1_unconstrained()
    n2, t2, m2, f2, u2, fx2, a2, body, charge = p2_variant()

    xs = np.linspace(2.25, 2.415, 12)
    ys = np.array([-0.30, -0.20, -0.10, 0.0, 0.10, 0.20, 0.30])
    pts = np.array([[x, y] for x in xs for y in ys])

    p1, e1 = eval_p1(n1, t1, f1, pts)
    p2, e2 = eval_p2(n2, t2, f2, pts)
    ok = ~np.isnan(p1) & ~np.isnan(p2)
    d = p2[ok] - p1[ok]
    print(f"\n[C] probes={len(pts)} matched={int(ok.sum())}")
    print(f"    dphi (P2v - P1): max|.|={np.abs(d).max():.6e}  rms={np.sqrt(np.mean(d*d)):.6e}  "
          f"mean={d.mean():+.6e}")
    print(f"    max |phi| at probes: P1={np.nanmax(np.abs(p1)):.6f} P2v={np.nanmax(np.abs(p2)):.6f}")

    print(f"\n[D] charge: P1={CF.RHO_SI*a1[m1 == CF.SI].sum():.6f}  "
          f"P2v={A.RHO*a2[charge].sum():.6f}")
    print(f"    charged Si area: P1={a1[m1 == CF.SI].sum():.9f}  "
          f"P2v={a2[charge].sum():.9f}")
    print(f"    fixed nodes: P1={int(fx1.sum())}  P2v={int(fx2.sum())}")

    out = dict(p1=dict(overshoot=float((f1[np.unique(t1[m1 == CF.SI])]-V_D).max()),
                       min_u=float(u1[np.unique(t1[m1 == CF.SI])[~fx1[np.unique(t1[m1 == CF.SI])]]].min()),
                       charge=float(CF.RHO_SI*a1[m1 == CF.SI].sum()),
                       area=float(a1[m1 == CF.SI].sum())),
               p2v=dict(overshoot=float((f2[np.unique(t2[m2 == M.SI])]-V_D).max()),
                        min_u=float(u2[np.unique(t2[m2 == M.SI])[~fx2[np.unique(t2[m2 == M.SI])]]].min()),
                        charge=float(A.RHO*a2[charge].sum()),
                        area=float(a2[charge].sum())),
               probes=dict(n=int(ok.sum()), max_abs=float(np.abs(d).max()),
                           rms=float(np.sqrt(np.mean(d*d))), mean=float(d.mean())))
    json.dump(out, open("p1_p2_equiv.json", "w"), indent=2)
    print("wrote p1_p2_equiv.json")


if __name__ == "__main__":
    main()
