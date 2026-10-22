#!/usr/bin/env python3
"""M2 Step B/C: scaled matrix-free MINRES KKT polish from the OSQP checkpoint.

Active set is seeded from the OSQP dual (lam = -sK*y >= 0) as in the L1
pipeline, then repaired by bulk updates. No dense independence filter: the
Bernstein rows are element-local and linearly independent by construction.
Acceptance is computed on the ORIGINAL unscaled system.
"""
from __future__ import annotations
import json
import time
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

import p2_audit as A
from p2_bernstein import B
from p2_nested_consistency import assemble_p2
from p2_exact_donor import exact_load

M = A.M
V_D = A.V_DRAIN


def build():
    z = np.load("p2_B1_p2mesh.npz")
    nodes, tris, mat = z["nodes"], z["tris"], z["mat"]
    used = np.unique(tris.ravel())
    o2n = np.full(len(nodes), -1, dtype=int)
    o2n[used] = np.arange(len(used))
    nodes = nodes[used]; tris = o2n[tris]
    x = nodes[:, 0]
    cxs = nodes[tris][:, :3, 0].mean(axis=1)
    si = mat == M.SI
    body = si & (cxs <= A.BODY_BOTTOM + 1e-12)
    charged = si & ~body
    K, _f0 = assemble_p2(nodes, tris, mat, np.zeros(len(tris), bool))
    f, _ftot = exact_load(nodes, tris, mat)
    fixed = np.zeros(len(nodes), bool); val = np.zeros(len(nodes))
    si_nodes = np.unique(tris[si])
    gate = np.unique(tris[mat == M.GATE]); shield = np.unique(tris[mat == M.SHIELD])
    src = si_nodes[x[si_nodes] <= 1e-12]; drn = np.flatnonzero(x >= M.BOTTOM_DEPTH - 1e-9)
    for s in (np.unique(tris[body]), gate, shield, src):
        fixed[s] = True
    val[fixed] = V_D
    fixed[drn] = True; val[drn] = 0.0
    idx = np.flatnonzero(~fixed)
    cons = np.flatnonzero(charged); nC = 6*len(cons)
    C = sp.coo_matrix((np.tile(B.ravel(), len(cons)),
                       (np.repeat(np.arange(nC), 6), np.tile(tris[cons], (1, 6)).ravel())),
                      shape=(nC, len(nodes))).tocsr()
    C_F = C[:, idx].tocsc()
    cD = C[:, fixed] @ val[fixed]
    Kff = K[idx][:, idx].tocsc()
    g = -f[idx] - (K[idx][:, fixed] @ val[fixed])
    return dict(nodes=nodes, tris=tris, mat=mat, si=si, fixed=fixed, val=val,
                idx=idx, C_F=C_F, cD=cD, K=K, Kff=Kff, g=g, f=f, nfree=len(idx))


def make_op(Kff, CW, n):
    dK = np.maximum(Kff.diagonal(), 1e-12)
    du = 1.0/np.sqrt(dK)
    CWu = CW.multiply(du[None, :])
    rn = np.sqrt(np.asarray(CWu.multiply(CWu).sum(axis=1)).ravel())
    rn = np.maximum(rn, 1e-12)
    dl = 1.0/rn
    D = np.concatenate([du, dl])

    def Hmv(zt):
        z = D*zt
        uu, ll = z[:n], z[n:]
        return D*np.concatenate([Kff@uu - CW.T@ll, -CW@uu])

    m = CW.shape[0]
    return spla.LinearOperator((n+m, n+m), matvec=Hmv, dtype=float), D


def polish(d, W, u, lam, rtol=1e-8, maxiter=8000, nref=6, tag=""):
    Kff, C_F, cD, g = d["Kff"], d["C_F"], d["cD"], d["g"]
    n = d["nfree"]; m = len(W)
    CW = C_F[W]
    t0 = time.time()
    Hop, D = make_op(Kff, CW, n)
    b = np.concatenate([g, cD[W]])
    zt0 = np.concatenate([u/D[:n], lam/D[n:]])
    print(f"  {tag} MINRES start |W|={m} rtol={rtol:g}", flush=True)
    zt, info = spla.minres(Hop, D*b, x0=zt0, rtol=rtol, maxiter=maxiter)
    t1 = time.time()
    z = D*zt
    u = z[:n].copy(); lam = z[n:].copy()
    print(f"  {tag} MINRES done info={info} ({time.time()-t0:.0f}s)", flush=True)

    for k in range(nref):
        r = Kff@u - g
        stat = r - CW.T@lam
        slackW = CW@u + cD[W]
        rel = float(np.linalg.norm(stat)/max(np.linalg.norm(r), 1e-30))
        print(f"  {tag} refine {k}: stat_rel={rel:.3e} min slackW={slackW.min():+.3e}",
              flush=True)
        if rel < 5e-10:
            break
        # rho = b - Hx, with H[u;lam] = [K u - C'lam; -C u] and b = [g; cD]:
        #   rho_u = -(K u - g - C'lam) = -stat
        #   rho_l = cD + C u          = +slackW
        # The second block must carry +slackW; -slackW amplifies the
        # constraint residual every round instead of removing it.
        res = np.concatenate([-stat, slackW])
        # corrections need real convergence: the initial tolerance is measured
        # against ||b||, which is huge, so it reports success far too early.
        zc_t, cinfo = spla.minres(Hop, D*res, rtol=1e-8, maxiter=8000)
        print(f"  {tag}   corr {k} done info={cinfo} ({time.time()-t1:.0f}s)", flush=True)
        t1 = time.time()
        zc = D*zc_t
        u += zc[:n]; lam += zc[n:]
    return u, lam, int(info)


def gate(d, W, u, lam):
    Kff, C_F, cD, g = d["Kff"], d["C_F"], d["cD"], d["g"]
    CW = C_F[W]
    r = Kff@u - g
    stat = r - CW.T@lam
    rel = float(np.linalg.norm(stat)/max(np.linalg.norm(r), 1e-30))
    slack = C_F@u + cD
    lam_full = np.zeros(len(cD)); lam_full[W] = lam
    comp = float(np.max(np.abs(lam_full*slack)))
    u_full = d["val"].copy(); u_full[d["idx"]] = u
    overshoot = float(max(0.0, -u_full.min()))
    si_set = set(np.unique(d["tris"][d["si"]]).tolist())
    ox = [int(i) for i in np.unique(d["tris"][d["mat"] == M.OX]) if int(i) not in si_set]
    ox_free = np.intersect1d(np.array(ox, dtype=int), d["idx"])
    Ku = d["K"] @ u_full + d["f"]
    ox_res = float(np.abs(Ku[ox_free]).max()) if len(ox_free) else 0.0
    viol = int((slack < -1e-9).sum())
    return dict(rel=rel, slack=slack, lam_full=lam_full, comp=comp,
                ox=ox_res, viol=viol, overshoot=overshoot, u_full=u_full,
                lam_min=float(lam_full.min()))


def ok(g):
    return (g["rel"] < 1e-8 and g["comp"] < 1e-10 and
            g["slack"].min() >= -1e-9 and g["lam_min"] >= -1e-9 and
            g["ox"] < 1e-9 and g["overshoot"] < 1e-9)


def main():
    d = build()
    chk = np.load("p2_B1X_osqp_checkpoint.npz")
    u = chk["x"].copy()
    sK = float(chk["sK"])
    lam_osqp = -sK*chk["y"]
    print(f"[B1X-B] checkpoint status={chk['status']} iters={int(chk['iters'])} "
          f"free={d['nfree']} sK={sK:.3e} lam_max={lam_osqp.max():.3e}", flush=True)

    C_F, cD = d["C_F"], d["cD"]
    W = np.flatnonzero(lam_osqp > 1e-8)
    lam = lam_osqp[W].copy()
    print(f"[B1X-B] seed |W|={len(W)} (lam>1e-8)", flush=True)

    best = None
    for outer in range(6):
        u, lam, info = polish(d, W, u, lam, tag=f"[{outer}]")
        g = gate(d, W, u, lam)
        print(f"  outer {outer}: |W|={len(W)} stat={g['rel']:.3e} "
              f"min slack={g['slack'].min():+.3e} lam_min={g['lam_min']:+.3e} "
              f"viol={g['viol']}", flush=True)
        key = (g["rel"], -g["slack"].min())
        if best is None or key < best[0]:
            best = (key, g, W.copy(), u.copy(), lam.copy())
        if ok(g):
            print("  GATE: PASS", flush=True)
            break
        add = np.flatnonzero(g["slack"] < -1e-12)
        drop = W[g["lam_full"][W] < -1e-9]
        if len(add) == 0 and len(drop) == 0:
            print("  active set stable", flush=True)
            break
        W = np.setdiff1d(np.union1d(W, add), drop)
        lam = np.zeros(len(W))
        print(f"  update: +{len(add)} -{len(drop)} -> |W|={len(W)}", flush=True)

    _, g, W, u, lam = best
    W = np.asarray(W, dtype=int)
    passed = ok(g)
    print(f"[B1X] final |W|={len(W)} stat={g['rel']:.3e} comp={g['comp']:.3e} "
          f"min slack={g['slack'].min():+.3e} viol={g['viol']} "
          f"lam_min={g['lam_min']:+.3e} oxide={g['ox']:.3e} "
          f"overshoot={g['overshoot']:+.3e}", flush=True)
    print(f"  GATE: {'PASS' if passed else 'FAIL'}", flush=True)

    np.savez_compressed("p2_B1X_scaled.npz", nodes=d["nodes"], tris=d["tris"],
                        mat=d["mat"], u=g["u_full"], lam=g["lam_full"], active=W)
    json.dump(dict(active=int(len(W)), stationarity=g["rel"],
                   complementarity=g["comp"],
                   min_slack=float(g["slack"].min()),
                   violations=int(g["viol"]),
                   lam_min=float(g["lam_min"]), oxide=g["ox"],
                   overshoot=float(g["overshoot"]), gate=bool(passed),
                   minres_info=int(info)), open("p2_B1X_scaled.json", "w"), indent=2)
    print("wrote p2_B1X_scaled.npz / .json")


if __name__ == "__main__":
    main()
