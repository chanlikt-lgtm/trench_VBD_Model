#!/usr/bin/env python3
"""CP-L2: S0-L2 production solve under the frozen manifest.

Chain: prelaunch RAM gate -> frozen-input verification -> OSQP ->
checkpoint save/reload verify -> scaled MINRES polish -> full-set F7.

Authorized scope is L0 only. No area postprocessing, no L1.
"""
from __future__ import annotations
import ctypes
import ctypes.wintypes as wt
import hashlib
import json
import sys
import time

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
import osqp

import p2_audit as A
from p2_bernstein import B
from p2_nested_consistency import assemble_p2
from p2_exact_donor import exact_load

M = A.M
TAG = "S0_L2"
MESH = "../p2_3_spatial_S0/levels/S0_L2_p2mesh.npz"
MESH_SHA = "7E6CF12400B0A1308DC23BC7DEC89DA96EF30766784C90D0CB433166523C9136"
EXPECT_DOFS = 913595
EXPECT_ROWS = 1413522
RAM_GATE_MB = 3000.0
RAM_STOP_MB = 1500.0
T0 = time.time()


class MEMSTATUS(ctypes.Structure):
    _fields_ = [("dwLength", wt.DWORD), ("dwMemoryLoad", wt.DWORD),
                ("ullTotalPhys", ctypes.c_ulonglong),
                ("ullAvailPhys", ctypes.c_ulonglong),
                ("ullTotalPageFile", ctypes.c_ulonglong),
                ("ullAvailPageFile", ctypes.c_ulonglong),
                ("ullTotalVirtual", ctypes.c_ulonglong),
                ("ullAvailVirtual", ctypes.c_ulonglong),
                ("ullAvailExtendedVirtual", ctypes.c_ulonglong)]


def ram_mb():
    ms = MEMSTATUS(); ms.dwLength = ctypes.sizeof(MEMSTATUS)
    ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(ms))
    return ms.ullTotalPhys/1024**2, ms.ullAvailPhys/1024**2


def rss_mb():
    import psutil
    return psutil.Process().memory_info().rss/1024**2


def peak_mb():
    import psutil
    return psutil.Process().memory_full_info().peak_wset/1024**2


def mark(stage, marks):
    tot, avail = ram_mb()
    marks[stage] = dict(t=round(time.time()-T0, 1), rss=round(rss_mb()),
                        avail=round(avail), total=round(tot))
    print(f"[{stage}] t={marks[stage]['t']}s rss={marks[stage]['rss']}MB "
          f"avail={marks[stage]['avail']}MB", flush=True)
    if avail < RAM_STOP_MB:
        print(f"STOP: available RAM {avail:.0f} MB < {RAM_STOP_MB} MB", flush=True)
        raise SystemExit(3)


def build():
    z = np.load(MESH)
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
    K, _ = assemble_p2(nodes, tris, mat, np.zeros(len(tris), bool))
    f, _ = exact_load(nodes, tris, mat)
    fixed = np.zeros(len(nodes), bool); val = np.zeros(len(nodes))
    si_nodes = np.unique(tris[si])
    gate_n = np.unique(tris[mat == M.GATE]); shield = np.unique(tris[mat == M.SHIELD])
    src = si_nodes[x[si_nodes] <= 1e-12]
    drn = np.flatnonzero(x >= M.BOTTOM_DEPTH - 1e-9)
    for s in (np.unique(tris[body]), gate_n, shield, src):
        fixed[s] = True
    val[fixed] = A.V_DRAIN
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
                idx=idx, C_F=C_F, cD=cD, K=K, Kff=Kff, g=g, f=f,
                nfree=len(idx), ncons=nC)


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
    return dict(rel=rel, slack=slack, lam_full=lam_full, comp=comp, ox=ox_res,
                viol=viol, overshoot=overshoot, u_full=u_full,
                lam_min=float(lam_full.min()))


def ok(g):
    return (g["rel"] < 1e-8 and g["comp"] < 1e-10 and g["slack"].min() >= -1e-9
            and g["lam_min"] >= -1e-9 and g["ox"] < 1e-9 and g["overshoot"] < 1e-9)


def main():
    marks = {}
    tot, avail = ram_mb()
    print(f"[C0] physical {tot:.0f} MB, available {avail:.0f} MB "
          f"(gate {RAM_GATE_MB:.0f} MB)", flush=True)
    if avail < RAM_GATE_MB:
        print("HOLD: insufficient available RAM at launch")
        json.dump(dict(status="HOLD", reason="prelaunch RAM below gate",
                       available_mb=avail, gate_mb=RAM_GATE_MB),
                  open("../p2_3_spatial_S0/CP-L2_L2.json", "w"), indent=2)
        return 2
    mark("start", marks)

    sha = hashlib.sha256(open(MESH, "rb").read()).hexdigest().upper()
    print(f"[C1] mesh sha256 {sha}")
    print(f"     match manifest: {sha == MESH_SHA}")
    if sha != MESH_SHA:
        print("HOLD: mesh hash mismatch")
        return 2

    d = build()
    mark("after_assembly", marks)
    print(f"[C1] dofs={d['nfree']} rows={d['ncons']} "
          f"dofs_ok={d['nfree']+int(d['fixed'].sum()) == EXPECT_DOFS} "
          f"rows_ok={d['ncons'] == EXPECT_ROWS}", flush=True)

    Kff, C_F, cD, g = d["Kff"], d["C_F"], d["cD"], d["g"]
    sK = float(Kff.diagonal().max())
    prob = osqp.OSQP()
    prob.setup(P=(Kff/sK).tocsc(), q=(-g/sK), A=C_F,
               l=-cD, u=np.full(len(cD), np.inf),
               verbose=False, eps_abs=1e-8, eps_rel=1e-8,
               max_iter=30000, polishing=False, scaling=10)
    res = prob.solve()
    mark("after_osqp", marks)
    pri = float(getattr(res.info, "prim_res", getattr(res.info, "pri_res", np.nan)))
    dua = float(getattr(res.info, "dual_res", getattr(res.info, "dua_res", np.nan)))
    print(f"[C2] status={res.info.status} iters={res.info.iter} "
          f"obj={res.info.obj_val:.6e} pri={pri:.3e} dua={dua:.3e}", flush=True)

    chk_path = "../p2_3_spatial_S0/CP-L2_L2_checkpoint.npz"
    np.savez_compressed(chk_path, x=res.x, y=res.y, idx=d["idx"], sK=sK,
                        status=str(res.info.status), iters=int(res.info.iter),
                        obj=float(res.info.obj_val), mesh_sha=sha)
    z = np.load(chk_path)
    same = (z["x"].shape == res.x.shape and z["y"].shape == res.y.shape
            and np.array_equal(z["x"], res.x) and np.array_equal(z["y"], res.y))
    print(f"[C3] checkpoint reload: x{z['x'].shape} y{z['y'].shape} identical={same}",
          flush=True)

    u = res.x.copy()
    lam_osqp = -sK*res.y
    W = np.flatnonzero(lam_osqp > 1e-8)
    lam = lam_osqp[W].copy()
    print(f"[C4] seed |W|={len(W)}", flush=True)
    repairs = []
    minres_info = []
    best = None
    for outer in range(6):
        CW = C_F[W]
        n = d["nfree"]
        Hop, D = make_op(Kff, CW, n)
        b = np.concatenate([g, cD[W]])
        z0 = np.concatenate([u/D[:n], lam/D[n:]])
        zt, info = spla.minres(Hop, D*b, x0=z0, rtol=1e-8, maxiter=8000)
        minres_info.append(int(info))
        zz = D*zt; u = zz[:n].copy(); lam = zz[n:].copy()
        for k in range(6):
            r = Kff@u - g
            stat = r - CW.T@lam
            slackW = CW@u + cD[W]
            rel = float(np.linalg.norm(stat)/max(np.linalg.norm(r), 1e-30))
            print(f"  [{outer}] refine {k}: stat_rel={rel:.3e}", flush=True)
            if rel < 5e-10:
                break
            res_c = np.concatenate([-stat, slackW])   # frozen corrected sign
            zc, _ = spla.minres(Hop, D*res_c, rtol=1e-8, maxiter=8000)
            zzc = D*zc
            u += zzc[:n]; lam += zzc[n:]
        gres = gate(d, W, u, lam)
        print(f"  [{outer}] |W|={len(W)} stat={gres['rel']:.3e} "
              f"slack={gres['slack'].min():+.2e} lam_min={gres['lam_min']:+.2e}",
              flush=True)
        if best is None or gres["rel"] < best[0]["rel"]:
            best = (gres, W.copy(), u.copy(), lam.copy())
        if ok(gres):
            break
        add = np.flatnonzero(gres["slack"] < -1e-12)
        drop = W[gres["lam_full"][W] < -1e-9]
        if len(add) == 0 and len(drop) == 0:
            break
        repairs.append(dict(added=int(len(add)), dropped=int(len(drop))))
        W = np.setdiff1d(np.union1d(W, add), drop)
        lam = np.zeros(len(W))
    mark("after_polish", marks)

    gres, W, u, lam = best
    passed = ok(gres)
    print(f"[C5] stat={gres['rel']:.3e} comp={gres['comp']:.3e} "
          f"slack={gres['slack'].min():+.3e} viol={gres['viol']} "
          f"lam_min={gres['lam_min']:+.3e} oxide={gres['ox']:.3e} "
          f"overshoot={gres['overshoot']:+.3e}")
    print(f"     F7: {'PASS' if passed else 'FAIL'}")

    out = dict(status="PASS" if passed else "FAIL", mesh_sha=sha,
               dofs=d["nfree"], rows=d["ncons"], fixed=int(d["fixed"].sum()),
               osqp=dict(status=str(res.info.status), iters=int(res.info.iter),
                         obj=float(res.info.obj_val), pri=pri, dua=dua),
               checkpoint=dict(path=chk_path, identical=bool(same),
                               x_shape=list(res.x.shape), y_shape=list(res.y.shape)),
               active=int(len(W)), repairs=repairs, minres_info=minres_info,
               f7=dict(stationarity=gres["rel"], complementarity=gres["comp"],
                       min_slack=float(gres["slack"].min()),
                       violations=int(gres["viol"]), lam_min=gres["lam_min"],
                       oxide=gres["ox"], overshoot=gres["overshoot"],
                       gate=bool(passed)),
               marks=marks, peak_rss_mb=round(peak_mb(), 1))
    json.dump(out, open("../p2_3_spatial_S0/CP-L2_L2.json", "w"), indent=2)
    np.savez_compressed("../p2_3_spatial_S0/CP-L2_L2_polished.npz",
                        u=gres["u_full"], lam=gres["lam_full"], active=W,
                        nodes=d["nodes"], tris=d["tris"], mat=d["mat"])
    print("wrote CP-L2_L2.json / _polished.npz")
    return 0


if __name__ == "__main__":
    sys.exit(main())
