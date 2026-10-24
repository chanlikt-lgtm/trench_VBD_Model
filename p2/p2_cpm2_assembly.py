#!/usr/bin/env python3
"""CP-M2-B: assembly-only measurement for one S0 level.

Builds the stiffness operator, the exact clipped donor load and the Bernstein
constraint matrix, records sparsity and finiteness, and reports peak RSS.
Deliberately does NOT import or call osqp, and performs no factorization,
polishing, F7 or area work.

Usage:  python p2_cpm2_assembly.py L0|L1
"""
from __future__ import annotations
import ctypes
import ctypes.wintypes as wt
import json
import os
import sys
import time

import numpy as np
import scipy.sparse as sp

import p2_audit as A
from p2_bernstein import B
from p2_nested_consistency import assemble_p2
from p2_exact_donor import exact_load

M = A.M


class PMC(ctypes.Structure):
    _fields_ = [("cb", wt.DWORD), ("PageFaultCount", wt.DWORD),
                ("PeakWorkingSetSize", ctypes.c_size_t),
                ("WorkingSetSize", ctypes.c_size_t),
                ("QuotaPeakPagedPoolUsage", ctypes.c_size_t),
                ("QuotaPagedPoolUsage", ctypes.c_size_t),
                ("QuotaPeakNonPagedPoolUsage", ctypes.c_size_t),
                ("QuotaNonPagedPoolUsage", ctypes.c_size_t),
                ("PagefileUsage", ctypes.c_size_t),
                ("PeakPagefileUsage", ctypes.c_size_t)]


def peak_rss_mb():
    try:
        import psutil
        return psutil.Process().memory_full_info().peak_wset / (1024**2)
    except Exception:
        pass
    pmc = PMC(); pmc.cb = ctypes.sizeof(PMC)
    k32 = ctypes.windll.kernel32
    k32.GetCurrentProcess.restype = ctypes.c_void_p
    fn = getattr(ctypes.windll.psapi, "GetProcessMemoryInfo", None)
    if fn is None:
        return float("nan")
    fn.argtypes = [ctypes.c_void_p, ctypes.POINTER(PMC), wt.DWORD]
    fn.restype = wt.BOOL
    if not fn(k32.GetCurrentProcess(), ctypes.byref(pmc), pmc.cb):
        return float("nan")
    return pmc.PeakWorkingSetSize / (1024**2)


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


def main():
    tag = sys.argv[1]
    total_mb, avail_mb = ram_mb()
    t0 = time.time()

    z = np.load(f"../p2_3_spatial_S0/levels/S0_{tag}_p2mesh.npz")
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

    tK = time.time()
    K, _ = assemble_p2(nodes, tris, mat, np.zeros(len(tris), bool))
    t_asm = time.time() - tK

    tf = time.time()
    f, ftot = exact_load(nodes, tris, mat)
    t_donor = time.time() - tf

    fixed = np.zeros(len(nodes), bool)
    val = np.zeros(len(nodes))
    si_nodes = np.unique(tris[si])
    gate = np.unique(tris[mat == M.GATE]); shield = np.unique(tris[mat == M.SHIELD])
    src = si_nodes[x[si_nodes] <= 1e-12]
    drn = np.flatnonzero(x >= M.BOTTOM_DEPTH - 1e-9)
    for s in (np.unique(tris[body]), gate, shield, src):
        fixed[s] = True
    val[fixed] = A.V_DRAIN
    fixed[drn] = True; val[drn] = 0.0
    idx = np.flatnonzero(~fixed)

    tC = time.time()
    cons = np.flatnonzero(charged); nC = 6*len(cons)
    C = sp.coo_matrix((np.tile(B.ravel(), len(cons)),
                       (np.repeat(np.arange(nC), 6), np.tile(tris[cons], (1, 6)).ravel())),
                      shape=(nC, len(nodes))).tocsr()
    C_F = C[:, idx].tocsc()
    t_bern = time.time() - tC

    Kff = K[idx][:, idx].tocsc()
    fin = dict(K=bool(np.isfinite(K.data).all()), f=bool(np.isfinite(f).all()),
               C=bool(np.isfinite(C.data).all()), Kff=bool(np.isfinite(Kff.data).all()))

    rec = dict(
        level=tag,
        nodes_raw=int(len(o2n)), nodes_compacted=int(len(nodes)),
        elements=int(len(tris)),
        dofs_total=int(len(nodes)), dofs_free=int(len(idx)),
        dofs_fixed=int(fixed.sum()),
        charged_elements=int(charged.sum()),
        bernstein_rows=int(nC),
        K=dict(shape=list(K.shape), nnz=int(K.nnz),
               mb=float((K.data.nbytes + K.indices.nbytes + K.indptr.nbytes)/1024**2)),
        Kff=dict(shape=list(Kff.shape), nnz=int(Kff.nnz),
                 mb=float((Kff.data.nbytes + Kff.indices.nbytes + Kff.indptr.nbytes)/1024**2)),
        C_F=dict(shape=list(C_F.shape), nnz=int(C_F.nnz),
                 mb=float((C_F.data.nbytes + C_F.indices.nbytes + C_F.indptr.nbytes)/1024**2)),
        finite=fin,
        donor_sum=float(ftot),
        seconds=dict(assembly=round(t_asm, 1), donor=round(t_donor, 1),
                     bernstein=round(t_bern, 1), total=round(time.time()-t0, 1)),
        peak_rss_mb=round(peak_rss_mb(), 1),
        ram=dict(total_mb=round(total_mb), avail_at_start_mb=round(avail_mb)),
    )
    print(json.dumps(rec, indent=2))
    json.dump(rec, open(f"../p2_3_spatial_S0/CP-M2-B_{tag}.json", "w"), indent=2)


if __name__ == "__main__":
    main()
