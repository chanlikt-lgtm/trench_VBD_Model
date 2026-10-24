p="p2_cpm2_assembly.py"; s=open(p,encoding="utf-8").read()
s=s.replace('''def peak_rss_mb():
    pmc = PMC(); pmc.cb = ctypes.sizeof(PMC)
    ctypes.windll.psapi.GetProcessMemoryInfo(
        ctypes.windll.kernel32.GetCurrentProcess(), ctypes.byref(pmc), pmc.cb)
    return pmc.PeakWorkingSetSize / (1024**2)''',
'''def peak_rss_mb():
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
    return pmc.PeakWorkingSetSize / (1024**2)''')
open(p,"w",encoding="utf-8").write(s)
print("patched RSS probe")
