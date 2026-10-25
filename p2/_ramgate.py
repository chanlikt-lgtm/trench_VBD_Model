import ctypes, ctypes.wintypes as wt
class MS(ctypes.Structure):
    _fields_=[("dwLength",wt.DWORD),("dwMemoryLoad",wt.DWORD),("ullTotalPhys",ctypes.c_ulonglong),
              ("ullAvailPhys",ctypes.c_ulonglong),("ullTotalPageFile",ctypes.c_ulonglong),
              ("ullAvailPageFile",ctypes.c_ulonglong),("ullTotalVirtual",ctypes.c_ulonglong),
              ("ullAvailVirtual",ctypes.c_ulonglong),("ullAvailExtendedVirtual",ctypes.c_ulonglong)]
m=MS(); m.dwLength=ctypes.sizeof(MS); ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(m))
print("CP-L1.0 physical %.0f MB  available %.0f MB  gate 3000 MB -> %s" %
      (m.ullTotalPhys/1024**2, m.ullAvailPhys/1024**2, "PASS" if m.ullAvailPhys/1024**2>=3000 else "HOLD"))
