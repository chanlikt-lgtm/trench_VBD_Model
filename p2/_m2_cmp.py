import json, numpy as np
A={"L1b":(0.021842101270,0.021845749070),"M2":(0.021839484478,0.021843096605),
   "L1":(0.022071011922,0.022075109085),"mid":(0.021789286073,0.021793581617)}
def diff(a,b): return (A[b][0]-A[a][1], A[b][1]-A[a][0])
print("area intervals")
for k,v in A.items(): print(f"  {k:4s} [{v[0]:.12f}, {v[1]:.12f}] width {v[1]-v[0]:.3e}")
for a,b in (("L1b","M2"),("M2","L1"),("L1b","L1"),("mid","L1b")):
    lo,up=diff(a,b); print(f"  {a}->{b}: [{lo:+.6e}, {up:+.6e}]  mid {(lo+up)/2:+.6e}")
print("\nfront lower edges at eps=1e-8")
cuts={
 "L1b":{"2.35":(-0.0774,0.0719),"2.41":(-0.2014,0.1946)},
 "M2" :{"2.35":(-0.0731,0.0708),"2.41":(-0.1922,0.1948)},
 "L1" :{"2.35":(-0.0719,0.0680),"2.41":(-0.1996,0.1940)},
}
for x in ("2.35","2.41"):
    print(f"  x={x}")
    for k in ("L1b","M2","L1"): print(f"    {k:4s} {cuts[k][x]}")
    print(f"    |M2-L1b| = {abs(cuts['M2'][x][0]-cuts['L1b'][x][0]):.4f}   "
          f"|L1-L1b| = {abs(cuts['L1'][x][0]-cuts['L1b'][x][0]):.4f}")
