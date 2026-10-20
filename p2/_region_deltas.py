import json
d=json.load(open("p2_area_regions.json"))
whole={"mid":(0.021789286073,0.021793581617),"L1b":(0.021842101270,0.021845749070),
       "M2":(0.021839484478,0.021843096605),"L1":(0.022071011922,0.022075109085)}
print("conservation: sum_lower <= whole_lower ; sum_upper >= whole_upper")
for m in ("mid","L1b","M2","L1"):
    s=d[m]["SUM"]; w=whole[m]
    print(f"  {m:4s} sum [{s[0]:.12f},{s[1]:.12f}]  whole [{w[0]:.12f},{w[1]:.12f}]  "
          f"dL={s[0]-w[0]:+.2e} dU={s[1]-w[1]:+.2e}  ok={s[0]<=w[0] and s[1]>=w[1]}")
def dv(reg,a,b):
    A=d[a][reg]; B=d[b][reg]
    return (B[0]-A[1], B[1]-A[0])
print("\nchange intervals per region")
for a,b in (("mid","L1b"),("L1b","M2"),("L1b","L1"),("M2","L1")):
    print(f"  {a}->{b}")
    for r in ("A","B","REM","SUM"):
        lo,up=dv(r,a,b); print(f"    {r:4s} [{lo:+.4e}, {up:+.4e}]")
print("\nshare of the L1b->L1 increase (mid of interval)")
tot=sum((dv(r,"L1b","L1")[0]+dv(r,"L1b","L1")[1])/2 for r in ("A","B","REM"))
for r in ("A","B","REM"):
    lo,up=dv(r,"L1b","L1"); m=(lo+up)/2
    base=(d["L1b"][r][0]+d["L1b"][r][1])/2
    print(f"  {r:4s} {m:+.4e}  share {100*m/tot:5.1f}%  relative {100*m/base:+.2f}%")
print("\nlayout redistribution L1b->M2")
for r in ("A","B","REM"):
    lo,up=dv(r,"L1b","M2"); print(f"  {r:4s} [{lo:+.3e}, {up:+.3e}]")
