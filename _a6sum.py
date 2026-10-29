import json
d=json.load(open("p2_3_spatial_S0/CP-L2-A6_arc_region.json"))
L=d["levels"]
def iv(t): return (L[t]["arc_lower"], L[t]["arc_upper"])
def cov(t): return (L[t]["comp_lower"], L[t]["comp_upper"])
def sub(a,b): return (a[0]-b[1], a[1]-b[0])
print("arc bounds (1e-13):")
for t in ("L0","L1","L2"):
    print(f"  {t} [{iv(t)[0]:.9e}, {iv(t)[1]:.9e}]")
print("\narc increments:")
for a,b in (("L0","L1"),("L1","L2")):
    x=sub(iv(b),iv(a)); print(f"  {a}->{b} [{x[0]:+.6e}, {x[1]:+.6e}] contains0={x[0]<=0<=x[1]}")
print("\ncomplement bounds:")
for t in ("L0","L1","L2"):
    print(f"  {t} [{cov(t)[0]:.12f}, {cov(t)[1]:.12f}]")
print("\ncomplement increments:")
for a,b in (("L0","L1"),("L1","L2")):
    x=sub(cov(b),cov(a)); print(f"  {a}->{b} [{x[0]:+.9e}, {x[1]:+.9e}]")
w=json.load(open("p2_3_spatial_S0/DA_S0_12.json"))
print("\nwhole-domain dA01 [%+.9e, %+.9e]" % tuple(w["dA01"]))
print("whole-domain dA12 [%+.9e, %+.9e]" % tuple(w["dA12"]))
rec=dict(arc={t:list(iv(t)) for t in ("L0","L1","L2")},
         comp={t:list(cov(t)) for t in ("L0","L1","L2")},
         arc_d01=list(sub(iv("L1"),iv("L0"))), arc_d12=list(sub(iv("L2"),iv("L1"))),
         comp_d01=list(sub(cov("L1"),cov("L0"))), comp_d12=list(sub(cov("L2"),cov("L1"))),
         whole_d01=w["dA01"], whole_d12=w["dA12"], AREA_TOL=1e-13)
json.dump(rec, open("p2_3_spatial_S0/CP-L2-A6_summary.json","w"), indent=2)
print("\nwrote CP-L2-A6_summary.json")
