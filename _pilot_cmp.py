import json
d=json.load(open("p2_3_spatial_S0/levels/pilot_summary.json"))
L0,L1=d[0],d[1]
print(f"{'':12s}{'L0':>14s}{'L1':>14s}{'ratio':>10s}")
for k in ("nodes","elements","si_area","exact_charged","centroid_charged","aspect_max"):
    a,b=L0[k],L1[k]
    r=b/a if a else float('nan')
    print(f"{k:12s}{a:14.9g}{b:14.9g}{r:10.4f}")
print("\nregion medians (p50/p90) and ratio")
for k in ("band<=0.04","trans0.04-0.06","tip","bodycut","background"):
    if k in L0["regions"] and k in L1["regions"]:
        a=L0["regions"][k]; b=L1["regions"][k]
        print(f"  {k:16s} L0 {a['p50']:.6f}/{a['p90']:.6f}  L1 {b['p50']:.6f}/{b['p90']:.6f}  "
              f"hL1/hL0 p50={b['p50']/a['p50']:.4f} p90={b['p90']/a['p90']:.4f}  n={b['n']}")
print(f"\nSi area drift  L1-L0 = {L1['si_area']-L0['si_area']:+.6e} "
      f"({100*(L1['si_area']-L0['si_area'])/L0['si_area']:+.2e} %)")
print(f"exact charged drift    = {L1['exact_charged']-L0['exact_charged']:+.6e}")
