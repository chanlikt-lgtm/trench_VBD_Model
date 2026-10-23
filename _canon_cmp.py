import json, numpy as np, sys, os
sys.path.insert(0, os.path.join(os.getcwd(),"p2"))
import mesh_stage2b as M
d=json.load(open("p2_3_spatial_S0/levels/pilot_summary.json"))
L0,L1=d[0],d[1]
print("CANONICAL PAIR")
for k in ("nodes","elements","si_area","exact_charged","centroid_charged","aspect_max","n_arc","n_arcs"):
    print(f"  {k:18s} L0 {L0[k]!s:>22s}   L1 {L1[k]!s:>22s}")
print(f"\n  Si delta      {L1['si_area']-L0['si_area']:+.3e}")
print(f"  charged delta {L1['exact_charged']-L0['exact_charged']:+.3e}")
print("\nREGIONAL h (p50 / p90) and ratio")
for k in ("band<=0.04","trans0.04-0.06","tip","bodycut","background"):
    a=L0["regions"].get(k); b=L1["regions"].get(k)
    if a and b:
        print(f"  {k:15s} {a['p50']:.6f}/{a['p90']:.6f}  ->  {b['p50']:.6f}/{b['p90']:.6f}"
              f"   r50={b['p50']/a['p50']:.4f} r90={b['p90']/a['p90']:.4f}  n1={b['n']}")
R = M.A_BOT - M.A_ARC_END
nch = 64
sag = R*(1-np.cos((np.pi/2)/(2*nch)))
print(f"\nCANONICAL GEOMETRY ERROR (65 nodes / 64 chords per arc)")
print(f"  arc radius R = {R:.6f} um   sagitta = {sag:.3e} um")
print(f"  Si area canonical {L0['si_area']:.12f} vs curved-arc L0 1.954559924360")
print(f"  delta = {L0['si_area']-1.954559924360:+.6e} um^2 "
      f"({100*(L0['si_area']-1.954559924360)/1.954559924360:+.2e} %)")
