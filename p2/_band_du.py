import numpy as np, json
import p2_front_conditioning as F
from p2_exact_metrics import cut_intervals
from p2_topology_clean import clean

MESH=[("mid","p2_p3_rmid_polished.npz"),("L1b","p2_L1b_scaled.npz"),("L1","p2_L1_scaled.npz")]
out={}
for x0 in (2.35,2.41):
    prof={}; fronts={}
    for n,f in MESH:
        z=np.load(f); o,_=F.detect_order(z["nodes"],z["tris"])
        ys,us,ds=F.profile(z["nodes"],z["tris"],z["u"],x0,o)
        prof[n]=(ys,us,ds)
        iv=cut_intervals(z["nodes"],z["tris"],z["mat"],z["u"],x0,1e-8,True)
        m,_s=clean(iv); fronts[n]=m
    grid=np.linspace(-0.25,0.25,2001)
    interp={n:np.interp(grid,prof[n][0],prof[n][1]) for n,_ in MESH}
    ends=[y for n,_ in MESH for iv in fronts[n] for y in iv]
    rec={}
    for a,b in (("L1b","L1"),("mid","L1b"),("mid","L1")):
        d=np.abs(interp[a]-interp[b])
        band=np.zeros_like(grid,dtype=bool)
        for y in ends: band |= np.abs(grid-y)<0.02
        rec[f"{a}-{b}"]=dict(max_all=float(d.max()),
                             max_band=float(d[band].max()) if band.any() else float("nan"),
                             y_band=float(grid[band][np.argmax(d[band])]) if band.any() else float("nan"),
                             p95_band=float(np.percentile(d[band],95)) if band.any() else float("nan"))
    out[f"{x0}"]=rec
    print(f"--- x={x0} ---")
    for k,v in rec.items():
        print(f"  {k:9s} max|du| all={v['max_all']:.3e}  in +-0.02 of front: {v['max_band']:.3e} "
              f"(p95 {v['p95_band']:.3e}) at y={v['y_band']:+.4f}")
json.dump(out,open("p2_front_du_band.json","w"),indent=2)
print("wrote p2_front_du_band.json")
