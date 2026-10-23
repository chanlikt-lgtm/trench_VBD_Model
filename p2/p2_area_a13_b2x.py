import json, time
import numpy as np
import p2_area_bounds as AB
AB.AREA_TOL = 1e-13
z = np.load("p2_B2X_scaled.npz")
t0=time.time()
lo,up,nch,nm = AB.bounds(z["nodes"], z["tris"], z["mat"], z["u"])
dt=time.time()-t0
print(f"[B2X] lower={lo:.15f} upper={up:.15f} width={up-lo:.3e} mixed={nm} ({dt:.0f}s)")
json.dump(dict(lower=float(lo),upper=float(up),width=float(up-lo),mixed=int(nm),
               seconds=round(dt,1)), open("p2_area_a13_b2x.json","w"), indent=2)
