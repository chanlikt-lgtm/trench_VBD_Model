import json, time, math
import numpy as np
import sys, os
sys.path.insert(0, "p2")
import p2_area_bounds as AB

AB.AREA_TOL = 1e-13
z = np.load("../p2_3_spatial_S0/CP-M2-C_L0_polished.npz")
nodes, tris, mat, u = z["nodes"], z["tris"], z["mat"], z["u"]
area = np.abs(__import__("p2_audit").element_geom(nodes, tris)[2])
amin = float(area.min())
depth_bound = math.log(amin/AB.AREA_TOL)/math.log(4)
t0=time.time()
lo, up, nch, nm = AB.bounds(nodes, tris, mat, u)
dt=time.time()-t0
rec=dict(solution="S0-L0", commit="a3a5248", AREA_TOL=AB.AREA_TOL,
         MAX_DEPTH=AB.MAX_DEPTH, min_element_area=amin,
         max_terminating_depth=round(depth_bound,2),
         depth_limit_binding=bool(depth_bound > AB.MAX_DEPTH),
         lower=float(lo), upper=float(up), width=float(up-lo),
         charged=int(nch), mixed=int(nm), seconds=round(dt,1))
print(json.dumps(rec, indent=2))
json.dump(rec, open("../p2_3_spatial_S0/A_S0-L0.json","w"), indent=2)
print("wrote A_S0-L0.json")

