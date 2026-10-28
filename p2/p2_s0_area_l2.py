import json, time, math, hashlib
import numpy as np
import p2_audit as A
import p2_area_bounds as AB

MESH = "../p2_3_spatial_S0/levels/S0_L2_p2mesh.npz"
SOL  = "../p2_3_spatial_S0/CP-L2_L2_polished.npz"
MESH_SHA = "7E6CF12400B0A1308DC23BC7DEC89DA96EF30766784C90D0CB433166523C9136"

def sha(p): return hashlib.sha256(open(p,"rb").read()).hexdigest().upper()

prov = dict(mesh_sha=sha(MESH), mesh_expected=MESH_SHA,
            mesh_match=(sha(MESH)==MESH_SHA), solution_sha=sha(SOL))
f7 = json.load(open("../p2_3_spatial_S0/CP-L2_L2.json"))
prov.update(f7_gate=f7["f7"]["gate"], f7_stationarity=f7["f7"]["stationarity"],
            free_dofs=f7["dofs"], compacted_dofs=913595)
print("CP-L2-A0", json.dumps(prov, indent=2), flush=True)
assert prov["mesh_match"] and prov["f7_gate"] and f7["dofs"]==638124

AB.AREA_TOL = 1e-13
z = np.load(SOL)
nodes, tris, mat, u = z["nodes"], z["tris"], z["mat"], z["u"]
area = np.abs(A.element_geom(nodes, tris)[2]); amin=float(area.min())
db = math.log(amin/AB.AREA_TOL)/math.log(4)
t0=time.time(); lo,up,nch,nm = AB.bounds(nodes,tris,mat,u); dt=time.time()-t0
rec = dict(solution="S0-L2", commit="6e8d9e1", AREA_TOL=AB.AREA_TOL, MAX_DEPTH=AB.MAX_DEPTH,
           min_element_area=amin, max_terminating_depth=round(db,2),
           depth_limit_binding=bool(db>AB.MAX_DEPTH),
           lower=float(lo), upper=float(up), width=float(up-lo),
           charged=int(nch), mixed=int(nm), seconds=round(dt,1),
           free_dofs=f7["dofs"], compacted_dofs=913595, provenance=prov)
print(json.dumps(rec, indent=2), flush=True)
json.dump(rec, open("../p2_3_spatial_S0/A_S0-L2.json","w"), indent=2)
print("wrote A_S0-L2.json")
