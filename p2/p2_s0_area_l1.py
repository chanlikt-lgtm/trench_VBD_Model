import json, time, math, hashlib
import numpy as np
import p2_audit as A
import p2_area_bounds as AB

MESH = "../p2_3_spatial_S0/levels/S0_L1_p2mesh.npz"
SOL  = "../p2_3_spatial_S0/CP-L1_L1_polished.npz"
MESH_SHA = "7FA1C69DA097612960A2668335A014894A0DA4D9FF5C96D623E10FCAE0A3D1BB"

def sha(p):
    return hashlib.sha256(open(p,"rb").read()).hexdigest().upper()

prov = dict(mesh_sha=sha(MESH), mesh_expected=MESH_SHA,
            mesh_match=(sha(MESH) == MESH_SHA), solution_sha=sha(SOL))
f7 = json.load(open("../p2_3_spatial_S0/CP-L1_L1.json"))
prov["f7_gate"] = f7["f7"]["gate"]
prov["f7_stationarity"] = f7["f7"]["stationarity"]
prov["free_dofs"] = f7["dofs"]
prov["compacted_dofs"] = 459288
print("CP-M4.0", json.dumps(prov, indent=2), flush=True)
assert prov["mesh_match"] and prov["f7_gate"]

AB.AREA_TOL = 1e-13
z = np.load(SOL)
nodes, tris, mat, u = z["nodes"], z["tris"], z["mat"], z["u"]
area = np.abs(A.element_geom(nodes, tris)[2])
amin = float(area.min())
db = math.log(amin/AB.AREA_TOL)/math.log(4)
t0 = time.time()
lo, up, nch, nm = AB.bounds(nodes, tris, mat, u)
dt = time.time()-t0
rec = dict(solution="S0-L1", commit="e1c8a2d", AREA_TOL=AB.AREA_TOL, MAX_DEPTH=AB.MAX_DEPTH,
           min_element_area=amin, max_terminating_depth=round(db,2),
           depth_limit_binding=bool(db > AB.MAX_DEPTH),
           lower=float(lo), upper=float(up), width=float(up-lo),
           charged=int(nch), mixed=int(nm), seconds=round(dt,1),
           free_dofs=f7["dofs"], compacted_dofs=459288, provenance=prov)
print(json.dumps(rec, indent=2), flush=True)
json.dump(rec, open("../p2_3_spatial_S0/A_S0-L1.json","w"), indent=2)
print("wrote A_S0-L1.json")
