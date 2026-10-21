import numpy as np, json
import p2_audit as A
from p2_bernstein import B
M=A.M; EPS=1e-8
REG={"A":(2.32,2.37,-0.13,0.13),"B":(2.395,2.420,-0.24,0.24)}
z=np.load("p2_B1_p2mesh.npz"); nodes,tris,mat,u=None,None,z["mat"],None
nodes,tris,mat=z["nodes"],z["tris"],z["mat"]
# boundary needs a solution; use the frozen L1b field only for coverage of the band geometry,
# so measure size of elements near the frozen contour instead:
band=np.load("p2_band_contour.npz")["pts"]
area=np.abs(A.element_geom(nodes,tris)[2]); P=nodes[tris][:,:3]
cx=P[:,:,0].mean(1); cy=P[:,:,1].mean(1)
si=mat==M.SI
h=np.sqrt(4*area/np.sqrt(3))
# distance from element centroids to band polyline (coarse: min over band points)
from scipy.spatial import cKDTree
tree=cKDTree(band)
d,_=tree.query(np.column_stack([cx,cy]))
near=(d<=0.045)&si
print(f"elements within 0.045 um of frozen band: {int(near.sum())}")
print(f"  h percentiles: p10={np.percentile(h[near],10):.5f} med={np.median(h[near]):.5f} "
      f"p90={np.percentile(h[near],90):.5f} max={h[near].max():.5f}")
print(f"  fraction with h<=0.0045: {100*np.mean(h[near]<=0.0045):.1f}%")
out={}
for k,r in REG.items():
    m=near&(cx>=r[0])&(cx<=r[1])&(cy>=r[2])&(cy<=r[3])
    if m.sum(): out[k]=dict(n=int(m.sum()),med=float(np.median(h[m])),p90=float(np.percentile(h[m],90)))
st=near&(cx>2.37)&(cx<2.395)
if st.sum(): out["strip"]=dict(n=int(st.sum()),med=float(np.median(h[st])),p90=float(np.percentile(h[st],90)))
for k,v in out.items(): print(f"  {k:5s} n={v['n']:5d} med={v['med']:.5f} p90={v['p90']:.5f}")
json.dump(dict(near=int(near.sum()),band_h=dict(p10=float(np.percentile(h[near],10)),
      med=float(np.median(h[near])),p90=float(np.percentile(h[near],90)),
      frac_le_00045=float(np.mean(h[near]<=0.0045))),regions=out),
      open("p2_B1_coverage.json","w"),indent=2)
print("wrote p2_B1_coverage.json")
