import numpy as np, json
import p2_audit as A
from scipy.spatial import cKDTree
M=A.M
band=np.load("p2_band_contour.npz")["pts"]
tree=cKDTree(band)
def stats(fn,tag):
    z=np.load(fn); nodes,tris,mat=z["nodes"],z["tris"],z["mat"]
    area=np.abs(A.element_geom(nodes,tris)[2]); P=nodes[tris][:,:3]
    cx=P[:,:,0].mean(1); cy=P[:,:,1].mean(1)
    si=mat==M.SI
    h=np.sqrt(4*area/np.sqrt(3))
    d,_=tree.query(np.column_stack([cx,cy]))
    groups={"near<=0.045":si&(d<=0.045),"mid 0.045-0.15":si&(d>0.045)&(d<=0.15),
            "far>0.15":si&(d>0.15)}
    out={"tag":tag,"si_n":int(si.sum()),"si_area":float(area[si].sum())}
    for k,m in groups.items():
        if m.sum()<5: continue
        out[k]=dict(n=int(m.sum()),area=float(area[m].sum()),
                    med=float(np.median(h[m])),p10=float(np.percentile(h[m],10)),
                    p90=float(np.percentile(h[m],90)))
    xs=[]
    for a in (0.0,1.0,2.0,2.20,2.30,2.32,2.42,2.50):
        m=si&(cx>=a)&(cx<a+0.10)
        if m.sum()>20: xs.append((round(a,2),int(m.sum()),float(np.median(h[m]))))
    out["x_bands"]=xs
    return out
res=[stats("p2_L1b_scaled.npz","L1b"),stats("p2_B1_p2mesh.npz","B1")]
for r in res:
    print(f"\n=== {r['tag']}  Si elements {r['si_n']}  Si area {r['si_area']:.9f}")
    for k in ("near<=0.045","mid 0.045-0.15","far>0.15"):
        if k in r:
            v=r[k]; print(f"  {k:16s} n={v['n']:7d} area={v['area']:.6f} "
                          f"h p10/med/p90={v['p10']:.5f}/{v['med']:.5f}/{v['p90']:.5f}")
    print("  x band  n      h_med")
    for a,n,mh in r["x_bands"]:
        print(f"   {a:4.2f}  {n:6d}  {mh:.5f}")
json.dump(res,open("p2_b1_comparability.json","w"),indent=2)
print("\nwrote p2_b1_comparability.json")
