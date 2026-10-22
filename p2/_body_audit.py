import numpy as np, json
import p2_audit as A
from p2_charge_mask_audit import clip_area
M=A.M; XB=A.BODY_BOTTOM
print(f"RHO = {A.RHO:.6e}")
res={}
for tag,f in (("B0","p2_B0_p2mesh.npz"),("B1","p2_B1_p2mesh.npz")):
    z=np.load(f); nodes,tris,mat=z["nodes"],z["tris"],z["mat"]
    P=nodes[tris][:,:3]
    det=((P[:,1,0]-P[:,0,0])*(P[:,2,1]-P[:,0,1])-(P[:,2,0]-P[:,0,0])*(P[:,1,1]-P[:,0,1]))
    area=0.5*np.abs(det)
    cx=P[:,:,0].mean(1)
    si=mat==M.SI
    xmin=P[:,:,0].min(1); xmax=P[:,:,0].max(1)
    strad=si&(xmin<=XB)&(xmax>XB)
    charged=si&(cx>XB)
    above=si&(xmin>XB); below=si&(xmax<=XB)
    exact=area[above].sum()
    for k in np.flatnonzero(strad):
        exact+=clip_area([tuple(q) for q in P[k]],XB)
    over=0.0; under=0.0
    for k in np.flatnonzero(strad):
        ck=clip_area([tuple(q) for q in P[k]],XB)
        if cx[k]>XB: over+=area[k]-ck
        else: under+=ck
    cent=area[charged].sum()
    h=np.sqrt(4*area/np.sqrt(3))
    nb=np.abs(cx-XB)<=0.05
    res[tag]=dict(straddle_n=int(strad.sum()),straddle_area=float(area[strad].sum()),
                  exact=float(exact),centroid=float(cent),over=float(over),under=float(under),
                  net=float(cent-exact),h_near_med=float(np.median(h[nb])),
                  h_near_n=int(nb.sum()),sumf=float(A.RHO*cent))
    r=res[tag]
    print(f"\n{tag}: straddling {r['straddle_n']} area {r['straddle_area']:.6e}")
    print(f"  exact {r['exact']:.12f}  centroid {r['centroid']:.12f}  net {r['net']:+.3e}")
    print(f"  overcharge {r['over']:+.6e}  undercharge {r['under']:+.6e}  over-under {r['over']-r['under']:+.6e}")
    print(f"  h near body cut (|x-0.54|<=0.05): n={r['h_near_n']} med={r['h_near_med']:.5f}")
    print(f"  sum(f) = {r['sumf']:.9e}")
print(f"\ndelta sum(f) B1-B0 = {res['B1']['sumf']-res['B0']['sumf']:+.6e}")
json.dump(res,open("p2_body_cut_audit.json","w"),indent=2)
print("wrote p2_body_cut_audit.json")
