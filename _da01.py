import json
L0=json.load(open("p2_3_spatial_S0/A_S0-L0.json"))
L1=json.load(open("p2_3_spatial_S0/A_S0-L1.json"))
l0,u0=L0["lower"],L0["upper"]
l1,u1=L1["lower"],L1["upper"]
lo,up=l1-u0, u1-l0
print("A_S0-L0 [%.15f, %.15f] width %.3e" % (l0,u0,u0-l0))
print("A_S0-L1 [%.15f, %.15f] width %.3e" % (l1,u1,u1-l1))
print("dA      [%+.9e, %+.9e]  width %.3e" % (lo,up,up-lo))
print("positive:", lo>0, "| contains zero:", lo<=0<=up)
print("relative to L0 midpoint: %.4f %%" % (100*0.5*(lo+up)/(0.5*(l0+u0))))
rec=dict(L0=[l0,u0], L1=[l1,u1], delta=[lo,up], delta_width=up-lo,
         positive=bool(lo>0), contains_zero=bool(lo<=0<=up),
         relative_percent=100*0.5*(lo+up)/(0.5*(l0+u0)),
         L0_commit="76f6f42", L1_commit="e1c8a2d", AREA_TOL=1e-13)
json.dump(rec, open("p2_3_spatial_S0/DA_S0_01.json","w"), indent=2)
print("wrote DA_S0_01.json")
