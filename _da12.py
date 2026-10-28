import json
L0=json.load(open("p2_3_spatial_S0/A_S0-L0.json"))
L1=json.load(open("p2_3_spatial_S0/A_S0-L1.json"))
L2=json.load(open("p2_3_spatial_S0/A_S0-L2.json"))
l0,u0=L0["lower"],L0["upper"]; l1,u1=L1["lower"],L1["upper"]; l2,u2=L2["lower"],L2["upper"]
d01=(l1-u0, u1-l0); d12=(l2-u1, u2-l1)
sup12=max(abs(d12[0]),abs(d12[1])); inf01=(0 if d01[0]<=0<=d01[1] else min(abs(d01[0]),abs(d01[1])))
print("A_L0 [%.15f, %.15f]" % (l0,u0))
print("A_L1 [%.15f, %.15f]" % (l1,u1))
print("A_L2 [%.15f, %.15f]" % (l2,u2))
print("dA01 [%+.9e, %+.9e]  width %.3e" % (d01[0],d01[1],d01[1]-d01[0]))
print("dA12 [%+.9e, %+.9e]  width %.3e" % (d12[0],d12[1],d12[1]-d12[0]))
print("sup|dA12| %.9e  <  inf|dA01| %.9e  -> %s" % (sup12, inf01, sup12<inf01))
import math
r=2**-0.5
r_lo=d12[0]/d01[1]; r_hi=d12[1]/d01[0]
print("ratio dA12/dA01 [%.4f, %.4f]" % (r_lo,r_hi))
print("candidate p = ln(ratio)/ln(r): [%.3f, %.3f]" % (math.log(r_lo)/math.log(r), math.log(r_hi)/math.log(r)))
rec=dict(L0=[l0,u0],L1=[l1,u1],L2=[l2,u2],dA01=list(d01),dA12=list(d12),
         sup_abs_dA12=sup12, inf_abs_dA01=inf01, shrinkage=bool(sup12<inf01),
         ratio=[r_lo,r_hi], candidate_p=[math.log(r_lo)/math.log(r), math.log(r_hi)/math.log(r)])
json.dump(rec, open("p2_3_spatial_S0/DA_S0_12.json","w"), indent=2)
print("wrote DA_S0_12.json")
