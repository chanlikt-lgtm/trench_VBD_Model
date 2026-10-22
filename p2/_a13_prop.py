import json
d=json.load(open("p2_area_a13.json"))
C={"B0":(d["B0"]["lower"],d["B0"]["upper"]),"B1":(d["B1"]["lower"],d["B1"]["upper"])}
X={"B0X":(d["B0X"]["lower"],d["B0X"]["upper"]),"B1X":(d["B1X"]["lower"],d["B1X"]["upper"])}
def sub(a,b): return (a[0]-b[1], a[1]-b[0])
dC=sub(C["B1"],C["B0"]); dX=sub(X["B1X"],X["B0X"]); dsrc=sub(dX,dC)
print("eight endpoints:")
for k in ("B0","B1","B0X","B1X"):
    print(f"  {k:4s} L={d[k]['lower']:.15f}  U={d[k]['upper']:.15f}  w={d[k]['width']:.3e}  t={d[k]['seconds']}s")
print(f"\nDeltaA_C  [{dC[0]:+.9e}, {dC[1]:+.9e}]  width {dC[1]-dC[0]:.3e}")
print(f"DeltaA_X  [{dX[0]:+.9e}, {dX[1]:+.9e}]  width {dX[1]-dX[0]:.3e}")
print(f"Delta_src [{dsrc[0]:+.9e}, {dsrc[1]:+.9e}]  width {dsrc[1]-dsrc[0]:.3e}")
print(f"contains zero: {dsrc[0] <= 0 <= dsrc[1]}")
if dsrc[0]>0 or dsrc[1]<0:
    print(f"  excludes zero; correction in [{100*dsrc[0]/dC[1]:+.3f}%, {100*dsrc[1]/dC[0]:+.3f}%] of increment")
