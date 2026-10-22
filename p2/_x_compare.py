import json
C={"B0":(0.021644935001,0.021649266844),"B1":(0.021987652318,0.021991846424)}
X={"B0X":(0.021644934927,0.021649266770),"B1X":(0.021987627914,0.021991822020)}
def d(a,b): return (b[0]-a[1], b[1]-a[0])
dC=d(C["B0"],C["B1"]); dX=d(X["B0X"],X["B1X"])
print(f"DeltaA_centroid  [{dC[0]:+.9e}, {dC[1]:+.9e}]  width {dC[1]-dC[0]:.2e}")
print(f"DeltaA_exact     [{dX[0]:+.9e}, {dX[1]:+.9e}]  width {dX[1]-dX[0]:.2e}")
corr=(dX[0]-dC[1], dX[1]-dC[0])
print(f"correction dX-dC [{corr[0]:+.3e}, {corr[1]:+.3e}]  width {corr[1]-corr[0]:.2e}  straddles_zero={corr[0]<0<corr[1]}")
print(f"  as fraction of increment: [{100*corr[0]/dC[1]:+.2f}%, {100*corr[1]/dC[0]:+.2f}%]")
for m in ("B0","B1"):
    k=m+"X"; dm=(X[k][0]-C[m][1], X[k][1]-C[m][0])
    print(f"delta_{m}: [{dm[0]:+.3e}, {dm[1]:+.3e}]  width {dm[1]-dm[0]:.2e}  straddles_zero={dm[0]<0<dm[1]}")
print("\nF7:")
for t in ("B0X","B1X"):
    j=json.load(open(f"p2_{t}_scaled.json"))
    print(f"  {t}: stat {j['stationarity']:.3e} comp {j['complementarity']:.3e} "
          f"slack {j['min_slack']:+.2e} viol {j['violations']} oxide {j['oxide']:.2e} "
          f"overshoot {j['overshoot']:+.2e} gate {j['gate']}")
