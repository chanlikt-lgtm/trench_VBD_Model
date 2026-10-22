import json
A={"B0":(0.021644935001,0.021649266844),"B1":(0.021987652318,0.021991846424),
   "L1b":(0.021842101270,0.021845749070),"L1":(0.022071011922,0.022075109085)}
def d(a,b): return (A[b][0]-A[a][1], A[b][1]-A[a][0])
for a,b in (("B0","B1"),("L1b","B1"),("L1b","L1")):
    lo,up=d(a,b); print(f"{a}->{b}: [{lo:+.6e}, {up:+.6e}]")
    print(f"   relative to {a}: [{100*lo/A[a][1]:+.3f}%, {100*up/A[a][0]:+.3f}%]")
p0=json.load(open("p2_B0_post.json")); p1=json.load(open("p2_B1_post.json"))
print("\ncuts at 1e-8")
for x in ("2.35","2.41"):
    print(f"  x={x}  B0 {[(round(a,4),round(b,4)) for a,b in p0['cuts'][x]['1e-08']['main']]}")
    print(f"         B1 {[(round(a,4),round(b,4)) for a,b in p1['cuts'][x]['1e-08']['main']]}")
print("\ndonor: dSum(f)=0.3035 (1.54e-4 rel);  area: 1.58% rel  -> ratio ~103")
