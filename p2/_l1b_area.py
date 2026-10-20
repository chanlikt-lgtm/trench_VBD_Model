import json
A = {"r3": (0.020926837673, 0.020930251232),
     "r4": (0.021574089791, 0.021578224392),
     "mid": (0.021789286073, 0.021793581617),
     "L1": (0.022071011922, 0.022075109085)}
p = json.load(open("p2_L1b_post.json"))
A["L1b"] = (p["area"]["lower"], p["area"]["upper"])
for k in ("r3","r4","mid","L1","L1b"):
    print(f"{k:4s} [{A[k][0]:.12f}, {A[k][1]:.12f}]  width {A[k][1]-A[k][0]:.3e}")
def d(a,b):
    return (A[b][0]-A[a][1], A[b][1]-A[a][0])
for a,b in (("mid","L1b"),("L1b","L1"),("mid","L1"),("r4","mid")):
    lo,up = d(a,b)
    print(f"{a}->{b}: [{lo:.12e}, {up:.12e}]  mid {(lo+up)/2:.12e}")
r1 = d("L1b","L1"); r0 = d("mid","L1b")
print(f"ratio L1b->L1 / mid->L1b: [{r1[0]/r0[1]:.3f}, {r1[1]/r0[0]:.3f}]")
