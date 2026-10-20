import json
f7 = json.load(open("p2_L1_scaled.json"))
post = json.load(open("p2_L1_post.json"))
doc = open("P2_3_LATEST_SESSION_L1_F7_CERTIFIED_L1B_DECISION.md", encoding="utf-8").read()

checks = []
def ck(name, hay, needle, tol=None):
    ok = needle in doc
    checks.append((name, ok, needle))

ck("F7 stationarity", doc, "1.484e-10")
ck("F7 complementarity", doc, "5.778e-14")
ck("F7 min slack", doc, "4.222e-11")
ck("F7 oxide", doc, "4.766e-13")
ck("F7 overshoot", doc, "2.127e-11")
ck("L1 area lower", doc, "0.022071011922")
ck("L1 area upper", doc, "0.022075109085")
ck("L1 area width", doc, "0.000004097163")

a = {k: (v["lower"], v["upper"]) for k, v in post.items()} if False else None

A = {
 "r3":  (0.020926837673, 0.020930251232),
 "r4":  (0.021574089791, 0.021578224392),
 "mid": (0.021789286073, 0.021793581617),
 "L1":  (post["area"]["lower"], post["area"]["upper"]),
}
inc = {}
for a_, b_ in (("r3","r4"), ("r4","mid"), ("mid","L1")):
    inc[b_] = (A[b_][0]-A[a_][1], A[b_][1]-A[a_][0])

print("increments:")
for k, v in inc.items():
    print(f"  {k:>4}: [{v[0]:.12f}, {v[1]:.12f}]")
r1 = (inc["mid"][0]/inc["r4"][1], inc["mid"][1]/inc["r4"][0])
print(f"  ratio r4->mid / r3->r4 : [{r1[0]:.6f}, {r1[1]:.6f}]  (doc 0.324-0.341)")
r2 = (inc["L1"][0]/inc["mid"][1], inc["L1"][1]/inc["mid"][0])
print(f"  factor L1 / mid       : [{r2[0]:.6f}, {r2[1]:.6f}]  (doc 1.264-1.354)")

p1 = {"r3": 0.021196729, "r4": 0.021845237}
for k in ("r3","r4"):
    lo, up = A[k][0]/p1[k], A[k][1]/p1[k]
    print(f"  P2/P1 {k}: [{lo:.6f}, {up:.6f}]  ({100*(lo-1):+.2f}% to {100*(up-1):+.2f}%)")

print("\ndoc contains:")
for n, ok, nd in checks:
    print(f"  {'OK ' if ok else 'MISS'} {n}: {nd}")

print("\ncut x=2.35 eps=1e-8:", [[round(x,4) for x in iv] for iv in post["cuts"]["2.35"]["1e-08"]["main"]])
print("cut x=2.41 eps=1e-8:", [[round(x,4) for x in iv] for iv in post["cuts"]["2.41"]["1e-08"]["main"]])
print("cut x=2.41 eps=1e-5:", [[round(x,4) for x in iv] for iv in post["cuts"]["2.41"]["1e-05"]["main"]])
print("doc mentions merge at 1e-5:", "merge at `eps = 1e-5 V`" in doc)
