import json
d=json.load(open("../p2_3_spatial_S0/CP-L2_geometry.json"))
for t in ("L1","L2"):
    q=d[t]["quality"]
    print(f"{t}: min angle {q['min_angle']:.2f}  aspect max {q['aspect_max']:.3f}  "
          f"p99 {q['aspect_p99']:.3f}  >3 {q['n_gt3']}  >5 {q['n_gt5']}  zero {q['zero_area']}")
    print("   orientation:", json.dumps(q["orientation"]))
print("L2 nodes/elems:", d["L2"]["nodes"], d["L2"]["elements"], "Si", d["L2"]["si_area"])
print("regions:", json.dumps(d["L2"]["regions"], indent=1))
