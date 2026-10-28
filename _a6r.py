import json, os
p="p2_3_spatial_S0/CP-L2-A6_arc_region.json"
if not os.path.exists(p):
    print("NOT YET WRITTEN")
else:
    d=json.load(open(p))
    print("AREA_TOL", d["AREA_TOL"])
    for k,v in d["levels"].items():
        print(f"{k}: arc [{v['arc_lower']:.9e}, {v['arc_upper']:.9e}] w {v['arc_width']:.2e} cand {v['candidates']} {v['seconds']}s")
        print(f"     comp [{v['comp_lower']:.12f}, {v['comp_upper']:.12f}]")
