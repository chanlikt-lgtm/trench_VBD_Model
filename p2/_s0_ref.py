import json
for t in ("B0","B1","B2"):
    d=json.load(open(f"p2_{t}_mesh.json"))
    print(f"{t}: nodes {d['nodes']} elems {d['elements']} h_global {d['h_global']}")
