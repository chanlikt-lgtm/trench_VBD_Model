import json
for t in ("B0","B1"):
    d=json.load(open(f"p2_{t}_mesh.json"))
    print(f"{t}: nodes {d['nodes']} elems {d['elements']} Si {d['si_area']:.12f} "
          f"exact {d['exact_charged']:.12f} centroid {d['centroid_charged']:.12f} "
          f"near_band {d.get('near_band')}")
