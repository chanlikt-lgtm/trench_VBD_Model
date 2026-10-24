import json
a=json.load(open("p2_3_spatial_S0/CP-M2-A_manifest.json"))
print("manifest status:", a["status"], "| osqp", a["software"]["osqp"], "| gmsh", a["software"]["gmsh"])
for t in ("L0","L1"):
    b=json.load(open(f"p2_3_spatial_S0/CP-M2-B_{t}.json"))
    m=a["mesh"][t]
    ok = (b["dofs_total"]==m["compacted"] and b["charged_elements"]==m["charged"]
          and b["bernstein_rows"]==m["bernstein_rows"])
    print(f"{t}: A2 correspondence {ok} | peak {b['peak_rss_mb']:.0f} MB | "
          f"donor-exact {b['donor_sum']-1.5752528820380198:+.2e}")
tot=24467; thr=0.7*tot
print(f"70% threshold {thr:.0f} MB | projected L1 solve ~1240 MB | fits: {1240<thr}")
print(f"available at measurement: L0 4486 MB, L1 4213 MB | margin over 1600 MB projection: {4213-1600} MB")
