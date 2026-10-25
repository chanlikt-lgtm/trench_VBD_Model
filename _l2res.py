import json
r=json.load(open("p2_3_spatial_S0/CP-L2_resources.json"))
g=json.load(open("p2_3_spatial_S0/CP-L2_geometry.json"))
print("band_sha_match:", g["band_sha_match"])
print("mesh sha L2:", g["mesh_sha"]["L2"])
print("L2 dofs:", r["dofs_total"], "free:", r["dofs_free"], "fixed:", r["dofs_fixed"])
print("charged:", r["charged_elements"], "rows:", r["bernstein_rows"])
print("K nnz", r["K"]["nnz"], "%.1f MB" % r["K"]["mb"],
      "| Kff nnz", r["Kff"]["nnz"], "%.1f MB" % r["Kff"]["mb"],
      "| C_F nnz", r["C_F"]["nnz"], "%.1f MB" % r["C_F"]["mb"])
print("donor", r["donor_sum"], "| delta vs exact %.2e" % (r["donor_sum"]-1.5752528820380198))
print("assembly s:", r["seconds"], "peak RSS", r["peak_rss_mb"], "MB",
      "| RAM avail", r["ram"]["avail_at_start_mb"])
