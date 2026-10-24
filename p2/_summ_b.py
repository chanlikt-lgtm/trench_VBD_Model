import json
for t in ("L0","L1"):
    d=json.load(open(f"../p2_3_spatial_S0/CP-M2-B_{t}.json"))
    print(f"{t}: dofs {d['dofs_total']} free {d['dofs_free']} | K {d['K']['nnz']} nnz {d['K']['mb']:.1f} MB | "
          f"Kff {d['Kff']['nnz']} nnz {d['Kff']['mb']:.1f} MB | C_F {d['C_F']['nnz']} nnz {d['C_F']['mb']:.1f} MB")
    print(f"    peak RSS {d['peak_rss_mb']:.0f} MB | donor {d['donor_sum']:.12f} | "
          f"t={d['seconds']['total']}s | RAM total {d['ram']['total_mb']} avail {d['ram']['avail_at_start_mb']}")
