for tag in ("B0","B1"):
    for kind in ("osqp","polish"):
        src=f"p2_{tag}_{kind}.py"; dst=f"p2_{tag}X_{kind}.py"
        s=open(src,encoding="utf-8").read()
        s=s.replace("from p2_nested_consistency import assemble_p2",
                    "from p2_nested_consistency import assemble_p2\nfrom p2_exact_donor import exact_load")
        s=s.replace("    K, f = assemble_p2(nodes, tris, mat, charged)",
                    "    K, _f0 = assemble_p2(nodes, tris, mat, np.zeros(len(tris), bool))\n    f, _ftot = exact_load(nodes, tris, mat)")
        s=s.replace(f"p2_{tag}_osqp", f"p2_{tag}X_osqp").replace(f"p2_{tag}_scaled", f"p2_{tag}X_scaled")
        s=s.replace(f"[{tag}-A]", f"[{tag}X-A]").replace(f"[{tag}-B]", f"[{tag}X-B]").replace(f"[{tag}]", f"[{tag}X]")
        open(dst,"w",encoding="utf-8").write(s)
        print("wrote", dst)
