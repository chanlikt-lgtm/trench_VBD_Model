p="p2_mesh_S0.py"; s=open(p,encoding="utf-8").read()
s=s.replace('    os.makedirs("p2_3_spatial_S0/levels", exist_ok=True)\n    out = os.path.join("..","p2_3_spatial_S0","levels")',
            '    out = os.path.join("..","p2_3_spatial_S0","levels")\n    os.makedirs(out, exist_ok=True)')
open(p,"w",encoding="utf-8").write(s)
print("patched makedirs order")
