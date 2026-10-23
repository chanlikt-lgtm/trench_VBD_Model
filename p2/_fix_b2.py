p="p2_mesh_B2.py"; s=open(p,encoding="utf-8").read()
s=s.replace('gmsh.option.setNumber("Mesh.MeshSizeMin", H_TIP)',
            'gmsh.option.setNumber("Mesh.MeshSizeMin", min(H_TIP, H_BAND))')
s=s.replace('print("wrote p2_B1_p2mesh.npz / p2_B1_mesh.json")',
            'print(f"wrote {PREFIX}_p2mesh.npz / {PREFIX}_mesh.json")')
open(p,"w",encoding="utf-8").write(s)
print("patched MeshSizeMin floor and output message")
