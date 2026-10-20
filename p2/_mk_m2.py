p="p2_mesh_L1b.py"; s=open(p,encoding="utf-8").read()
s=s.replace("THICK = 0.012","THICK = 0.030")
s=s.replace('PREFIX = "p2_L1b"','PREFIX = "p2_M2"')
s=s.replace("L1b mesh:","M2 mesh:")
s=s.replace('print(f"[L1b]','print(f"[M2]')
s=s.replace("Lower-contrast diagnostic variant of L1 (H_LOCAL 0.0040 -> 0.0055).",
            "Mesh-layout control: same core size as L1b, transition spread 0.012 -> 0.030.")
open("p2_mesh_M2.py","w",encoding="utf-8").write(s)
print("wrote p2_mesh_M2.py")
