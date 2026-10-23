p="p2_s0_tests.py"; s=open(p,encoding="utf-8").read()
s=s.replace("I = D.int_N(nodes[tris[0]][:, :3], -1e9, area)",
            "Pv = nodes[tris[0][:3]]\n    I = D.int_N(Pv, -1e9, area)")
s=s.replace("I2 = D.int_N(nodes[tris[0]][:, :3], 1e9, area)",
            "I2 = D.int_N(Pv, 1e9, area)")
s=s.replace("    P = nodes[tris[0]][:, :3]\n    xb = 0.7 + 0.4",
            "    P = Pv\n    xb = 0.7 + 0.4")
open(p,"w",encoding="utf-8").write(s)
print("patched point slicing")
