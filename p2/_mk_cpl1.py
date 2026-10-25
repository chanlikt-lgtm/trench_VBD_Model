s = open("p2_cpm2c_l0.py", encoding="utf-8").read()
rep = [
 ('TAG = "S0_L0"', 'TAG = "S0_L1"'),
 ('MESH = "../p2_3_spatial_S0/levels/S0_L0_p2mesh.npz"',
  'MESH = "../p2_3_spatial_S0/levels/S0_L1_p2mesh.npz"'),
 ('MESH_SHA = "58ABC85E86E3D17CD6A76ADE85AE90137FC12C58D291EF2C07B1A4759712D241"',
  'MESH_SHA = "7FA1C69DA097612960A2668335A014894A0DA4D9FF5C96D623E10FCAE0A3D1BB"'),
 ('EXPECT_DOFS = 231274', 'EXPECT_DOFS = 459288'),
 ('EXPECT_ROWS = 355932', 'EXPECT_ROWS = 709740'),
 ('CP-M2-C_L0.json', 'CP-L1_L1.json'),
 ('CP-M2-C_L0_checkpoint.npz', 'CP-L1_L1_checkpoint.npz'),
 ('CP-M2-C_L0_polished.npz', 'CP-L1_L1_polished.npz'),
 ('"""CP-M2-C: S0-L0 production solve under the frozen manifest.',
  '"""CP-L1: S0-L1 production solve under the frozen manifest.'),
]
for a,b in rep:
    assert a in s, a[:40]
    s = s.replace(a,b)
open("p2_cpl1.py","w",encoding="utf-8").write(s)
print("wrote p2_cpl1.py")
