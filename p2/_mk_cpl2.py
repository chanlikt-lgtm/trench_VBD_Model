s = open("p2_cpl1.py", encoding="utf-8").read()
rep = [
 ('TAG = "S0_L1"', 'TAG = "S0_L2"'),
 ('MESH = "../p2_3_spatial_S0/levels/S0_L1_p2mesh.npz"',
  'MESH = "../p2_3_spatial_S0/levels/S0_L2_p2mesh.npz"'),
 ('MESH_SHA = "7FA1C69DA097612960A2668335A014894A0DA4D9FF5C96D623E10FCAE0A3D1BB"',
  'MESH_SHA = "7E6CF12400B0A1308DC23BC7DEC89DA96EF30766784C90D0CB433166523C9136"'),
 ('EXPECT_DOFS = 459288', 'EXPECT_DOFS = 913595'),
 ('EXPECT_ROWS = 709740', 'EXPECT_ROWS = 1413522'),
 ('CP-L1_L1.json', 'CP-L2_L2.json'),
 ('CP-L1_L1_checkpoint.npz', 'CP-L2_L2_checkpoint.npz'),
 ('CP-L1_L1_polished.npz', 'CP-L2_L2_polished.npz'),
 ('"""CP-L1: S0-L1 production solve', '"""CP-L2: S0-L2 production solve'),
]
for a,b in rep:
    assert a in s, a[:45]
    s = s.replace(a,b)
open("p2_cpl2.py","w",encoding="utf-8").write(s)
print("wrote p2_cpl2.py")
