p="p2_cpm2c_l0.py"; s=open(p,encoding="utf-8").read()
s=s.replace("    osqp.matrix_type = osqp.CSC\n","")
open(p,"w",encoding="utf-8").write(s)
print("removed matrix_type line")
