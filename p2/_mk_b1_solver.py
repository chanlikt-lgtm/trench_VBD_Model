for src,dst in (("p2_M2_osqp.py","p2_B1_osqp.py"),("p2_M2_polish.py","p2_B1_polish.py")):
    s=open(src,encoding="utf-8").read()
    s=s.replace("p2_M2_p2mesh.npz","p2_B1_p2mesh.npz").replace("p2_M2_osqp","p2_B1_osqp").replace("p2_M2_scaled","p2_B1_scaled")
    s=s.replace("[M2-A]","[B1-A]").replace("[M2-B]","[B1-B]").replace("[M2]","[B1]")
    open(dst,"w",encoding="utf-8").write(s)
print("wrote p2_B1_osqp.py, p2_B1_polish.py")
