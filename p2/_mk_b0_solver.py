for src,dst in (("p2_B1_osqp.py","p2_B0_osqp.py"),("p2_B1_polish.py","p2_B0_polish.py")):
    s=open(src,encoding="utf-8").read()
    s=s.replace("p2_B1_p2mesh.npz","p2_B0_p2mesh.npz").replace("p2_B1_osqp","p2_B0_osqp").replace("p2_B1_scaled","p2_B0_scaled")
    s=s.replace("[B1-A]","[B0-A]").replace("[B1-B]","[B0-B]").replace("[B1]","[B0]")
    open(dst,"w",encoding="utf-8").write(s)
print("wrote p2_B0_osqp.py, p2_B0_polish.py")
