for src,dst in (("p2_L1b_osqp.py","p2_M2_osqp.py"),("p2_L1b_polish.py","p2_M2_polish.py")):
    s=open(src,encoding="utf-8").read()
    s=s.replace("p2_L1b_p2mesh.npz","p2_M2_p2mesh.npz").replace("p2_L1b_osqp","p2_M2_osqp").replace("p2_L1b_scaled","p2_M2_scaled")
    s=s.replace("[L1b-A]","[M2-A]").replace("[L1b-B]","[M2-B]").replace("[L1b]","[M2]")
    s=s.replace("L1b Step A","M2 Step A").replace("L1b Step B/C","M2 Step B/C")
    open(dst,"w",encoding="utf-8").write(s)
print("wrote p2_M2_osqp.py, p2_M2_polish.py")
