for src,dst in (("p2_B1X_osqp.py","p2_B2X_osqp.py"),("p2_B1X_polish.py","p2_B2X_polish.py")):
    s=open(src,encoding="utf-8").read()
    s=s.replace("p2_B1_p2mesh.npz","p2_B2_p2mesh.npz").replace("p2_B1X_osqp","p2_B2X_osqp").replace("p2_B1X_scaled","p2_B2X_scaled")
    s=s.replace("[B1X-A]","[B2X-A]").replace("[B1X-B]","[B2X-B]").replace("[B1X]","[B2X]")
    open(dst,"w",encoding="utf-8").write(s)
print("wrote p2_B2X_osqp.py, p2_B2X_polish.py")
