s=open("p2_B1_post.py",encoding="utf-8").read()
s=s.replace("p2_B1_scaled.npz","p2_B0_scaled.npz").replace("p2_B1_post.json","p2_B0_post.json").replace("[B1 area]","[B0 area]")
open("p2_B0_post.py","w",encoding="utf-8").write(s)
print("wrote p2_B0_post.py")
