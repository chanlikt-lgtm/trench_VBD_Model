s=open("p2_L1_post.py",encoding="utf-8").read()
s=s.replace("p2_L1_scaled.npz","p2_M2_scaled.npz").replace("p2_L1_post.json","p2_M2_post.json").replace("[L1 area]","[M2 area]")
open("p2_M2_post.py","w",encoding="utf-8").write(s)
print("wrote p2_M2_post.py")
