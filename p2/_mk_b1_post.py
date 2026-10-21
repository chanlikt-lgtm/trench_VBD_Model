s=open("p2_M2_post.py",encoding="utf-8").read()
s=s.replace("p2_M2_scaled.npz","p2_B1_scaled.npz").replace("p2_M2_post.json","p2_B1_post.json").replace("[M2 area]","[B1 area]")
open("p2_B1_post.py","w",encoding="utf-8").write(s)
print("wrote p2_B1_post.py")
