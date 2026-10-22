s=open("p2_B1_post.py",encoding="utf-8").read()
s=s.replace("p2_B1_scaled.npz","p2_B0X_scaled.npz").replace("p2_B1_post.json","p2_B0X_post.json").replace("[B1 area]","[B0X area]")
open("p2_B0X_post.py","w",encoding="utf-8").write(s)
print("wrote p2_B0X_post.py")
