p="p2_plot_trench.py"
s=open(p,encoding="utf-8").read()
s=s.replace('X0, X1, Y0, Y1 = 2.24, 2.52, -0.34, 0.34','X0, X1, Y0, Y1 = 2.08, 2.56, -0.34, 0.34')
s=s.replace('("L1b  (h_local = 0.0055)", "p2_L1b_scaled.npz")','("L1b, h_loc 0.0055", "p2_L1b_scaled.npz")')
s=s.replace('("L1   (h_local = 0.0040)", "p2_L1_scaled.npz")','("L1, h_loc 0.0040", "p2_L1_scaled.npz")')
s=s.replace('figsize=(6.8, 8.0)','figsize=(7.2, 7.6)')
s=s.replace('axes[0].legend(fontsize=7.5, loc="upper left", framealpha=0.85)','axes[0].legend(fontsize=7.5, loc="lower left", framealpha=0.85)')
open(p,"w",encoding="utf-8").write(s)
print("patched")
