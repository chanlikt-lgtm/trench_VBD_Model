p="p2_plot_trench.py"
s=open(p,encoding="utf-8").read()
s=s.replace('    fig, axes = plt.subplots(2, 1, figsize=(8.2, 7.4), sharex=True, sharey=True)',
            '    fig, axes = plt.subplots(2, 1, figsize=(6.8, 8.0), sharex=True, sharey=True,\n                             constrained_layout=True)')
s=s.replace('        ax.set_aspect("equal")\n        ax.set_title(f"potential φ in the trench window — {lab}", fontsize=9.5)',
            '        ax.set_aspect("equal", adjustable="box", anchor="C")\n        ax.set_title(f"φ in the trench window — {lab}", fontsize=9.5)')
s=s.replace('    axes[-1].set_xlabel("x (µm)")\n    fig.tight_layout()\n    fig.savefig("fig6_trench_phi.png")',
            '    axes[-1].set_xlabel("x (µm)")\n    fig.savefig("fig6_trench_phi.png")')
s=s.replace('    axes[0].legend(fontsize=8, loc="lower right")',
            '    axes[0].legend(fontsize=7.5, loc="upper left", framealpha=0.85)')
s=s.replace('    fig.suptitle("Obstacle variable across the cut (element samples)", fontsize=10)\n    fig.tight_layout()',
            '    fig.suptitle("Obstacle variable across the cut (element samples)", fontsize=10)\n    fig.tight_layout(rect=(0, 0, 1, 0.94))')
open(p,"w",encoding="utf-8").write(s)
print("patched")
