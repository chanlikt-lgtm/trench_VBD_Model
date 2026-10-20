import re
p = "p2_make_figs.py"
s = open(p, encoding="utf-8").read()

# fig1: give the log axis room so the L1 label is not clipped
s = s.replace('    ax.set_xscale("log")\n    ax.invert_xaxis()\n    ax.set_xlabel("realised element size',
              '    ax.set_xscale("log")\n    ax.invert_xaxis()\n    ax.margins(x=0.18)\n    ax.set_xlabel("realised element size')

# fig3/4: draw each interval as a proper broken bar with a legend
old = s[s.index("def fig_front("):s.index("def fig_edges(")]
new = '''def fig_front(x0, fname, title):
    data = cuts(x0, 1e-8)
    names = [m[0] for m in MESHES]
    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    for i, name in enumerate(names):
        spans = [(a, b - a) for (a, b) in data[name]]
        ax.broken_barh(spans, (i - 0.28, 0.56), facecolors=COL[name],
                       edgecolors="k", linewidth=0.6, alpha=0.55)
        lo = min(a for a, b in data[name]); hi = max(b for a, b in data[name])
        ax.plot([i, i], [lo, hi], color=COL[name], lw=1.0, ls=":")
        ax.annotate(f"{len(spans)} interval{'s' if len(spans) > 1 else ''}",
                    (i, hi), textcoords="offset points", xytext=(0, 5),
                    ha="center", fontsize=7.5)
    ax.set_xticks(range(len(names)))
    ax.set_xticklabels([f"{n}\\nh={m[2]:.4f}" for n, m in zip(names, MESHES)])
    ax.set_xlim(-0.6, len(names) - 0.4)
    ax.set_ylabel("y at the cut (µm)")
    ax.set_title(title, fontsize=10)
    ax.axhline(0, color="k", lw=0.6, ls=":")
    ax.plot([], [], "s", color="0.45", alpha=0.55, label="neutral interval (u < ε)")
    ax.legend(fontsize=7.5, loc="lower right")
    fig.tight_layout()
    fig.savefig(fname)
    print(fname)


'''
s = s.replace(old, new)

# fig5: mesh labels become axis ticks instead of crowded annotations
s = s.replace('''        for h, v, n in zip(hs, lower, names):
            ax.annotate(n, (h, v), textcoords="offset points", xytext=(4, -11),
                        fontsize=7.5, color="#1f6fb4")
        ax.set_xscale("log")
        ax.invert_xaxis()''', '''        ax.set_xscale("log")
        ax.invert_xaxis()
        ax.set_xticks(hs)
        ax.set_xticklabels([f"{n}\\n{h:.4f}" for n, h in zip(names, hs)], fontsize=7.5)
        ax.margins(x=0.12)''')
open(p, "w", encoding="utf-8").write(s)
print("patched")
