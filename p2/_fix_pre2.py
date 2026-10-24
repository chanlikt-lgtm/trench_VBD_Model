p="p2_cpm2_pre.py"; s=open(p,encoding="utf-8").read()
old = '''        print(f"[{t}] midside dev {out[t]['midside_max_dev']:.2e}  "
              f"min angle {out[t]['quality']['min_angle']:.2f}  "
              f"aspect max {out[t]['quality']['aspect_max']:.3f}  "
              f"inverted {out[t]['quality']['n_inverted']}  "
              f"orphans {out[t]['orphan_nodes']}")
        print(f"     P2 dofs {out[t]['p2_dofs']}  charged {out[t]['charged_elements']}  "
              f"Bernstein rows {out[t]['bernstein_rows']}  donor delta {d:+.2e}")'''
new = '''        print(f"[{t}] midside dev {out[t]['midside_max_dev']:.2e}  "
              f"min angle {out[t]['quality']['min_angle']:.2f}  "
              f"aspect max {out[t]['quality']['aspect_max']:.3f}  "
              f"orphans {out[t]['orphan_nodes']}")
        print(f"     dofs raw {out[t]['p2_dofs_raw']} compacted {out[t]['p2_dofs_compacted']}"
              f"  charged {out[t]['charged_elements']}  rows {out[t]['bernstein_rows']}"
              f"  donor delta {d:+.2e}")
        for k, v in out[t]["orientation"].items():
            print(f"     orientation mat {k}: n={v['n']} pos={v['pos']} neg={v['neg']} zero={v['zero']}")'''
assert old in s
s = s.replace(old, new)
open(p,"w",encoding="utf-8").write(s)
print("patched print block")
