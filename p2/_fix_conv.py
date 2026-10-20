import io
# 1) canonical extractor: pass barycentric weights of v0 and v1
p="p2_exact_metrics.py"
s=open(p,encoding="utf-8").read()
old="            if quad:\n                return float(np.dot(u[tris[e]], A.shape_vals(l1, l2)))"
new="            if quad:\n                # shape_vals expects the weights of vertices 0 and 1;\n                # (l1,l2) are the coefficients of (P1-P0),(P2-P0).\n                return float(np.dot(u[tris[e]], A.shape_vals(1.0-l1-l2, l1)))"
assert old in s, "extractor anchor not found"
s=s.replace(old,new)
open(p,"w",encoding="utf-8").write(s)

# 2) diagnostic evaluator: same correction, plus derivative chain
p="p2_front_conditioning.py"
s=open(p,encoding="utf-8").read()
old="    L1 = L[:, 0]; L2 = L[:, 1]; L3 = 1.0 - L1 - L2"
new="    a = L[:, 0]; b = L[:, 1]\n    L1 = 1.0 - a - b; L2 = a; L3 = b"
assert old in s, "evaluator anchor not found"
s=s.replace(old,new)
old2="    dL1 = Minv[0, 1]; dL2 = Minv[1, 1]; dL3 = -(dL1 + dL2)"
new2="    da = Minv[0, 1]; db = Minv[1, 1]\n    dL1 = -(da + db); dL2 = da; dL3 = db"
assert old2 in s, "derivative anchor not found"
s=s.replace(old2,new2)
open(p,"w",encoding="utf-8").write(s)
print("patched both evaluation sites")
