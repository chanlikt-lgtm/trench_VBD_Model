p="p2_s0_tests.py"; s=open(p,encoding="utf-8").read()

# FIX 1: the fixture triangle had all mid-side dofs aliased to node 0.
s = s.replace('''    nodes = np.array([[at, 0.0], [at+1.0, 0.0], [at, 1.0]])
    tris = np.array([[0, 1, 2, 0, 0, 0]])
    mat = np.array([M.SI])''',
'''    nodes = np.array([[at, 0.0], [at+1.0, 0.0], [at, 1.0],
                      [at+0.5, 0.0], [at+0.5, 0.5], [at, 0.5]])
    tris = np.array([[0, 1, 2, 3, 4, 5]])
    mat = np.array([M.SI])''')

# FIX 2: P2 vertex shape functions are sign-changing (N = L(2L-1)), so
# int N_i over a partial region may legitimately be negative. Assert the
# conserved total and finiteness instead of non-negativity.
s = s.replace('''    neg = float(min(0.0, I3.min()))

    RES["donor"] = dict(full_identity_max_err=e1, body_side_max=I2.max(),
                        cut_sum_err=float(e3), cut_min_component=neg,
                        cut_area=Aclip, passed=bool(e1 < 1e-12 and I2.max() == 0.0
                                                    and e3 < 1e-12 and neg >= 0.0))''',
'''    neg = float(I3.min())

    RES["donor"] = dict(full_identity_max_err=e1, body_side_max=I2.max(),
                        cut_sum_err=float(e3), cut_min_component=neg,
                        cut_area=Aclip,
                        note="negative vertex components are expected: N=L(2L-1) is sign-changing",
                        passed=bool(e1 < 1e-12 and I2.max() == 0.0
                                    and e3 < 1e-12 and np.isfinite(I3).all()))''')
open(p,"w",encoding="utf-8").write(s)
print("fixtures corrected")
