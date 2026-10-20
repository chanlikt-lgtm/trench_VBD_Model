p="p2_front_conditioning.py"
s=open(p,encoding="utf-8").read()
s=s.replace('''    dw = np.column_stack([
        (4*L1-1)*dL1, (4*L2-1)*dL2, (4*L3-1)*dL3,
        4*(L1*dL2 + L2*dL1), 4*(L2*dL3 + L3*dL2), 4*(L3*dL1 + L1*dL3)])''','''    dw = np.column_stack([
        (4*L1-1)*dL1, (4*L2-1)*dL2, (4*L3-1)*dL3,
        4*(L1*dL2 + L2*dL1), 4*(L2*dL3 + L3*dL2), 4*(L3*dL1 + L1*dL3)])''')
s=s.replace('''    N = np.column_stack([
        L1*(2*L1-1), L2*(2*L2-1), L3*(2*L3-1),
        4*L1*L2, 4*L2*L3, 4*L3*L1])
    dN = np.column_stack([
        (4*L1-1)*dL1, (4*L2-1)*dL2, (4*L3-1)*dL3,
        4*(L1*dL2 + L2*dL1), 4*(L2*dL3 + L3*dL2), 4*(L3*dL1 + L1*dL3)])
    # reorder the local dofs to the detected edge permutation
    vals = ue.copy()
    return (N @ vals)[inside], (dN @ vals)[inside], ys[inside]''','''    Nc = {(0, 1): 4*L1*L2, (1, 2): 4*L2*L3, (2, 0): 4*L3*L1}
    dNc = {(0, 1): 4*(L1*dL2 + L2*dL1), (1, 2): 4*(L2*dL3 + L3*dL2),
           (2, 0): 4*(L3*dL1 + L1*dL3)}
    N = np.column_stack([L1*(2*L1-1), L2*(2*L2-1), L3*(2*L3-1),
                         Nc[order[0]], Nc[order[1]], Nc[order[2]]])
    dN = np.column_stack([(4*L1-1)*dL1, (4*L2-1)*dL2, (4*L3-1)*dL3,
                          dNc[order[0]], dNc[order[1]], dNc[order[2]]])
    return (N @ ue)[inside], (dN @ ue)[inside], ys[inside]''')
open(p,"w",encoding="utf-8").write(s)
print("patched ordering")
