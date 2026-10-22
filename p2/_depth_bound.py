import math
for tag,amin in (("B0",1.641e-05),("B1",3.664e-06),("B0X",1.641e-05),("B1X",3.664e-06)):
    d=math.log(amin/1e-13)/math.log(4)
    print(f"{tag}: min element area {amin:.3e} -> max terminating depth {d:.2f} (MAX_DEPTH 24)")
