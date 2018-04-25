import numpy as np
RHS_True=np.array([7+0.j,20.,-15,-3.,16,-27])
print(RHS_True.dtype)

RHS_True.tofile('tofile.dat')
