import numpy as np
from scipy import linalg
from scipy.sparse import diags

A = diags([1,-4,6,-4,1],[-2,-1,0,1,2], shape=(100,100)).toarray()

A[0][0]=9

A[98][98]=5
A[98][99]=-2

A[99][98]=-2
A[99][99]=1


b = np.full((100,1),1)

lu, piv = linalg.lu_factor(A)

x = linalg.lu_solve((lu, piv), b)

print(x)
