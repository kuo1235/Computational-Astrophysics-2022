import numpy as np
from scipy import linalg
from scipy.sparse import diags
import timeit

A = diags([1,-4,6,-4,1],[-2,-1,0,1,2], shape=(100,100)).toarray()

A[0][0]=9

A[98][98]=5
A[98][99]=-2

A[99][98]=-2
A[99][99]=1


b = np.full((100,1),1)

def solver_bandmatrix(A,b,N):
    ab=np.full((5,N), 0)
    for i in np.arange(N):
        if i+2<N:
            ab[0][i+2]=A[i][i+2]
            ab[4][i]=A[i+2][i]
        if i+1<N:
            ab[1][i+1]=A[i][i+1]
            ab[3][i]=A[i+1][i]    
        ab[2][i]=A[i][i]
    
    x = linalg.solve_banded((2, 2), ab, b)
    
    return x    

%timeit solver_bandmatrix(A,b,100)

x = solver_bandmatrix(A,b,100)

print(x)
