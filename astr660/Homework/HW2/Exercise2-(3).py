'''
This code is used to calculate $\int_0^2 x^2 dx$ using Monte Carlo method.
Author: Lin Yen-Hsing (NTHU) 2022.10.07
'''
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


def cal_xsquare(N):
    np.random.seed(111)
    x = np.random.random(N)*2
    np.random.seed(222)
    y = np.random.random(N)*4
    result = np.sum(y<x**2)/np.shape(x)[0]*8
    print(result)
    return result

def err(N):
    ERR = abs((cal_xsquare(N) - 8/3)/(8/3))
    return ERR

def expect(N,a):
    return a/np.sqrt(N)

N = np.logspace(1, 7, 20).astype(int)
Errs = list(map(err, N))

popt, pcov = curve_fit(expect, N, Errs)

plt.loglog(N, Errs, label='Simulated data', linestyle='solid')
plt.loglog(N, expect(N, popt[0]), label='Expected $1/\sqrt{N}$ law', linestyle='dashed')

plt.xlabel('$N$')
plt.ylabel('$\epsilon$')
plt.legend()
plt.tight_layout()
plt.savefig('pi_MC.png', dpi=300)

plt.show()
