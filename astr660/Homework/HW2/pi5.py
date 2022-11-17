'''
This code is used to calculate pi using Monte Carlo method.
Author: Lin Yen-Hsing (NTHU) 2022.10.07
'''
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


def cal_pi(N):
    np.random.seed(111)
    x = np.random.random(N)
    np.random.seed(222)
    y = np.random.random(N)
    r = np.sqrt(x**2 + y**2)
/cluster/home/yhkuo/Computational-Astrophysics-2022/astr660/Homework/HW2/    PI = np.sum(r<1)/np.shape(r)[0]*4
    print(PI)
    return PI

def err(N):
    ERR = abs((cal_pi(N) - np.pi)/np.pi)
    return ERR

def expect(N,a):
    return a/np.sqrt(N)

N = np.logspace(1, 7, 20).astype(int)    #np.logspace(start, stop, number of sample)
Errs = list(map(err, N))                 #map(function, iterable) 

popt, pcov = curve_fit(expect, N, Errs)  # curve.fit(function, xdata, ydata), returns  

plt.loglog(N, Errs, label='Simulated data', linestyle='solid')                              # Make a plot with log scaling on both the x and y axis 
plt.loglog(N, expect(N, popt[0]), label='Expected $1/\sqrt{N}$ law', linestyle='dashedZZ')

plt.xlabel('$N$')
plt.ylabel('$\epsilon$')
plt.legend()
plt.tight_layout()                   #Adjust the location of axes/labels/titles
plt.savefig('pi_MC.png', dpi=300)

plt.show()








