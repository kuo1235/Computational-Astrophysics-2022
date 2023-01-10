import pandas as pd
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

#---------------pick file first---------------------------
df = pd.read_excel("N200_1n100_L10000.xlsx")
#---------------------------------------------------------
S = df['S_total_r']
R = df['R']

S_total_list = np.array(S)
R_total_list = np.array(R)


def fit_func(x,a,b,c):
    return a * 4*np.pi * x**2 + b * np.log(x**2) + c

def plot(N, n_ini, n_fin, l_max):
 
    S_list = np.array([])
    R_list = np.array([])

    for i in range(n_ini,n_fin+1):
        S_list = np.append(S_list, S_total_list[i-1])
        R_list = np.append(R_list, R_total_list[i-1])


    popt, pcov = curve_fit(fit_func, R_list, S_list)
    print(popt)

    plt.scatter(R_list, S_list, label='data')     
    plt.plot(R_list, fit_func(R_list, *popt), label='fit: a=%6.4f, b=%6.4f, c=%6.4f' % tuple(popt), color='r')

    plt.title('N'+str(N)+'_' + str(n_ini) + '<n<' + str(n_fin) + ', l_max=' + str(l_max))
    plt.xlabel('Area')
    plt.ylabel('S')
    plt.legend()
    plt.savefig('N'+str(N)+'_'+str(n_ini)+'n'+str(n_fin)+'_L'+str(l_max))

    plt.show()


plot(200,1,64,10000)
