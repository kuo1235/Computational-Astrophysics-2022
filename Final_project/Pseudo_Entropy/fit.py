import pandas as pd
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

#---------------pick file first---------------------------
df1 = pd.read_excel("N40_C1_10_C2_1_m1_1e-5_m2_1e-5.xlsx")
#df2 = pd.read_excel("N40_C1_1_C2_1_m1_5e-5_m2_1e-5.xlsx")
#df3 = pd.read_excel("N20_C1_1_C2_1_m1_1e-4_m2_1e-5.xlsx")

#---------------------------------------------------------
S1 = df1['Pseudo_Entropy']
#S2 = df2['Pseudo_Entropy']
#S3 = df3['Pseudo_Entropy']
R = df1['n_sub']

S1_total_list = np.array(S1)
#S2_total_list = np.array(S2)
#S3_total_list = np.array(S3)
R_total_list = np.array(R)


def fit_func(n_sub, a, b):

    return a * (40/np.pi) * np.sin( (np.pi * n_sub / 40)) + b
    #return a * (20/np.pi)**2 * np.sin( (np.pi * n_sub / 20))**2 + b

    #return a * np.log((40/np.pi) * np.sin( np.pi * n_sub / 40)) + b
    #return a * np.log(n_sub) + b

def r_squared(ydata,calc_ydata):
###计算拟合度
    ydata = [float(x) for x in ydata]
    calc_ydata = [float(x) for x in calc_ydata]
    residual = np.array(ydata)-np.array(calc_ydata)
#总平方和sst=ssr+sse
    sst=np.sum((ydata-np.mean(ydata))**2)
#回归平方和
    ssr=np.sum((calc_ydata-np.mean(ydata))**2)
#残差平方和
    sse=np.sum(residual**2)
#print(ssr,sse,sst)
    r_squared = ssr/sst
    return r_squared

def plot(N, n_ini, n_fin, C1, C2, m1, m2):
 
    ps1_array = np.array([])
    #ps2_array = np.array([])
    #ps3_array = np.array([])
    n_array = np.array([])

    for i in range(n_ini,n_fin+1):
        ps1_array = np.append(ps1_array, S1_total_list[i-1])
        #ps2_array = np.append(ps2_array, S2_total_list[i-1])
        #ps3_array = np.append(ps3_array, S3_total_list[i-1])

        n_array = np.append(n_array, R_total_list[i-1])

    #print(ps_array)
    #print(n_array)

    popt1, pcov1 = curve_fit(fit_func, n_array, ps1_array)
    #popt2, pcov2 = curve_fit(fit_func, n_array, ps2_array)
    #popt3, pcov3 = curve_fit(fit_func, n_array, ps3_array)

    plt.scatter(n_array, ps1_array, label='data1')
    #plt.scatter(n_array, ps2_array, label='data2')
    #plt.scatter(n_array, ps3_array, label='data3')

    plt.plot(n_array, fit_func(n_array, *popt1), label='fit: a=%6.4f, b=%6.4f' % tuple(popt1), color='r')
    #plt.plot(n_array, fit_func(n_array, *popt2), label='fit: a=%6.4f, b=%6.4f' % tuple(popt2), color='r')
    #plt.plot(n_array, fit_func(n_array, *popt3), label='fit: a=%6.4f, b=%6.4f' % tuple(popt3), color='r')

    plt.title('N=' + str(N) + '_C1=' + str(C1) + '_C2=' + str(C2)) 
    #plt.title('N=' + str(N) + '_m1=' + str(m1) + '_m2=' + str(m2))

    plt.xlabel('N_sub')
    plt.ylabel('Pseudo Entropy')
    plt.legend()
    plt.savefig('N=' + str(N) + '_C1=' + str(C1) + '_C2=' + str(C2) +  '_m1=' + str(m1) + '_m2=' + str(m2) )

    plt.show()

plot(40,9,16,10,1,1e-5,1e-5)
