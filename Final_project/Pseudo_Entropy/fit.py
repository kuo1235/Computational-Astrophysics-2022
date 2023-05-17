import pandas as pd
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

#---------------pick file first---------------------------
#df1 = pd.read_excel("N40_C1_1_C2_1_m1_1e1_m2_1e1.xlsx")
#df2 = pd.read_excel("N40_C1_1_C2_1_m1_1e-0_m2_1e-0.xlsx")
#df3 = pd.read_excel("N40_C1_1_C2_1_m1_1e-1_m2_1e-1.xlsx")
df4 = pd.read_excel("N40_C1_1_C2_100_m1_1e-2_m2_1e-2.xlsx")
df5 = pd.read_excel("N40_C1_1_C2_100_m1_1e-3_m2_1e-3.xlsx")
df6 = pd.read_excel("N40_C1_1_C2_100_m1_1e-4_m2_1e-4.xlsx")
df7 = pd.read_excel("N40_C1_1_C2_100_m1_1e-5_m2_1e-5.xlsx")
df8 = pd.read_excel("N40_C1_1_C2_100_m1_1e-6_m2_1e-6.xlsx")
#df9 = pd.read_excel("N40_C1_1_C2_1_m1_1e-7_m2_1e-7.xlsx")

#---------------------------------------------------------
#S1 = df1['Pseudo_Entropy']
#S2 = df2['Pseudo_Entropy']
#S3 = df3['Pseudo_Entropy']
S4 = df4['Pseudo_Entropy']
S5 = df5['Pseudo_Entropy']
S6 = df6['Pseudo_Entropy']
S7 = df7['Pseudo_Entropy']
S8 = df8['Pseudo_Entropy']
#S9 = df9['Pseudo_Entropy']

R = df4['n_sub']

#S1_total_list = np.array(S1)
#S2_total_list = np.array(S2)
#S3_total_list = np.array(S3)
S4_total_list = np.array(S4)
S5_total_list = np.array(S5)
S6_total_list = np.array(S6)
S7_total_list = np.array(S7)
S8_total_list = np.array(S8)
#S9_total_list = np.array(S9)

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
    ps2_array = np.array([])
    ps3_array = np.array([])
    ps4_array = np.array([])
    ps5_array = np.array([])
    ps6_array = np.array([])
    ps7_array = np.array([])
    ps8_array = np.array([])
    ps9_array = np.array([])

    n_array = np.array([])

    for i in range(n_ini,n_fin+1):
        #ps1_array = np.append(ps1_array, S1_total_list[i-1])
        #ps2_array = np.append(ps2_array, S2_total_list[i-1])
        #ps3_array = np.append(ps3_array, S3_total_list[i-1])
        ps4_array = np.append(ps4_array, S4_total_list[i-1])
        ps5_array = np.append(ps5_array, S5_total_list[i-1])
        ps6_array = np.append(ps6_array, S6_total_list[i-1])
        ps7_array = np.append(ps7_array, S7_total_list[i-1])
        ps8_array = np.append(ps8_array, S8_total_list[i-1])
        #ps9_array = np.append(ps9_array, S9_total_list[i-1])

        n_array = np.append(n_array, R_total_list[i-1])

    #print(ps_array)
    #print(n_array)

    #popt1, pcov1 = curve_fit(fit_func, n_array, ps1_array)
    #popt2, pcov2 = curve_fit(fit_func, n_array, ps2_array)
    #popt3, pcov3 = curve_fit(fit_func, n_array, ps3_array)
    popt4, pcov4 = curve_fit(fit_func, n_array, ps4_array)
    popt5, pcov5 = curve_fit(fit_func, n_array, ps5_array)
    popt6, pcov6 = curve_fit(fit_func, n_array, ps6_array)
    popt7, pcov7 = curve_fit(fit_func, n_array, ps7_array)
    popt8, pcov8 = curve_fit(fit_func, n_array, ps8_array)
    #popt9, pcov9 = curve_fit(fit_func, n_array, ps9_array)


    #plt.scatter(n_array, ps1_array, label='data1')
    #plt.scatter(n_array, ps2_array, label='data2')
    #plt.scatter(n_array, ps3_array, label='data3')
    plt.scatter(n_array, ps4_array, label='data4')
    plt.scatter(n_array, ps5_array, label='data5')
    plt.scatter(n_array, ps6_array, label='data6')
    plt.scatter(n_array, ps7_array, label='data7')
    plt.scatter(n_array, ps8_array, label='data8')
    #plt.scatter(n_array, ps9_array, label='data9')

    #plt.plot(n_array, fit_func(n_array, *popt1), label='m1=m2=1e1, fit: a=%6.4f, b=%6.4f' % tuple(popt1))
    #plt.plot(n_array, fit_func(n_array, *popt2), label='m1=m2=1e0, fit: a=%6.4f, b=%6.4f' % tuple(popt2))
    #plt.plot(n_array, fit_func(n_array, *popt3), label='m1=m2=1e-1, fit: a=%6.4f, b=%6.4f' % tuple(popt3))
    plt.plot(n_array, fit_func(n_array, *popt4), label='m1=m2=1e-2, fit: a=%6.4f, b=%6.4f' % tuple(popt4))
    plt.plot(n_array, fit_func(n_array, *popt5), label='m1=m2=1e-3, fit: a=%6.4f, b=%6.4f' % tuple(popt5))
    plt.plot(n_array, fit_func(n_array, *popt6), label='m1=m2=1e-4, fit: a=%6.4f, b=%6.4f' % tuple(popt6))
    plt.plot(n_array, fit_func(n_array, *popt7), label='m1=m2=1e-5, fit: a=%6.4f, b=%6.4f' % tuple(popt7))
    plt.plot(n_array, fit_func(n_array, *popt8), label='m1=m2=1e-6, fit: a=%6.4f, b=%6.4f' % tuple(popt8))
    #plt.plot(n_array, fit_func(n_array, *popt9), label='m1=m2=1e-7, fit: a=%6.4f, b=%6.4f' % tuple(popt9))

    plt.title('N=' + str(N) + '_C1=' + str(C1) + '_C2=' + str(C2)) 
    #plt.title('N=' + str(N) + '_m1=' + str(m1) + '_m2=' + str(m2))

    plt.xlabel('N_sub')
    plt.ylabel('Pseudo Entropy')
    plt.legend()
    #plt.savefig('N=' + str(N) + '_C1=' + str(C1) + '_C2=' + str(C2))
    #plt.savefig('N=' + str(N) + '_m1=' + str(m1) + '_m2=' + str(m2) )
    plt.savefig('N=' + str(N) + '_C1=1_C2=100')

    plt.show()

plot(40,9,16,1,100,1e-01,1e-05)
