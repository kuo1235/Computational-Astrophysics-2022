import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import linalg
from scipy.optimize import curve_fit
import time

#m1 = 1.0 * 10**(-3)
#m2 = 2.5 * 10**(-4)

#C1 = 1
#C2 = 10


def eigenK(coupling_matrix):
    eigenvalK, eigenvecK = linalg.eig(coupling_matrix)
    return eigenvalK, eigenvecK

def omega(C, m, k): #k can be n dimensional 
    
    omega_square = k**2 + C*(m**2)
    omega = omega_square**(1/2)

    return omega

def matrix(N, N_sub, C1, C2, m1, m2):

    #--------------------------------X, P, R matrix------------------------------------------- 
    X = np.zeros((N_sub, N_sub), dtype=complex)
    P = np.zeros((N_sub, N_sub), dtype=complex)
    R = np.zeros((N_sub, N_sub), dtype=complex)
    for i in range(1, N_sub+1):        
        for j in range(1, N_sub+1):
            x = 0
            p = 0
            r = 0

            for k in range(0, N):
                x_k = (0.5/N)*(2/(omega(C1, m1, k) + omega(C2, m2, k))) * np.cos(2 * np.pi * k * (i - j) / N)
                p_k = (0.5/N)*(2*omega(C1,m1,k)*omega(C2,m2,k))/(omega(C1,m1,k)+omega(C2,m2,k))*np.cos(2 * np.pi * k * (i-j) / N)
                r_k = (0.5j/N)*(omega(C2,m2,k)-omega(C1,m1,k))/(omega(C1,m1,k)+omega(C2,m2,k))*np.cos(2 * np.pi * k * (i-j) / N)
                
                #x_k = (0.5/N) * (2/(omega(m1, N, k) + omega(m2, N, k))) * np.exp(2j * np.pi * k * (i-j) / N)
                #p_k = (0.5/N) * (2*omega(m1, N, k)*omega(m2, N, k)) / (omega(m1, N, k) + omega(m2, N, k)) * np.exp(2j * np.pi * k * (i-j) / N)
                #r_k = (0.5/N) * (omega(m2, N, k)-omega(m1, N, k)) / (omega(m1, N, k) + omega(m2, N, k)) * np.exp(2j * np.pi * k * (i-j) / N)

                x += x_k 
                p += p_k
                r += r_k

            X[i-1][j-1] = x
            P[i-1][j-1] = p
            R[i-1][j-1] = r
    
    #print(X)
    #print(P)
    #print(R)
    
    #-------------------------------Gamma matrix---------------

    h1 = np.hstack((X,R))

    R_t = np.transpose(R)

    h2 = np.hstack((R_t,P))

    Gamma = np.vstack((h1,h2))

    #--------------------------------J matrix----------------

    o = np.zeros((N_sub,N_sub))

    I = np.identity(N_sub)

    h1 = np.hstack((o,I))
    h2 = np.hstack((-I,o))

    J = np.vstack((h1,h2))
    
    #print(J)

    #print(Gamma)

    #---------------------------iJG matrix-----------------------

    iJG = 1j * np.matmul(J,Gamma)

    #print(iJG)

    eigenvalC, eigenvecC = eigenK(iJG)
    
   # print(eigenvalC)
    
    eigenvalC = eigenvalC[eigenvalC >= 0]

    #print(eigenvalC)


#matrix(10,5)

#def mode_entropy(eigenvalC):
    #size = np.size(eigenvalC)

    mode_entropy=0

    for i in range(0,N_sub):
        #if abs(np.real(eigenvalC[i])-0.5) < 10**-10:
            #modeS = 0
        #else:
        modeS = (eigenvalC[i]+0.5) * np.log(eigenvalC[i]+0.5) - (eigenvalC[i]-0.5) * np.log(eigenvalC[i]-0.5)
        #modeS = (np.real(eigenvalC[i]+0.5)) * np.log(np.real(eigenvalC[i])+0.5) - (np.real(eigenvalC[i])-0.5) * np.log(np.real(eigenvalC[i])-0.5)
 
        #print(modeS) 
         
        mode_entropy += modeS
    
    #print(mode_entropy)

    return mode_entropy

#matrix(200,40, 1, 10, 5.0*10**(-8), 1.0*10**(-8)) 

def fit_func(n_sub, a, b, c):
    #return a * np.log((200/np.pi) * np.sin( np.pi * n_sub/200)) + b
    return a * 4 * np.pi * n_sub**2 + b * np.log(n_sub**2) + c 

def plot(N, n_ini, n_fin, step, C1, C2, m1, m2): 
    
    d = int((n_fin - n_ini)/step)
    #print(d)

    ps_array = np.array([])
    n_array  = np.array([])

    start_time = time.time()

    for i in range(0, d+1):

        ps = matrix(N, n_ini + step*i, C1, C2, m1, m2)

        ps_array = np.append(ps_array, ps)
        n_array  = np.append(n_array, n_ini + step*i) 

    popt, pcov = curve_fit(fit_func, n_array, ps_array)
    print(popt)

    print('Time used: {} sec'.format(time.time()-start_time))

    col1 = "n_sub"
    col2 = "Pseudo_Entropy"
    col3 = "a"
    col4 = "b"
    col5 = "c"
    
    arealaw = pd.DataFrame({col1:n_array, col2:ps_array, col3:popt[0], col4:popt[1], col5:popt[2]})
    arealaw.to_excel('N200.xlsx', sheet_name='sheet1', index=False)

    plt.scatter(4* np.pi * n_array**2 * 0.001, ps_array, label='data')     
    plt.plot(4* np.pi * n_array**2 * 0.001, fit_func(n_array, *popt), label='fit: a=%6.4f, b=%6.4f, c=%6.4f' % tuple(popt), color='r')

    plt.title('N=' + str(N) + ', C1=' + str(C1) + ', C2=' + str(C2) + ', m1=' + str(m1) + ', m2=' + str(m2))
    plt.xlabel('Area/1000')
    plt.ylabel('Pseudo Entropy')
    plt.legend()
    plt.savefig('QFT_minimal_N'+str(N)+'_m1_'+str(m1)+', m2_'+str(m2))

    plt.show()


plot(200, 5, 195, 5, 1, 10, 1, 1)

def S(m1, m2, l):
    
    l_list = np.array([])
    S_list = np.array([])

    for i in range (1,l+1):

        l_list = np.append(l_list, i)

        S = 0.5 * np.log( -(m1**2 * np.log(m1*(1/i)) - m2**2 * np.log(m2*(1/i))) / (m1**2 - m2**2) )    
        S_list = np.append(S_list, S)
        
        print("l="+str(i)+": " +str(S))

    


    #plt.scatter(l_list, S_list, label='Analytical Solution')     
    plt.plot(l_list, S_list, label='Analytical Solution')

    plt.title('Pseudo Entropy')
    plt.xlabel('l')
    plt.ylabel('PS')
    plt.legend()

    plt.show()

#S(m1, m2, 200)




