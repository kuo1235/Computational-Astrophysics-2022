import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import linalg
from scipy.optimize import curve_fit
import time


def eigenK(coupling_matrix):
    eigenvalK, eigenvecK = linalg.eig(coupling_matrix)
    return eigenvalK, eigenvecK

def coupling_matrix(N, L, C, m):

    #Coupling Matrix---------------

    #epsilon = 1

    couplingK = np.zeros((N, N))
    for i in range(1, N):
        couplingK[i-1][i-1] = 0.5*( (C * m**2 + np.pi**2 * L*(L+1)) / (N**2 * np.sin( np.pi * (i+0.5)/N )**2) + np.sin( np.pi * (i-0.5)/N )**2 / np.sin( np.pi * (i+0.5)/N)**2  + 1)
        couplingK[i-1][i] = - 0.5 * np.sin( np.pi * (i-0.5)/N ) / np.sin( np.pi * (i+0.5)/N ) 
        couplingK[i][i-1] = couplingK[i-1][i]
        
    couplingK[N-1][N-1] = 0.5 * ( (C * m**2 + np.pi**2 * L*(L+1)) / (N**2 * np.sin( np.pi * (N+0.5)/N )**2) + np.sin( np.pi * (N-0.5)/N )**2 / np.sin( np.pi * (N+0.5)/N)**2 + 1 )

    return couplingK   

def omega(N, L, C, m):

    #Eigenvalue------------------------

    eigenvalK, eigenvecK = eigenK(coupling_matrix(N, L, C, m))

    print(eigenvalK)

    return eigenvalK


#omega(5, 10 ,1, 1)

def matrix(N, N_sub, L_max, C1, C2, m1, m2):
        
    
    #-------------------------------OmegaK---------------------------------------------------

    #eigenvalK1 = eigenK(coupling_matrix(N , L, C1, m1)
    
    #eigenvalK2 = eigenK(coupling_matrix(N , L, C2, m2)

    #--------------------------------X, P, R matrix------------------------------------------- 
    X = np.zeros((N_sub, N_sub), dtype=complex)
    P = np.zeros((N_sub, N_sub), dtype=complex)
    R = np.zeros((N_sub, N_sub), dtype=complex)
    for i in range(1, N_sub+1):        
        for j in range(1, N_sub+1):
            x = 0
            p = 0
            r = 0

            for L in range(0, L_max+1):
                    
                eigenvalK1, eigenvecK1 = eigenK(coupling_matrix(N , L, C1, m1))
                eigenvalK2, eigenvecK2 = eigenK(coupling_matrix(N , L, C2, m2))
              
                x_l = 0
                p_l = 0
                r_l = 0

                for k in range(0, N):

                    omega1 = eigenvalK1[k]
                    omega2 = eigenvalK2[k]

                    x_k = (0.5/N) * (2 / (omega1 + omega2)) * np.cos(2 * np.pi * k * (i - j) / N)
                    p_k = (0.5/N) * (2 * omega1 * omega2) / (omega1 + omega2) * np.cos(2 * np.pi * k * (i-j) / N)
                    r_k = (0.5j/N) * (omega2 - omega1) / (omega1 + omega2) * np.cos(2 * np.pi * k * (i-j) / N)

                x_l += x_k 
                p_l += p_k
                r_l += r_k

            x += x_l
            p += p_l
            r += r_l

            X[i-1][j-1] = x
            P[i-1][j-1] = p
            R[i-1][j-1] = r

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
    
    #---------------------------------Gamma matrix-------------

    iJG = 1j * np.matmul(J,Gamma)

    eigenvalC, eigenvecC = eigenK(iJG)
    
    eigenvalC = eigenvalC[eigenvalC >= 0]


    #-----------------(angular)mode_entropy-------------------------------

    mode_entropy=0

    for i in range(0,N_sub):
        modeS = (eigenvalC[i]+0.5) * np.log(eigenvalC[i]+0.5) - (eigenvalC[i]-0.5) * np.log(eigenvalC[i]-0.5)
         
        mode_entropy += modeS
        
        mode_entropy = np.real(mode_entropy)

    #-----------------total_entropy----------------------------------------

    S_total = (2*l+1) * mode_entropy

    S_t += S_total

    return S_t


def fit_func(n_sub, a, b, c):
    #return a * np.log((200/np.pi) * np.sin( np.pi * n_sub/200)) + b
    return a * 4 * np.pi * n_sub**2 + b * np.log(n_sub**2) + c 

def plot(N, n_ini, n_fin, step, l_max, C1, C2, m1, m2): 
    
    d = int((n_fin - n_ini)/step)
    #print(d)

    ps_array = np.array([])
    n_array  = np.array([])

    start_time = time.time()

    for i in range(0, d+1):

        ps = matrix(N, n_ini + step*i, l_max, C1, C2, m1, m2)
 
        print('Progress: ' + str(n_ini + step*i) + '/' +str(N))

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
    arealaw.to_excel('K_N200_C1_1_C2_1_m1_2_m1_2.xlsx', sheet_name='sheet1', index=False)

    plt.scatter(4* np.pi * n_array**2 * 0.001, ps_array, label='data')     
    #plt.plot(4* np.pi * n_array**2 * 0.001, fit_func(n_array, *popt), label='fit: a=%6.4f, b=%6.4f, c=%6.4f' % tuple(popt), color='r')

    plt.title('N=' + str(N) + ', C1=' + str(C1) + ', C2=' + str(C2) + ', m1=' + str(m1) + ', m2=' + str(m2))
    plt.xlabel('Area/1000')
    plt.ylabel('Pseudo Entropy')
    plt.legend()
    plt.savefig('K_minimal_N'+str(N)+'_m1_'+str(m1)+', m2_'+str(m2)+'C_'+str(C1)+'_'+str(C2))

    plt.show()


plot(10, 1, 10, 1, 30, 1, 1, 2, 2)

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




