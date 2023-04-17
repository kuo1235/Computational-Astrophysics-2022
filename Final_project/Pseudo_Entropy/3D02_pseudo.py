import numpy as np
import matplotlib.pyplot as plt
from scipy import linalg
from scipy.optimize import curve_fit
import time
import pandas as pd


#m1 = 1.0 * 10**(-3)
#m2 = 2.5 * 10**(-4)

def eigenK(coupling_matrix):
    eigenvalK, eigenvecK = linalg.eig(coupling_matrix)
    return eigenvalK, eigenvecK

def omega(m, N, kx, ky, kz): #k can be n dimensional 
    
    omega_square = m**2 + (2 * np.sin(np.pi * kx / N))**2 + (2 * np.sin(np.pi * ky / N))**2  + (2 * np.sin(np.pi * kz / N))**2 
    omega = omega_square**(1/2)

    return omega

def matrix(N, N_sub, m1, m2):

    #--------------------------------X, P, R matrix------------------------------------------- 
    X = np.zeros((N_sub**3, N_sub**3), dtype=complex)
    P = np.zeros((N_sub**3, N_sub**3), dtype=complex)
    R = np.zeros((N_sub**3, N_sub**3), dtype=complex)
    
    II=0
    JJ=0
    KK=0
    
    
    for i1 in range(1, N_sub+1):        
        for i2 in range(1, N_sub+1):
            for i3 in range(1, N_sub+1):
                
                II += 1
                JJ = 0
                
                #print("(i1,i2,i3)=" + str(i1)+str(i2)+str(i3) + "_II=" + str(II))
                for j1 in range(1, N_sub+1):
                    for j2 in range(1, N_sub+1):
                        for j3 in range(1, N_sub+1): 
                                               
                            JJ += 1
                                                                                                                     
                            x = 0
                            p = 0
                            r = 0

                            #print("(j1,j2,j3)=" + str(j1)+str(j2)+str(j3) + "_JJ=" + str(JJ))

                            for kx in range(0, N):
                                for ky in range(0, N):
                                    for kz in range(0,N):
                                        
                                        x_k = (0.5/N**3) * (2/(omega(m1, N, kx, ky, kz) + omega(m2, N, kx, ky, kz))) * np.cos(2 * np.pi * (kx) * (II - JJ) / N) * np.cos(2 * np.pi * (ky) * (II - JJ) / N) * np.cos(2 * np.pi * (kz) * (II - JJ) / N)
                                        p_k = (0.5/N**3) * (2*omega(m1, N, kx, ky, kz)*omega(m2, N, kx, ky, kz)) / (omega(m1, N, kx, ky, kz) + omega(m2, N, kx, ky, kz)) * np.cos(2 * np.pi * (kx) * (II-JJ) / N) * np.cos(2 * np.pi * (ky) * (II - JJ) / N) * np.cos(2 * np.pi * (kz) * (II - JJ) / N)
                                        r_k = (0.5j/N**3) * (omega(m2, N, kx, ky, kz)-omega(m1, N, kx, ky, kz)) / (omega(m1, N, kx, ky, kz) + omega(m2, N, kx, ky, kz)) * np.cos(2 * np.pi * (kx) * (II-JJ) / N) * np.cos(2 * np.pi * (ky) * (II - JJ) / N) * np.cos(2 * np.pi * (kz) * (II - JJ) / N)
    
                #x_k = (0.5/N) * (2/(omega(m1, N, k) + omega(m2, N, k))) * np.exp(2j * np.pi * k * (i-j) / N)
                #p_k = (0.5/N) * (2*omega(m1, N, k)*omega(m2, N, k)) / (omega(m1, N, k) + omega(m2, N, k)) * np.exp(2j * np.pi * k * (i-j) / N)
                #r_k = (0.5/N) * (omega(m2, N, k)-omega(m1, N, k)) / (omega(m1, N, k) + omega(m2, N, k)) * np.exp(2j * np.pi * k * (i-j) / N)

                            x += x_k 
                            p += p_k
                            r += r_k
                            
                            
                    X[II-1][JJ-1] = x
                    P[II-1][JJ-1] = p
                    R[II-1][JJ-1] = r 
                    
                    
                    #X[i1-1][j1-1][1] = x
                    #P[i1-1][j1-1][1] = p
                    #R[i1-1][j1-1][1] = r
                    
                    #X[i2-1][j2-1][1] = x
                    #P[i2-1][j2-1][1] = p
                    #R[i2-1][j2-1][1] = r
    #print(X)
    #print(P)
    #print(R)
    
    #-------------------------------Gamma matrix---------------

    h1 = np.hstack((X,R))

    R_t = np.transpose(R)

    h2 = np.hstack((R_t,P))

    Gamma = np.vstack((h1,h2))

    #--------------------------------J matrix----------------

    o = np.zeros((N_sub**3,N_sub**3))

    I = np.identity(N_sub**3)

    h1 = np.hstack((o,I))
    h2 = np.hstack((-I,o))

    J = np.vstack((h1,h2))
    
    #print(J)

    #print(Gamma)

    #---------------------------iJG matrix-----------------------

    iJG = 1j * np.matmul(J,Gamma)

    #print(iJG)

    eigenvalC, eigenvecC = eigenK(iJG)
    
    print(eigenvalC)
    
    eigenvalC = eigenvalC[eigenvalC > 0.5]

    #print(np.shape(eigenvalC))


#matrix(10,5)

#def mode_entropy(eigenvalC):
    #size = np.size(eigenvalC)

    mode_entropy=0
     
    size = np.size(eigenvalC)
    
    for i in range(0,size):
        #if abs(np.real(eigenvalC[i])-0.5) < 10**-10:
            #modeS = 0
        #else:
        modeS = (eigenvalC[i]+0.5) * np.log(eigenvalC[i]+0.5) - (eigenvalC[i]-0.5) * np.log(eigenvalC[i]-0.5)
        #modeS = (np.real(eigenvalC[i]+0.5)) * np.log(np.real(eigenvalC[i])+0.5) - (np.real(eigenvalC[i])-0.5) * np.log(np.real(eigenvalC[i])-0.5)
 
        #print(modeS) 
         
        mode_entropy += modeS
    
    print(mode_entropy)

    return mode_entropy

#matrix(10,4,1,1) 

def fit_func(n_sub, a, b):
    return a * (10/np.pi)**2 * np.sin( (np.pi * n_sub / 30))**2 + b
    #return a * n_sub**2 + b * np.log(n_sub**2) + c

def plot(N, n_ini, n_fin, step, m1, m2): 
    
    d = int((n_fin - n_ini)/step)

    ps_array = np.array([])
    n_array  = np.array([])

    start_time1 = time.time()

    for i in range(0, d+1):

        start_time2 = time.time()

        ps = matrix(N, n_ini + step*i ,m1, m2)
        
        ps = np.real(ps)
        
        print('Progress: ' + str(n_ini + step*i) + '/' +str(N))
        print('Time used: {} sec'.format(time.time()-start_time2))
        
        ps_array = np.append(ps_array, ps)
        n_array  = np.append(n_array, n_ini + step*i) 
    
    print(ps_array)
    print(n_array)
    
    popt, pcov = curve_fit(fit_func, n_array, ps_array)
    print(popt)

    print('Total Time used: {} sec'.format(time.time()-start_time1))
    
    col1 = "n_sub"
    col2 = "Pseudo_Entropy"
    col3 = "a"
    col4 = "b"
    col5 = "c"
    
    arealaw = pd.DataFrame({col1:n_array, col2:ps_array, col3:popt[0], col4:popt[1]})
    arealaw.to_excel('3D_N10_m1_1e-5_m2_1e-5.xlsx', sheet_name='sheet1', index=False)

    plt.scatter(n_array, ps_array, label='data')     
    plt.plot(n_array, fit_func(n_array, *popt), label='fit: a=%6.4f, b=%6.4f' % tuple(popt), color='r')

    plt.title('N=' + str(N) + ', m1=' + str(m1) + ', m2=' + str(m2))
    plt.xlabel('N_sub')
    plt.ylabel('Pseudo Entropy')
    plt.legend()
    plt.savefig('3D_N'+str(N)+'_m1_'+str(m1)+'_m2_'+str(m2))

    plt.show()

plot(10, 2, 4, 2, 1.0*10**(-5), 1.0*10**(-5))


