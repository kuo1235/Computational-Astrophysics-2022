import numpy as np
import matplotlib.pyplot as plt
from scipy import linalg
from resource import getrusage, RUSAGE_SELF
import time
import pandas as pd
from multiprocessing import Process, Queue

def eigenK(coupling_matrix):
    eigenvalK, eigenvecK = linalg.eig(coupling_matrix)
    return eigenvalK, eigenvecK

def omega(t, m, N, kx, ky): #k can be n dimensional 
    
    omega_square = np.exp(t) * m**2 + (2 * np.sin(np.pi * kx / N))**2 + (2 * np.sin(np.pi * ky / N))**2 
    omega = omega_square**(1/2)

    return omega

def X_matrix(N, N_sub, t1, t2, m1, m2,q1):
    
    start_timeX = time.time()
    
    X = np.zeros((N_sub**2, N_sub**2), dtype=complex)
        
    II=0
    JJ=0
   
    for i1 in range(1, N_sub+1):        
        for i2 in range(1, N_sub+1):
                
            II += 1
            JJ = 0
                
            for j1 in range(1, N_sub+1):
                for j2 in range(1, N_sub+1):
                        
                    JJ += 1
                        
                    x = 0
                              
                    for kx in range(0, N):
                        for ky in range(0, N):
                    
                            x_k = (0.5/N**2) * (2/(omega(t1, m1, N, kx, ky) + omega(t2, m2, N, kx, ky))) * np.cos(2 * np.pi * (kx) * (II - JJ) / N) * np.cos(2 * np.pi * (ky) * (II - JJ) / N)
                            

                            x += x_k 
                         
                    X[II-1][JJ-1] = x
    
                    
    #print('XTime used: {} sec'.format(time.time()-start_timeX))
    #q.put(X)
    return q1.put(X) 
 
def P_matrix(N, N_sub, t1, t2, m1, m2, q2):
    
    start_timeP = time.time()
    
    P = np.zeros((N_sub**2, N_sub**2), dtype=complex)
 
    II=0
    JJ=0
   
    for i1 in range(1, N_sub+1):        
        for i2 in range(1, N_sub+1):
                
            II += 1
            JJ = 0
                
            for j1 in range(1, N_sub+1):
                for j2 in range(1, N_sub+1):
                        
                    JJ += 1
                                           
                    p = 0
          
                    for kx in range(0, N):
                        for ky in range(0, N):
                                               
                            p_k = (0.5/N**2) * (2*omega(t1, m1, N, kx, ky)*omega(t2, m2, N, kx, ky)) / (omega(t1, m1, N, kx, ky) + omega(t2, m2, N, kx, ky)) * np.cos(2 * np.pi * (kx) * (II-JJ) / N) * np.cos(2 * np.pi * (ky) * (II - JJ) / N)
                                           
                            p += p_k
                                               
                    P[II-1][JJ-1] = p

    #print('PTime used: {} sec'.format(time.time()-start_timeP))
    #print(P)  
    #q.put(P)
    return q2.put(P) 
                             
def R_matrix(N, N_sub, t1, t2, m1, m2, q3):
    
    start_timeR = time.time()
    
    R = np.zeros((N_sub**2, N_sub**2), dtype=complex)
    
    II=0
    JJ=0
   
    for i1 in range(1, N_sub+1):        
        for i2 in range(1, N_sub+1):
                
            II += 1
            JJ = 0
                
            for j1 in range(1, N_sub+1):
                for j2 in range(1, N_sub+1):
                        
                    JJ += 1
                                   
                    r = 0

                    #print("i=" + str(II) + ", j=" + str(JJ))
            
                    for kx in range(0, N):
                        for ky in range(0, N):
                    
                            r_k = (0.5j/N**2) * (omega(t2, m2, N, kx, ky)-omega(t1, m1, N, kx, ky)) / (omega(t1, m1, N, kx, ky) + omega(t2, m2, N, kx, ky)) * np.cos(2 * np.pi * (kx) * (II-JJ) / N) * np.cos(2 * np.pi * (ky) * (II - JJ) / N)
                
                            r += r_k
                    R[II-1][JJ-1] = r
    #print('RTime used: {} sec'.format(time.time()-start_timeR))
    #q.put(R)
    return q3.put(R) 

def Entropy(N, N_sub, t1, t2, m1, m2):

    #--------------------------------X, P, R matrix------------------------------------------- 
    #X = X_matrix(N, N_sub, m1, m2)
    #P = P_matrix(N, N_sub, m1, m2)
    #R = R_matrix(N, N_sub, m1, m2)
    
    #---------------------------------Multiprocessing-------------------------
    if __name__ == '__main__':
    
        #q1 = []

        q1 = Queue()
        q2 = Queue()   
        q3 = Queue()

        pX = Process(target=X_matrix, args=(N, N_sub, t1, t2, m1, m2,q1))
        pP = Process(target=P_matrix, args=(N, N_sub, t1, t2, m1, m2,q2))
        pR = Process(target=R_matrix, args=(N, N_sub, t1, t2, m1, m2,q3))
        
        #print(q1.full())


        pX.start()
        pP.start()
        pR.start()

        #print(000)

        #pX.join()
        #pP.join()
        #pR.join()
        
        #?????????????????
        #print(123)

        X = q1.get()
        P = q2.get()

        R = q3.get()
        
        #print(456)
        #print(X)
        #print(R)

        #-------------------------------Gamma matrix---------------
    
        h1 = np.hstack((X,R))
    
        R_t = np.transpose(R)
    
        h2 = np.hstack((R_t,P))
    
        Gamma = np.vstack((h1,h2))
    
        #--------------------------------J matrix----------------
    
        o = np.zeros((N_sub**2,N_sub**2))
    
        I = np.identity(N_sub**2)
    
        h1 = np.hstack((o,I))
        h2 = np.hstack((-I,o))
    
        J = np.vstack((h1,h2))
        
        #print(J)
    
        #print(Gamma)
    
        #---------------------------iJG matrix-----------------------
    
        iJG = 1j * np.matmul(J,Gamma)
    
        #print(iJG)
    
        eigenvalC, eigenvecC = eigenK(iJG)
        
        #print(eigenvalC)
        
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
     
            print(modeS) 
             
            mode_entropy += modeS
            #mode_entropy = np.real(mode_entropy)

        print(mode_entropy)
 
        return mode_entropy

#Entropy(20,9,1e-5,1e-5) 

#def fit_func(n_sub, a, b):
    #return a * (20/np.pi) * np.sin( (np.pi * n_sub / 20)) + b
    #return a * n_sub**2 + b * np.log(n_sub**2) + c

def plot(N, n_ini, n_fin, step, C1, C2, m1, m2): 

    print('N=' + str(N) + '_C1=' + str(C1) + '_C2=' + str(C2) + '_m1=' + str(m1) + '_m2=' + str(m2))

    d = int((n_fin - n_ini)/step)

    ps_iarray = np.array([])

    ps_array = np.array([])
    n_array  = np.array([])

    start_time1 = time.time()

    for i in range(0, d+1):

        start_time2 = time.time()

        ps = Entropy(N, n_ini + step*i, C1, C2, m1, m2)
        
        ps_iarray = np.append(ps_iarray, ps)

        ps = np.real(ps)
        
        print('Progress: ' + str(n_ini + step*i) + '/' +str(N))
        print('Time used: {} sec'.format(time.time()-start_time2))
        
        ps_array = np.append(ps_array, ps)
        n_array  = np.append(n_array, n_ini + step*i) 
    
        print("peak memory:", getrusage(RUSAGE_SELF).ru_maxrss / 1000 / 1000, "MB")

    print(ps_array)
    print(n_array)
    
    #popt, pcov = curve_fit(fit_func, n_array, ps_array)
    #print(popt)

    print('Total Time used: {} sec'.format(time.time()-start_time1))
    
    #print("peak memory:", getrusage(RUSAGE_SELF).ru_maxrss / 1000 / 1000, "MB")

    col2 = "n_sub"
    col3 = "Pseudo_Entropy"
    col1 = "i_PS"
    #col4 = "b"
    #col5 = "c"
  
    #arealaw = pd.DataFrame({col1:ps_iarray, col2:n_array, col3:ps_array})
    #arealaw.to_excel('N20_C1_1_C2_10000+10000i_m1_1e-5_m2_1e-5.xlsx', sheet_name='sheet1', index=False)

    plt.scatter(n_array, ps_array, label='data')     
    #plt.plot(n_array, fit_func(n_array, *popt), label='fit: a=%6.4f, b=%6.4f' % tuple(popt), color='r')

    plt.title('2D_N=' + str(N) + ', C1=' + str(C1) + ', C2=' + str(C2) + ', m1=' + str(m1) + ', m2=' + str(m2))
    plt.xlabel('N_sub')
    plt.ylabel('Pseudo Entropy')
    plt.legend()
    #plt.savefig('2D_N'+str(N) + ', C1=' + str(C1) + ', C2=' + str(C2) + '_m1_'+str(m1)+'_m2_'+str(m2))

    plt.show()

plot(20, 1, 10, 1, 10-1j, 10+1j, 1.0*10**(-5), 1.0*10**(-5))




