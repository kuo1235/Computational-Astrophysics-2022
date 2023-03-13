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
    couplingK = np.zeros((N, N))
    for i in range(1, N):
        couplingK[i-1][i-1] = 2 + ( 1 / i**2 )*(0.5 + L*(L+1))
        couplingK[i-1][i] = -( (i + 0.5)**2  / (i*(i+1)) )
        couplingK[i][i-1] = couplingK[i-1][i]
   
    couplingK[0][0] = 9/4 + L*(L+1) 
    couplingK[N-1][N-1] =  2 + ( 1 / (N)**2 )*(0.5 + L*(L+1))
   
    return couplingK
#omega(10, 10 ,1, 1)

def S_L(N, N_sub, L, C1, C2, m1, m2):
        
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
                       
            #print("i=" + str(i) + ", j=" + str(j))

            #S_t = 0

            #for L in range(0, L_max+1):
                
            #for k in range(-N+1, N):    
            for k in range(0, N):
        
                eigenvalK1, eigenvecK1 = eigenK(coupling_matrix(N , L, C1, m1))
                eigenvalK2, eigenvecK2 = eigenK(coupling_matrix(N , L, C2, m2))
                    
                #eigenvalK1 = sorted(eigenvalK1)
                #eigenvalK2 = sorted(eigenvalK2)

                #omega1 = eigenvalK1[k]
                #omega2 = eigenvalK2[k]

                omega1 = np.sqrt(eigenvalK1[k])
                omega2 = np.sqrt(eigenvalK2[k])

                #print('omega1 = ' + str(omega1) )
                #print('omega2 = ' + str(omega2) )
                
                x_k = (0.5/N) * (2 / (omega1 + omega2)) * np.cos(2 * np.pi * k * (i - j) / N)
                p_k = (0.5/N) * (2 * omega1 * omega2) / (omega1 + omega2) * np.cos(2 * np.pi * k * (i-j) / N)
                r_k = (0.5j/N) * (omega2 - omega1) / (omega1 + omega2) * np.cos(2 * np.pi * k * (i-j) / N)

                x += x_k
                p += p_k
                r += r_k
            
            X[i-1][j-1] = x
            P[i-1][j-1] = p
            R[i-1][j-1] = r
                
    return X, P, R

def submatrix(matrix, row_e, col_e):

    subrow = range(0, row_e)
    subcol = range(0, col_e)

    result = matrix[np.ix_(subrow,subcol)]

    return result


def mode_entropy(eigenvalC):
    size = np.size(eigenvalC)

    mode_entropy=0

    for i in range(0,size):
        if abs(np.real(eigenvalC[i])-0.5) < 10**-10:
            modeS = 0
        else: 
            modeS = (eigenvalC[i]+0.5) * np.log(eigenvalC[i]+0.5) - (eigenvalC[i]-0.5) * np.log(eigenvalC[i]-0.5)

        mode_entropy += modeS

    return mode_entropy

def fit_func(x,a,b,c):
    return a * 4*np.pi * x**2 + b * np.log(x**2) + c
    #return a * 4*np.pi * x**2

def total_entropy(l_max, N, n_ini, n_fin, C1, C2, m1, m2):
    
    R_list = np.array([])
    log_R_list = np.array([])
    area_list = np.array([])
    S_list = np.array([])  
    
    #R_list = []
    #log_R_list = []
    #area_list = []
    #S_list = []
    
    start_time = time.time()

    for r in range (n_ini, n_fin+1):
        S_t = 0
        for l in range (0, l_max+1):
            
            Xmatrix, Pmatrix, Rmatrix = S_L(N, n_ini + r, l, C1, C2, m1, m2)

            Xreduce = submatrix(Xmatrix, r , r)
            Preduce = submatrix(Pmatrix, r , r)

            eigenvalC, eigenvecC = eigenK(np.matmul(Xreduce,Preduce))
            #XPhalf = np.matmul(eigenvecC, np.matmul(np.diag((eigenvalC**(1/2))), linalg.inv(eigenvecC)))
            
            #eigenvalC, eigenvecC = eigenK(XPhalf)   

            eigenvalC = (eigenvalC)**(1/2)

            #print("eigenvalueC is:" + str(eigenvalC))

            S_mode = mode_entropy(eigenvalC)
            #print("mode_entropy is:" + str(S_mode))

            S_total = (2*l+1) * np.real(S_mode)

            S_t += S_total
            
        R_list = np.append(R_list, r+0.5) 
        log_R_list = np.append(log_R_list, np.log((r+0.5)**2))
        area_list = np.append(area_list, 4 * np.pi * (r+0.5)**2)    
        S_list = np.append(S_list, S_t) 
        
        #R_list.append(r+0.5)
        #log_R_list.append(np.log(r+0.5))
        #area_list.append(4 * np.pi * (r+0.5)**2)
        #S_list.append(S_t)

        print("r=" + str(r) +", S_total_r="+str(S_t))    
        #print("eigenvalueC is:" + str(eigenvalC))

        #print("{} | ".format(r), end="") 

    #return S_t
    print('Time used: {} sec'.format(time.time()-start_time))    
    #path = 'Massless_Scalar_Field.txt'
    #f = open(path, 'w')
    #print(list1, file=f)
    #print(list2, file=f)
    #f.close()

    #dataset = pd.read_csv("path")
    #x = dataset.iloc[:, 1].values.reshape(-1,1)
    #y = dataset.iloc[:, 2].values

    #X = area_list.reshape(-1,1)
    #y = S_list

    popt, pcov = curve_fit(fit_func, R_list, S_list)
    print(popt)
    
    col1 = "R"
    col2 = "S_total_r"
    col3 = "a"
    col4 = "b"
    col5 = "c"

    #arealaw = pd.DataFrame({col1:R_list, col2:S_list, col3:popt[0], col4:popt[1], col5:popt[2]})
    #arealaw.to_excel('N_n_L.xlsx', sheet_name='sheet1', index=False)


    plt.scatter(R_list, S_list, label='data')     
    plt.plot(R_list, fit_func(R_list, *popt), label='fit: a=%6.4f, b=%6.4f, c=%6.4f' % tuple(popt), color='r')

    plt.title('N=' + str(N) + ', ' + str(n_ini) + '<n<' + str(n_fin) + ', l_max=' + str(l_max))
    plt.xlabel('Area')
    plt.ylabel('S')
    plt.legend()
    #plt.savefig('N'+str(N)+'_'+str(n_ini)+'n'+str(n_fin)+'_L'+str(l_max))

    plt.show()

areaLaw = total_entropy(20, 20, 1, 10, 1, 1, 1, 1)  

