import numpy as np
import matplotlib.pyplot as plt
from scipy import linalg
from scipy.optimize import curve_fit


#m1 = 1.0 * 10**(-3)
#m2 = 2.5 * 10**(-4)



def eigenK(coupling_matrix):
    eigenvalK, eigenvecK = linalg.eig(coupling_matrix)
    return eigenvalK, eigenvecK

def omega(C, m, c, R, k): #k can be n dimensional 
    
    omega_square = k**2 + C*(m**2 + c * R)
    omega = omega_square**(1/2)

    return omega

def matrix(N, N_sub, C1, C2, m1, m2, xi, r):  #do not sure why we cam mot use "R" (generate array???)

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
                x_k = (0.5/N)*(2/(omega(C1, m1, xi, r, k) + omega(C2, m2, xi, r, k))) * np.cos(2 * np.pi * k * (i - j) / N)
                p_k = (0.5/N)*(2*omega(C1,m1,xi,r,k)*omega(C2,m2,xi,r,k))/(omega(C1,m1,xi,r,k)+omega(C2,m2,xi,r,k))*np.cos(2 * np.pi * k * (i-j) / N)
                r_k = (0.5j/N)*(omega(C2,m2,xi,r,k)-omega(C1,m1,xi,r,k))/(omega(C1,m1,xi,r,k)+omega(C2,m2,xi,r,k))*np.cos(2 * np.pi * k * (i-j) / N)
                

                #print(x_k)
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

#matrix(200,40, 1, 10, 5.0*10**(-8), 1.0*10**(-8), 1/6, 100) 

def fit_func(n_sub, a, b):
    return a * np.log((200/np.pi) * np.sin( np.pi * n_sub/200)) + b


def plot(N, n_ini, n_fin, step, C1, C2, m1, m2, xi, R): 
    
    d = int((n_fin - n_ini)/step)
    #print(d)

    ps_array = np.array([])
    n_array  = np.array([])

    for i in range(0, d+1):

        ps = matrix(N, n_ini + step*i, C1, C2, m1, m2, xi, R)

        ps_array = np.append(ps_array, ps)
        n_array  = np.append(n_array, n_ini + step*i) 

    popt, pcov = curve_fit(fit_func, n_array, ps_array)
    print(popt)

    plt.scatter(n_array, ps_array, label='data, fit: a=%6.4f, b=%6.4f' % tuple(popt))     
    #plt.plot(R_list, fit_func(R_list, *popt), label='fit: a=%6.4f, b=%6.4f, c=%6.4f' % tuple(popt), color='r')

    plt.title('N=' + str(N) + ', C1=' + str(C1) + ', C2=' + str(C2) + ', m1=' + str(m1) + ', m2=' + str(m2))
    plt.xlabel('N_sub')
    plt.ylabel('Pseudo Entropy')
    plt.legend()
    #plt.savefig('QFT_conformal_N'+str(N)+'_'+str(n_ini)+'n'+str(n_fin))

    plt.show()


plot(200, 20, 180, 2, 1, 10, 1, 1, 1/6, 100)

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




