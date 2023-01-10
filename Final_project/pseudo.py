import numpy as np
import matplotlib.pyplot as plt
from scipy import linalg

m1 = 1.0 * 10**(-3)
m2 = 2.5 * 10**(-4)

def eigenK(coupling_matrix):
    eigenvalK, eigenvecK = linalg.eig(coupling_matrix)
    return eigenvalK, eigenvecK

def omega(m, N, k): #k can be n dimensional 
    
    omega_square = m**2 + (2 * np.sin(np.pi * k / N))**2
    omega = omega_square**(1/2)

    return omega

def matrix(N, n_ini, n_fin):

    #--------------------------------X, P, R matrix------------------------------------------- 
    X = np.zeros((N, N))
    P = np.zeros((N, N))
    R = np.zeros((N, N), dtype=complex)
    for i in range(0, N):        
        for j in range(0, N):
            x = 0
            p = 0
            r = 0

            for k in range(n_ini, n_fin):
                x = (0.5/N) * (2/(omega(m1, N, k) + omega(m2, N, k))) * np.cos(2 * np.pi * k * (i-j) / N)
                p = (0.5/N) * (2 * omega(m1, N, k) * omega(m2, N, k)) / (omega(m1, N, k) + omega(m2, N, k)) * np.cos(2 * np.pi * k * (i-j) / N)
                r = (0.5/N) * (omega(m2, N, k)-omega(m1, N, k)) / (omega(m1, N, k) + omega(m2, N, k)) * np.cos(2 * np.pi * k * (i-j) / N)
                
                #x = (0.5/N) * (2/(omega(m1, N, k) + omega(m2, N, k))) * np.exp(2j * np.pi * k * (i-j) / N)
                #p = (0.5/N) * (2 * omega(m1, N, k) * omega(m2, N, k)) / (omega(m1, N, k) + omega(m2, N, k)) * np.exp(2j * np.pi * k * (i-j) / N)
                #r = (0.5/N) * (omega(m2, N, k)-omega(m1, N, k)) / (omega(m1, N, k) + omega(m2, N, k)) * np.exp(2j * np.pi * k * (i-j) / N)

                x += x 
                p += p
                r += r

            X[i][j] = x
            P[i][j] = p
            R[i][j] = r*1j
    
    #print(X)
    #print(P)
    #print(R)
    
    h1 = np.hstack((X,R))

    #R_t = np.transpose(R)

    h2 = np.hstack((R,P))

    Gamma = np.vstack((h1,h2))

    #--------------------------------J matrix----------------

    o = np.zeros((N,N))

    I = np.identity(N)

    h1 = np.hstack((o,I))
    h2 = np.hstack((-I,o))

    J = np.vstack((h1,h2))
    
    #print(J)

    #print(Gamma)

    #---------------------------iJG matrix-----------------------

    iJG = 1j * np.matmul(J,Gamma)

    eigenvalC, eigenvecC = eigenK(iJG)
    
    print(eigenvalC)
    #------------------------------------------------------------

    eigenvalCC, eigenvecCC =eigenK(Gamma)

    #print(eigenvalCC)

    return eigenvalC

#matrix(10,1,5)

def mode_entropy(eigenvalC):
    size = np.size(eigenvalC)

    mode_entropy=0

    for i in range(0,size):
        if abs(np.real(eigenvalC[i])-0.5) < 10**-10:
            modeS = 0
        else:
            modeS = (eigenvalC[i]+0.5) * np.log(eigenvalC[i]+0.5) - (eigenvalC[i]-0.5) * np.log(eigenvalC[i]-0.5)

        mode_entropy += modeS
    
    print(mode_entropy)

    return mode_entropy


#mode_entropy( matrix(5) )

def total_entropy(n_ini, n_fin):

    s_total = 0

    for i in range(n_ini, n_fin+1):
  
        eigenvalC = matrix(i, n_ini, n_fin)
        s_total += mode_entropy(eigenvalC)   

    print(s_total)

total_entropy(1,5)


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




