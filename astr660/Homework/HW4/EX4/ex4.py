import numpy as np
import matplotlib.pyplot as plt

h = 0.005 #step size

a = 1.0 #upper limit
b = 0.0 #lower limit

N = int((a-b)/h)

x_anl = np.array([])
t = np.array([])
t_list = np.linspace(0,1,N+1)
y_list = np.zeros([N+1])

y_list[0] = 1
y_list[N] = 1

def fdm_matrix(row_num):
    couplingK = np.zeros((row_num, row_num))
    for i in range(0, row_num-1):
        couplingK[i][i] = -2
        couplingK[i][i+1] = 1
        couplingK[i+1][i] = couplingK[i][i+1]
    couplingK[row_num-1][row_num-1] = -2    
    return couplingK

M = fdm_matrix(N)
#print(M)

for i in range(1,N+1):
    t = np.append(t, 6 * i * h**3)        

t[0] = t[0]-1
t[N-1] = t[N-1]-1

#print(t)

y = np.linalg.solve(M, t)

y = np.insert(y, [0], 1)

#print(y)


for i in range(0,N+1):
    x_anl = np.append(x_anl, t_list[i]**3 - t_list[i] + 1)

plt.xlabel("t")
plt.ylabel("y")
plt.title("Finite difference method, h=0.05")

plt.scatter(t_list, y, label="data")
plt.plot(t_list, x_anl, color="r", label="analytic solution")

plt.legend()
plt.savefig("EX4_2")

plt.show()
