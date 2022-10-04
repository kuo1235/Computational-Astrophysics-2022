import numpy as np
import math
import matplotlib.pyplot as plt

result=0.
x=100.
N=100
result1=[]
n=[]
result2 = np.exp(-100)
Error=[]

for i in range(0,N+1):
    result += (-x)**i/math.factorial(i)

    n= np.append(n,i)
    result1= np.append(result1, result)
    
    plt.subplot(1,2,1)
    plt.plot(n, result1)
    plt.title("exp^(-x) for x=1")
    plt.xlabel("N")
    plt.ylabel("exp^(-1)")

    error = math.log10(abs((result-result2)/result2))
    Error=np.append(Error, error)

    plt.subplot(1,2,2)
    plt.plot(n, Error)
    plt.title("Error")        
    plt.xlabel("N")
    plt.ylabel("log_10(error)")

#print(type(result))
plt.suptitle("Method1")
plt.show()







