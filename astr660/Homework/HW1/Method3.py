import numpy as np
import math
import matplotlib.pyplot as plt

def Method3(x, N):

    result=0.
    result1=[]
    n=[]
    result2 = np.exp(-x)
    Error=[]

    for i in range(0,N+1):
        result +=   (x)**i/math.factorial(i)

        n= np.append(n,i)
        # result1= np.append(result1,1/result)

        error = math.log10(abs(((1/result)-result2)/result2))
        Error=np.append(Error, error)
        
    plt.plot(n, Error, label='x='+str(x))
    plt.title("Error of different x")
    plt.xlabel("N")
    plt.ylabel("log_10(error)")

Method3(1,100)
Method3(2,100)
Method3(5,100)

Method3(10,100)
Method3(100,100)

plt.legend(loc='upper right')
plt.show()

