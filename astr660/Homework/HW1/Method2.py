import numpy as np
import math
import matplotlib.pyplot as plt

element=0
SUM=0
result=0
x=1.
N=100
result1=[]
n=[]
result2 = np.exp(-1)
Error=[]
y=[]

for i in range(0,N+1):
    if i==0:
        element=1 #第0個數為1

        y= np.append(y,element) #A0, A1, A2... 的數列
        n= np.append(n,i)

        SUM=np.sum(y) #A0, A0+A1, A0+A1+A2... 
        
        error = math.log10(abs((SUM-result2)/result2))
        Error=np.append(Error, error)

       # plt.plot(n, Error, 'b', label='Method2')

    else:
        element = (-x/i)*element
       
        
        SUM = element+ np.sum(y)
        y= np.append(y,element) #A0, A1 ,A2... 的數列

        n= np.append(n,i)

        error = math.log10(abs((SUM-result2)/result2))
        Error=np.append(Error, error)

       # plt.plot(n, Error, 'b', label='Method2')

Result=0
Result1=[]
nn=[]
Result2 = np.exp(-1)
EError=[]

for j in range(0,N+1):
    Result += (-x)**j/math.factorial(j)

    nn= np.append(nn,j)
    Result1= np.append(Result1, Result)

    eerror = math.log10(abs((Result-Result2)/Result2))
    EError=np.append(EError, eerror)

   # plt.plot(nn, EError, 'r', label='Method1')

print(Error)
print(EError)

plt.plot(nn, EError, 'r', label='Method1', alpha = 1.0)
plt.plot(n, Error, 'b', label='Method2', alpha = 0.5)
plt.title("Error of Method1 and Method2")
plt.xlabel("N")
plt.ylabel("log_10(error)")
plt.legend(loc='lower left')
plt.show()



