import numpy as np
import matplotlib.pyplot as plt

def forward(fx,h):                         
   
    dfdx = np.zeros(len(fx))
    for i,f in enumerate(fx[:-1]):           #enumerate(sequence, [start=0])
        dfdx[i] = (fx[i+1]-fx[i])/h
    return dfdx

def center(fx,h):
    
    dfdx = np.zeros(len(fx))
    for i in range(1,len(fx)-1):
        dfdx[i] = (fx[i+1]-fx[i-1])/(2.*h)
    return dfdx

func  = lambda x: np.cos(x)         # Differential target 
dfunc = lambda x: -1 * np.sin(x)    # True answer

Error_forward = np.array([])
Error_center = np.array([])
H = np.array([])

for i in range (1,10):

    h  = 10**(-1 * i)
    x = np.linspace(0,1,(int(1/h)+1))
    y = func(x) 
    dx = x[2]-x[1]                      # 

    dydx = dfunc(x)
    dydx_forward  = forward(func(x),h)
    dydx_center   = center(func(x),h)

    error_forward = abs((dydx-dydx_forward)/dydx)
    error_center = abs((dydx-dydx_center)/dydx)

    Error_forward = np.append(Error_forward, np.log10(error_forward)/np.log10(h))
    Error_center = np.append(Error_center, np.log10(error_center)/np.log10(h))
    H = np.append(H,-i)

#print("dydx_forward at x=0.1 is:" + str(dydx_forward[1]), "dydx_forward at x=1.0 is:" + str(dydx_forward[-1]))
#print("dydx_center at x=0.1 is:" + str(dydx_center[1]), "dydx_center at x=1.0 is:" + str(dydx_center[-1]))

#plt.plot(x,dydx,'k-',label="True answer")
plt.plot(H,Error_forward,label="Forward Error")
plt.plot(H,Error_center,label="Centered Error")
plt.xlabel("log10_h")
plt.ylabel("log10_Error")
plt.legend(loc="best")

plt.show()
