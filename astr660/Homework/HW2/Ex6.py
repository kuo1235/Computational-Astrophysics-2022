import numpy as np
import matplotlib.pyplot as plt

# Define the equations

def FDM(x, h):
    return (np.cos(x+h) - np.cos(x))/h

def CDM(x, h):
    return (np.cos(x+h) - np.cos(x-h))/2/h

# Define relative errors
# Note that d/dx cos(x) = -sin(x)

def EPS_FDM(x, h):
    err = abs((FDM(x, h) + np.sin(x))/np.sin(x))
    return err

def EPS_CDM(x, h):
    err = abs((CDM(x, h) + np.sin(x))/np.sin(x))
    return err

# Define the range of h we want to use
N = np.logspace(1,8,30)
h = 1/N

# Visualize the results
plt.plot(np.log10(h), np.log10(EPS_FDM(0.1,h)), label='FDM $x=0.1$', linestyle='solid' , color='C0')
plt.plot(np.log10(h), np.log10(EPS_CDM(0.1,h)), label='CDM $x=0.1$', linestyle='dashed', color='C0')
plt.plot(np.log10(h), np.log10(EPS_FDM(1,h)), label='FDM $x=1$', linestyle='solid' , color='C1')
plt.plot(np.log10(h), np.log10(EPS_CDM(1,h)), label='CDM $x=1$', linestyle='dashed', color='C1')
plt.plot(np.log10(h), np.log10(EPS_FDM(100,h)), label='FDM $x=100$', linestyle='solid' , color='C2')
plt.plot(np.log10(h), np.log10(EPS_CDM(100,h)), label='CDM $x=100$', linestyle='dashed', color='C2')
plt.xlim(np.log10(h[0]), np.log10(h[-1]))
plt.xlabel('$\log_{10}h$')
plt.ylabel('$\log_{10}\epsilon$')
plt.legend()
plt.tight_layout()
plt.savefig('Derivative.png', dpi=300)
