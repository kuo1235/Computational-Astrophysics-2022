import pandas as pd
import matplotlib.pyplot as plt 


D0 = pd.read_fwf('advection_00000.d')

D2 = pd.read_fwf('advection_00020.d')
D4 = pd.read_fwf('advection_00040.d')
D6 = pd.read_fwf('advection_00060.d')
D8 = pd.read_fwf('advection_00080.d')

D20 = pd.read_fwf('advection_00200.d')
D40 = pd.read_fwf('advection_00400.d')
D60 = pd.read_fwf('advection_00600.d')
D80 = pd.read_fwf('advection_00800.d')




plt.plot(D0['x'], D0['U(x)'], label='n=0', color='C0')

plt.plot(D2['x'], D2['U(x)'], label='n=2', color='C1', alpha=0.8)
plt.plot(D4['x'], D4['U(x)'], label='n=4', color='C2', alpha=0.6)
plt.plot(D6['x'], D6['U(x)'], label='n=6', color='C3', alpha=0.4)
plt.plot(D8['x'], D8['U(x)'], label='n=8', color='C4', alpha=0.2)


#plt.plot(D20['x'], D20['U(x)'], label='n=20', color='C1')
#plt.plot(D40['x'], D40['U(x)'], label='n=40', color='C2')
#plt.plot(D60['x'], D60['U(x)'], label='n=60', color='C3')
#plt.plot(D80['x'], D80['U(x)'], label='n=80', color='C4')


#plt.title('Upwind method, CFL=1.2')
plt.title('FTCS method')
plt.xlabel('x')
plt.ylabel('u(x)')
plt.legend()

#plt.savefig("Upwind_method")
plt.savefig("FTCS_method")

plt.show()

