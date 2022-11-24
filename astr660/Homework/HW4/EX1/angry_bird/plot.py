import pandas as pd
import matplotlib.pyplot as plt 
import numpy as np

euler = pd.read_fwf('angry_euler.dat')
RK2 = pd.read_fwf('angry_RK2.dat')
RK4 = pd.read_fwf('angry_RK4.dat')

euler1 = pd.read_fwf('angry_euler1.dat')
RK2_1 = pd.read_fwf('angry_RK2_1.dat')
RK4_1 = pd.read_fwf('angry_RK4_1.dat')

euler2 = pd.read_fwf('angry_euler2.dat')
RK2_2 = pd.read_fwf('angry_RK2_2.dat')
RK4_2 = pd.read_fwf('angry_RK4_2.dat')

euler3 = pd.read_fwf('angry_euler3.dat')
RK2_3 = pd.read_fwf('angry_RK2_3.dat')
RK4_3 = pd.read_fwf('angry_RK4_3.dat')


#plt.plot(euler1['posx'], euler1['posy'], label='Euler dt=0.1', color='C0')
#plt.plot(RK2_1['posx'], RK2_1['posy'], label='RK2 dt=0.1', color='C1')
#plt.plot(RK4_1['posx'], RK4_1['posy'], label='RK4 dt=0.1', color='C2')

#plt.plot(euler['N'], euler['posy'], label='Euler dt=1.0', color='C0')
#plt.plot(euler1['posx'], euler1['posy'], label='Euler dt=0.1', color='C1')
#plt.plot(euler2['posx'], euler2['posy'], label='Euler dt=0.01', color='C2')
#plt.plot(euler3['posx'], euler3['posy'], label='Euler dt=0.001', color='C3')

Euler_Error = np.sum(euler['Newerr_y'])
Euler1_Error = np.sum(euler1['Newerr_y'])
Euler2_Error = np.sum(euler2['Newerr_y'])
Euler3_Error = np.sum(euler3['Newerr_y'])

RK2_Error = np.sum(RK2['Newerr_y'])
RK2_1_Error = np.sum(RK2_1['Newerr_y'])
RK2_2_Error = np.sum(RK2_2['Newerr_y'])
RK2_3_Error = np.sum(RK2_3['Newerr_y'])

RK4_Error = np.sum(RK4['Newerr_y'])
RK4_1_Error = np.sum(RK4_1['Newerr_y'])
RK4_2_Error = np.sum(RK4_2['Newerr_y'])
RK4_3_Error = np.sum(RK4_3['Newerr_y'])

plt.plot([Euler_Error,Euler1_Error,Euler2_Error,Euler3_Error], [0,-1,-2,-3], label='Euler', color='C0')
plt.plot([RK2_Error,RK2_1_Error,RK2_2_Error,RK2_3_Error], [0,-1,-2,-3], label='RK2', color='C1')
plt.plot([RK4_Error,RK4_1_Error,RK4_2_Error,RK4_3_Error], [0,-1,-2,-3], label='RK4', color='C2')


plt.title('Error')
plt.xlabel('log_ε')
plt.ylabel('log_dt')
plt.legend()


plt.savefig("Error")

plt.show()

