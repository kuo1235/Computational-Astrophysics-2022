import pandas as pd
import matplotlib.pyplot as plt 


D1 = pd.read_fwf('binary_001.dat')
D2 = pd.read_fwf('binary_002.dat')



plt.plot(D1['posx [code]'], D1['posy [code]'], label='star1', color='C0')
plt.plot(D2['posx [code]'], D2['posy [code]'], label='star2', color='C1')

plt.title('star1&star2 trajectories, dt=0.01year, RK4')
plt.xlabel('x')
plt.ylabel('y')
plt.legend()

#plt.savefig("Upwind_method")
plt.savefig("binary_trajectory")

plt.show()

