import numpy as np

amin = 1
amax = 101

x = range(amin, amax)
xarray = np.array(x)

#print(xarray)
#print(len(xarray))

y = []

for i in range(len(xarray)):
		
	summation = xarray[i]
	if i > 0:
		summation += y[i-1]
	y.append(summation)
print("Summation from 1 to 100 = " + str(summation))

