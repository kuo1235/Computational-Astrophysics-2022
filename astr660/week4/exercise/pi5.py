import numpy as np

print("please enter the value N you want:")

N = 10

seed1 = 10000
seed2 = 5000


np.random.seed(seed1)
print(np.random.random(N))

np.random.seed(seed2)
print(np.random.random(N))









