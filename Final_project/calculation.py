import numpy as np

for i in range(0,51):
    S = 0.0235 * 4 * np.pi * (i+0.5)**2 - 0.0055 * np.log((i+0.5)**2) - 0.0354
    print(S)

