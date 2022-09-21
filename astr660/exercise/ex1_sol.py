'''
This code can calcuate the sum of integers between amin and amax
Author: Lin Yen-Hsing (NTHU) 2022.09.15

Run this code with:
$ python ex1.py amin amax
e.g.
$ python ex1.py 1 100 
then you will get:
Summation from 1 to 100 = 5050
'''

from sys import argv
import numpy as np

# Define amin and amax
AMIN = int(argv[1])
AMAX = int(argv[2])

# An array of integers from amin to amax
num = np.arange(AMIN, AMAX+1, 1)

# Sum up the integers
SUM = 0
for i in num:
    SUM += i

# Print our the results
print(f"Summation from {AMIN} to {AMAX} = {SUM}")

