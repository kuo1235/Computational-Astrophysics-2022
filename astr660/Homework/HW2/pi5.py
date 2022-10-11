import numpy as np
import matplotlib.pyplot as plt

def circle_function(x):
    y = np.sqrt(1. - x**2)

    return y

def monte_carlo(seed1, seed2, a, b, N):

    np.random.seed(seed1)                      #  
    sample_x = a + (b-a) * np.random.rand(N)   #
                                               #
    np.random.seed(seed2)                      #
    sample_y = np.random.rand(N)               #Construct random point (x, y)
    
    function_y = circle_function(sample_x)     #y for the real circle
       
    count = 0
    
    for i in range(N):
        if (function_y[i] >= sample_y[i]):
            count += 1
  
    answer = (count / N) * (b-a)**2
    
    return answer

print(monte_carlo(1000, 10000, -1, 1, 10000))








