import numpy as np
import matplotlib.pyplot as plt

real_answer = 8/3

def square_function(x):
    y = x**2

    return y

def monte_carlo(seed1, seed2, a, b, N):

    np.random.seed(seed1)                      #
    sample_x = a + (b-a) * np.random.rand(N)   #
                                               #
    np.random.seed(seed2)
                                               #
    function_y = square_function(sample_x)     #y for the real circle
    max_y = max(function_y)
    
    sample_y = max_y * np.random.rand(N)

    count = 0

    for i in range(N):
        if (function_y[i] >= sample_y[i]):
            count += 1

    answer = (count / N) * ((b-a) * max_y)

    return answer

print(monte_carlo(10, 100, 0, 2, 1000000))

error = 0
Error = np.array([])
n = np.array([])
Error_expect = np.array([])
Error_test = np.array([])

for N in range (1,8):
    error = abs((real_answer - monte_carlo(10, 100, 0, 2, 10**N))/real_answer)
    error_expect = 1 / np.sqrt(10**N)

    n = np.append(n, N)
    Error = np.append(Error, np.log10(error)/np.log10(10**N))
    Error_expect = np.append(Error_expect, error_expect/np.log10(10**N))

    Error_test = np.append(Error_test, error_expect/10**N)

    print("The relative error:" "N="+str(10**N),  "Error="  + str(error))

plt.subplot(1,2,1)
plt.plot(n, Error, 'r', label='Monte Carlo')
plt.plot(n, Error_expect, 'b', label='N^(-1/2)', linestyle='--')
plt.legend(loc='upper right')
plt.xlabel('log_10(N)')
plt.ylabel('log_10(Error)')
plt.title('log_10(Error) v.s. log_10(N)')

plt.subplot(1,2,2)
plt.plot(10**n, Error_test, 'r', label='Monte Carlo')
plt.plot(10**n, Error_test, 'b', label='N^(-1/2)', alpha = 0.5)
plt.legend(loc='upper right')
plt.xlabel('N')
plt.ylabel('Relative error')
plt.title('Error v.s. N compared to N^(-1/2)')

plt.show()


