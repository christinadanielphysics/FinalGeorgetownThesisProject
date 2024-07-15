from Exact.Green import *
from Exact.Problem import *
from Exact.Parameters import *

from math import cos, exp, pi
from scipy.integrate import quad



U = 1
t_1 = -1
N_up = 2
N_down = 2
V = 4
chemical_potential = U/2
k = 0
operator_spin = "up"
t_0 = 0
t_2 = 0


parameters = Parameters(t_0,t_1,t_2,U,chemical_potential,V)
problem = Problem(parameters,operator_spin,N_up,N_down)
green = Green(problem)
tol = 1e-2

# moment
def mu(k,n):
    greater_angular_frequency_differences, greater_weights, lesser_angular_frequency_differences, lesser_weights = green.spectral_function_up_for_one_wavevector(k)
    greater_sum = 0
    for g,greater_weight in enumerate(greater_weights):
        greater_sum = greater_sum + greater_weight * (greater_angular_frequency_differences[g])**n
    lesser_sum = 0
    for l,lesser_weight in enumerate(lesser_weights):
        lesser_sum = lesser_sum + lesser_weight * (lesser_angular_frequency_differences[l])**n
    return (lesser_sum + greater_sum) 



# a = -10
# b = 10
# call quad to integrate mu from a to b
#res, err = quad(mu, -10, 10)

k = 3 * pi / 2# one allowed wavevector
# n = 1 # first moment
n = 2 # second moment
print("The second moment is",mu(k,n))

#print("The numerical result is", res, "+/-", err)