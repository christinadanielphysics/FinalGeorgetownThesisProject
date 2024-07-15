from plotting_exact import exact as exact_real_w
import matplotlib.pyplot as plt
from matsubara_exact import exact as exact_matsubara
from plotting_approximate import go as approximate

plt.figure(1)

U_values = [1,2,4,8]
t = 1

n = 500
sampling_period = 0.1

for U in U_values:
    exact_real_w(U, t)
    exact_matsubara(U, t)
    approximate(U,t,n,sampling_period)


plt.ylabel(r"$E_0$")
plt.xlabel(r"$U$")
plt.show()