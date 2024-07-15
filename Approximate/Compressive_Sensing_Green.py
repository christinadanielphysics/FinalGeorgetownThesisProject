from sympy import symbols,lambdify
import matplotlib.pyplot as plt
from numpy import linspace,zeros
from Exact.tolerances import *
from cmath import pi

class Compressive_Sensing_Green:
    def __init__(self,problem,greater_weights_for_all_k_n,lesser_weights_for_all_k_n,w_values_greater,w_values_lesser):
        self.greater_weights_for_all_k_n = greater_weights_for_all_k_n
        self.lesser_weights_for_all_k_n = lesser_weights_for_all_k_n
        self.w_values_greater = w_values_greater # w_g - w_0 = -W_g
        self.w_values_lesser = w_values_lesser # w_0 - w_l = -W_l
        self.allowed_wavevectors = problem.allowed_wavevectors
        self.computational_zero = TOLERANCE_computational_zero_for_denominator_approximate # for weights close to zero
    def real_part_in_k_space(self,integer_n):
        w = symbols('w',real=True)
        Re_G = 0
        for g,weight in enumerate(self.greater_weights_for_all_k_n[integer_n]):
            if abs(weight) >= self.computational_zero:
                w_g_values = self.w_values_greater[integer_n]
                Re_G = Re_G - weight / (w - w_g_values[g])
        for l,weight in enumerate(self.lesser_weights_for_all_k_n[integer_n]):
            if abs(weight) >= self.computational_zero:
                w_l_values = self.w_values_lesser[integer_n]
                Re_G = Re_G + weight / (w - w_l_values[l])
        return Re_G 
    def plot_real_part_in_k_space(self):

        w = symbols('w',real=True)
        w_values = linspace(-5,5,1000)
        for integer_n,k_n in enumerate(self.allowed_wavevectors):
            Re_G = self.real_part_in_k_space(integer_n)
            Re_G_lamdified = lambdify(w, Re_G, modules=['numpy'])
            Re_G_values = Re_G_lamdified(w_values)
            if Re_G == 0:
                print("Re_G is zero!")
                Re_G_values = zeros(len(w_values))

            my_linestyle = "solid"
            if integer_n == 4 or integer_n == 5:
                my_linestyle = "dashed"
            
            plt.figure(3)
            plt.plot( w_values, Re_G_values, label="$n=$"+str(integer_n), linestyle = my_linestyle)
            plt.xlabel(r"$\omega$")
            plt.ylabel(r"$\mathfrak{Re} [ \lim_{T \rightarrow 0^+} G_{k_n k_n \sigma}(\omega) ]$")
            plt.title("from compressive sensing")
            plt.legend()
            plt.savefig("Approximate/Figures/green_real_part_in_k_space.png",dpi=800)
        plt.show()
        
    


          
    

