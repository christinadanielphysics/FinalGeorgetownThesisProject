from cmath import exp
import matplotlib.pyplot as plt
from numpy import zeros

class TimeEvolution:
    def __init__(self,problem,t_values):
        self.problem = problem
        self.minus_dimension = problem.hubbard_minus.dimension
        self.plus_dimension = problem.hubbard_plus.dimension
        self.green = problem.green
        self.plus_eigenvalues = problem.hubbard_plus.eigenvalues
        self.minus_eigenvalues = problem.hubbard_minus.eigenvalues
        self.ground_state_energy = problem.hubbard_N.ground_state_energy
        self.t_values = t_values
        self.V = problem.parameters.V
    def lesser_up(self,i,j,t):
        w_0 = self.ground_state_energy
        sum = 0
        for l in range(self.minus_dimension):
            w_l = self.minus_eigenvalues[l]
            sum = sum + 1j * exp(1j * (w_l - w_0 ) * t) * self.green.lesser_weight_up(i,j,l) 
        return sum
    def greater_up(self,i,j,t):
        w_0 = self.ground_state_energy
        sum = 0
        for g in range(self.plus_dimension):
            w_g = self.plus_eigenvalues[g]
            sum = sum + ( -1j ) * exp(1j * (w_0 - w_g) * t) * self.green.greater_weight_up(i,j,g) 
        return sum
    def time_evolve_lesser_up(self,i,j):
        L = len(self.t_values)
        g_real = zeros(L)
        g_imag = zeros(L)
        for index,t in enumerate(self.t_values):
            g_real[index] = self.lesser_up(i,j,t).real
            g_imag[index] = self.lesser_up(i,j,t).imag
        return g_real, g_imag
    def time_evolve_greater_up(self,i,j):
        L = len(self.t_values)
        g_real = zeros(L)
        g_imag = zeros(L)
        for index,t in enumerate(self.t_values):
            g_real[index] = self.greater_up(i,j,t).real 
            g_imag[index] = self.greater_up(i,j,t).imag
        return g_real, g_imag
    def plot_lesser_up(self,i,j):
        g_real, g_imag = self.time_evolve_lesser_up(i,j)
        plt.figure(1)
        plt.plot(self.t_values,g_real,label="real part",color="darkred")
        plt.plot(self.t_values,g_imag,label="imaginary part",color="coral")
        plt.title("Lesser")
        plt.legend()
        plt.xlabel("$t$")
        plt.savefig("Approximate/Figures/lesser_"+str(i)+"_"+str(j)+".png",dpi=800)
        plt.show()
    def plot_greater_up(self,i,j):
        g_real, g_imag = self.time_evolve_greater_up(i,j)
        plt.figure(2)
        plt.plot(self.t_values,g_real,label="real part",color="darkblue")
        plt.plot(self.t_values,g_imag,label="imaginary part",color="lightblue")
        plt.title("Greater")
        plt.legend()
        plt.xlabel("$t$")
        plt.savefig("Approximate/Figures/greater_"+str(i)+"_"+str(j)+".png",dpi=800)
        plt.show()
    def constant_greater(self,k_n,g,t):
        w_g = self.plus_eigenvalues[g]
        w_0 = self.ground_state_energy
        weight = 0
        for i in range(self.V):
            for j in range(self.V):
                weight = weight + (1/self.V) * exp(1j * k_n * (i-j)) * self.green.greater_weight_up(i,j,g) * exp(1j * (w_0 - w_g) * t) * (-1j)
        return weight
    def constant_lesser(self,k_n,l,t):
        w_l = self.minus_eigenvalues[l]
        w_0 = self.ground_state_energy
        weight = 0
        for i in range(self.V):
            for j in range(self.V):
                weight = weight + (1/self.V) * exp(1j * k_n * (i-j)) * self.green.lesser_weight_up(i,j,l) * exp(1j * (w_l - w_0) * t) * (1j)
        return weight
    def sum_greater(self,k_n,t):
        sum = 0
        for g in range(self.plus_dimension):
            sum = sum + self.constant_greater(k_n,g,t)
        return sum
    def sum_lesser(self,k_n,t):
        sum = 0
        for l in range(self.minus_dimension):
            sum = sum + self.constant_lesser(k_n,l,t)
        return sum
    def time_evolve_lesser_in_k_space(self,k_n):
        L = len(self.t_values)
        g_real = zeros(L)
        g_imag = zeros(L)
        for index,t in enumerate(self.t_values):
            real_value = self.sum_lesser(k_n,t).real
            imaginary_value = self.sum_lesser(k_n,t).imag
            if abs(real_value) > 1e-16:
                g_real[index] = real_value
            if abs(imaginary_value) > 1e-16:
                g_imag[index] = imaginary_value
        return g_real, g_imag
    def time_evolve_greater_in_k_space(self,k_n):
        L = len(self.t_values)
        g_real = zeros(L)
        g_imag = zeros(L)
        for index,t in enumerate(self.t_values):
            real_value = self.sum_greater(k_n,t).real
            imaginary_value = self.sum_greater(k_n,t).imag
            if abs(real_value) > 1e-16:
                g_real[index] = real_value
            if abs(imaginary_value) > 1e-16:
                g_imag[index] = imaginary_value
        return g_real, g_imag
    def plot_lesser_up_in_k_space(self,k_n):
        g_real, g_imag = self.time_evolve_lesser_in_k_space(k_n)
        plt.figure(1)
        plt.plot(self.t_values,g_real,label="real part",color="darkred")
        plt.plot(self.t_values,g_imag,label="imaginary part",color="coral")
        plt.title("Lesser")
        plt.legend()
        plt.xlabel("$t$")
        plt.savefig("Approximate/Figures/lesser_"+str(k_n)+".png",dpi=800)
        plt.show()
    def plot_greater_up_in_k_space(self,k_n):
        g_real, g_imag = self.time_evolve_greater_in_k_space(k_n)
        plt.figure(1)
        plt.plot(self.t_values,g_real,label="real part",color="darkred")
        plt.plot(self.t_values,g_imag,label="imaginary part",color="coral")
        plt.title("Greater")
        plt.legend()
        plt.xlabel("$t$")
        plt.savefig("Approximate/Figures/lesser_"+str(k_n)+".png",dpi=800)
        plt.show()


