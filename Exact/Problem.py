from .SelfEnergy import *
from .Green import *
from .FermionicOperations.Hubbard import *
from .FermionicOperations.Storage import *
from .FermionicOperations.combine import *
from .Root_Finding.Polynomial_in_w import *

import matplotlib.pyplot as plt
from sympy.plotting import plot as sympy_plot
from numpy import linspace
from sympy import lambdify,symbols,diff

class Problem:
    def __init__(self,parameters,operator_spin,N_up,N_down):
        
        self.operator_spin = operator_spin
        self.N_up = N_up
        self.N_down = N_down

        self.parameters = parameters

        self.hubbard_N = Hubbard(N_up,N_down,parameters)
        self.hubbard_N.form_hamiltonian_matrix()
        self.hubbard_N.compute_eigenvalues_and_eigenvectors()

        self.allowed_wavevectors = self.hubbard_N.allowed_wavevectors()

        self.hubbard_plus = None
        self.hubbard_minus = None
        if operator_spin == "up":
            self.hubbard_plus = Hubbard(N_up+1, N_down, parameters)
            self.hubbard_minus = Hubbard(N_up-1, N_down, parameters)
        elif operator_spin == "down":
            self.hubbard_plus = Hubbard(N_up, N_down+1, parameters)
            self.hubbard_minus = Hubbard(N_up, N_down-1, parameters)
        else:
            print("Invalid spin")
        self.hubbard_plus.form_hamiltonian_matrix()
        self.hubbard_plus.compute_eigenvalues_and_eigenvectors()
        self.hubbard_minus.form_hamiltonian_matrix()
        self.hubbard_minus.compute_eigenvalues_and_eigenvectors()

        self.storage = Storage(self)
        if operator_spin == "up":
            self.storage.form_c_dagger_up_ground_eigenstates()
            self.storage.form_c_up_ground_eigenstates()
        elif operator_spin == "down":
            self.storage.form_c_dagger_down_ground_eigenstates()
            self.storage.form_c_down_ground_eigenstates()
        else:
            print("Invalid spin")

        self.green = Green(self)
        self.self_energy = SelfEnergy(self.green)

        self.denominator = Denominator(self.green)

        self.denominator_for_k_n = None
        self.numerator_for_k_n = None
        self.derivative_for_k_n = None
    def get_spectral_function(self,i,j):
        if self.operator_spin == "up":
            return self.green.spectral_function_up(i,j)
        elif self.operator_spin == "down":
            return self.green.spectral_function_down(i,j)
        else:
            print("Invalid operator_spin")
    def plot_spectral_function(self,i,j):

        greater_angular_frequency_differences, greater_weights, lesser_angular_frequency_differences, lesser_weights = self.get_spectral_function(i,j)
        
        print("sum before combine",sum(lesser_weights), sum(greater_weights))
        lesser_angular_frequency_differences,lesser_weights = combine_exact_data(lesser_angular_frequency_differences,lesser_weights)
        greater_angular_frequency_differences,greater_weights = combine_exact_data(greater_angular_frequency_differences,greater_weights)

        # plt.scatter(lesser_angular_frequency_differences, lesser_weights,color="lightgreen")
        # plt.scatter(greater_angular_frequency_differences, greater_weights,color="darkgreen")
        # plt.show()

        plt.vlines(x = lesser_angular_frequency_differences, ymin = zeros(len(lesser_angular_frequency_differences)), ymax = lesser_weights, colors = 'darkblue')
        plt.vlines(x = greater_angular_frequency_differences, ymin = zeros(len(greater_angular_frequency_differences)), ymax = greater_weights, colors = 'darkblue')
        plt.xlabel(r"$\omega$")
        plt.ylabel(r"$A_{ij\sigma}(\omega)$")
        plt.savefig("Exact/Figures/spectral_function.png",dpi=800)
        plt.show()

        print("sum after combine",sum(lesser_weights), sum(greater_weights))
    def get_spectral_function_for_one_wavevector(self,k_n):
        if self.operator_spin == "up":
            return self.green.spectral_function_up_for_one_wavevector(k_n)
        elif self.operator_spin == "down":
            return self.green.spectral_function_down_for_one_wavevector(k_n)
        else:
            print("Invalid operator_spin")
    def plot_spectral_function_for_wavevectors(self):
        H_test_1 = 0
        H_test_2 = 0

        for n,k_n in enumerate(self.allowed_wavevectors):
            greater_angular_frequency_differences, greater_weights, lesser_angular_frequency_differences, lesser_weights = self.get_spectral_function_for_one_wavevector(k_n)
            
            epsilon_of_k_n = self.hubbard_N.band_structure(k_n)
            for g,greater_weight in enumerate(greater_weights):
                H_test_1 = H_test_1 - greater_weight * (greater_angular_frequency_differences[g]  ) / (2 * 2)
            for l,lesser_weight in enumerate(lesser_weights):
                H_test_2 = H_test_2 + lesser_weight * (lesser_angular_frequency_differences[l]  ) / (2 * 2)


            if n==0:
                color = "pink"
            elif n==1:
                color = "orange"
            elif n==2:
                color = "blue"
            elif n==3:
                color = "purple"
            else:
                color = "black"
            

            plt.figure(1)
            plt.scatter(lesser_angular_frequency_differences, lesser_weights,label=str(n),marker='*',c=color)
            plt.scatter(greater_angular_frequency_differences, greater_weights,label=str(n),marker='^',c=color)
            plt.xlabel(r"$\omega$")
            plt.ylabel(r"$A_{\sigma}(\omega)$")
            plt.legend()
            plt.savefig("Exact/Figures/spectral_function_for_wavevectors.png",dpi=800)
            print("sum",sum(lesser_weights).real, sum(greater_weights).real)
        # plt.show()
        print("H_test_1",H_test_1)
        print("H_test_2",H_test_2)
        print("H_test", (H_test_1+H_test_2)*self.parameters.V + self.parameters.mu*(self.N_up + self.N_down) )
        return (H_test_1+H_test_2)*self.parameters.V + self.parameters.mu*(self.N_up + self.N_down)
    def plot_noninteracting_based_on_band_structure(self):
        for n,k_n in enumerate(self.allowed_wavevectors):
            G_0_of_k_n_sigma = self.green.noninteracting_from_bandstructure(k_n)

            w = symbols('w',real=True)
            G_0_of_k_n_sigma_lamdified = lambdify(w, G_0_of_k_n_sigma, modules=['numpy'])

            w_values = linspace(-5,5,1000)
            G_0_of_k_n_sigma_values = G_0_of_k_n_sigma_lamdified(w_values)

            my_linestyle = "solid"
            if n == 4 or n == 5:
                my_linestyle = "dashed"


            plt.figure(1)
            plt.plot( w_values, G_0_of_k_n_sigma_values, label="$k_n=$"+str(round(k_n,2)), linestyle = my_linestyle)
            plt.xlabel(r"$\omega$")
            plt.ylabel(r"$\mathfrak{Re} [ \lim_{T \rightarrow 0^+} G_{0_{k_n k_n \sigma}}(\omega) ]$")

        plt.legend()
        plt.savefig("Exact/Figures/noninteracting_based_on_band_structure.png",dpi=800)
        plt.show()    
    def get_green_real_part_in_k_space(self,k_n):
        if self.operator_spin == "up":
            return self.green.real_part_in_k_space_up(k_n)
        elif self.operator_spin == "down":
            return self.green.real_part_in_k_space_down(k_n)
        else:
            print("Invalid operator_spin")
    def plot_green_real_part_in_k_space(self):

        for n,k_n in enumerate(self.allowed_wavevectors):

            Re_G = self.get_green_real_part_in_k_space(k_n)

            w = symbols('w',real=True)
            Re_G_lamdified = lambdify(w, Re_G, modules=['numpy'])

            w_values = linspace(-5,5,1000)
            Re_G_values = Re_G_lamdified(w_values)

            my_linestyle = "solid"
            if n == 4 or n == 5:
                my_linestyle = "dashed"
            
            plt.figure(2)
            plt.plot( w_values, Re_G_values, label="$k_n=$"+str(round(k_n,2)), linestyle = my_linestyle)
            plt.xlabel(r"$\omega$")
            plt.ylabel(r"$\mathfrak{Re} [ \lim_{T \rightarrow 0^+} G_{k_n k_n \sigma}(\omega) ]$")
            plt.legend()
            plt.savefig("Exact/Figures/green_real_part_in_k_space.png",dpi=800)
        plt.show() 
    def get_constant(self, k_n):
        if self.operator_spin == "up":
            return self.self_energy.constant_term_up(k_n)
        elif self.operator_spin == "down":
            return self.self_energy.constant_term_down(k_n)
        else:
            print("Invalid operator_spin")
    def print_constant(self):
        for n,k_n in enumerate(self.allowed_wavevectors):
            print("k_n = "+str(k_n)+" and constant = "+str(self.get_constant(k_n)))
    def form_denominator(self,k_n):
        if self.operator_spin == "up":
            return self.denominator.form_up(k_n)
        elif self.operator_spin == "down":
            return self.denominator.form_down(k_n)
        else:
            print("Invalid operator_spin")
    def form_numerator_and_denominator_and_derivative(self,k_n):
        w = symbols('w',real=True)
        product_l = self.denominator.product_l(k_n)
        product_g = self.denominator.product_g(k_n)

        epsilon_of_k_n = self.green.hubbard_N.band_structure(k_n)
        denominator = self.form_denominator(k_n)
        numerator = (w - epsilon_of_k_n) * denominator - product_l * product_g
        derivative = None #Polynomial_in_w(denominator).compute_symbolic_derivative_in_w()

        return numerator,denominator,derivative
    def set_numerator_and_denominator_and_derivative(self,numerator,denominator,derivative):
        self.denominator_for_k_n = denominator
        self.numerator_for_k_n = numerator
        self.derivative_for_k_n = derivative
    def get_denominator(self,symbolic_variable):
        w = symbols('w',real=True)
        return self.denominator_for_k_n.subs(w,symbolic_variable)
    def get_derivative_of_denominator(self,symbolic_variable):
        w = symbols('w',real=True)
        return self.derivative_for_k_n.subs(w,symbolic_variable)
    def get_numerator(self,symbolic_variable):
        w = symbols('w',real=True)
        return self.numerator_for_k_n.subs(w,symbolic_variable)

        
