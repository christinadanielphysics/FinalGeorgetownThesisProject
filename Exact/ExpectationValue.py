from Exact.Parameters import *
from Exact.Problem import *

from sympy import symbols,simplify,lambdify
import math
import numpy

class ExpectationValue:
    def __init__(self,all_self_energy_roots,all_self_energy_weights,problem,all_numerators,all_denominators):
        self.all_self_energy_roots = all_self_energy_roots
        self.all_self_energy_weights = all_self_energy_weights
        self.allowed_wavevectors = problem.allowed_wavevectors
        self.factor_for_spin = 2
        self.hubbard_N = problem.hubbard_N
        self.green = problem.green
        self.lesser_dimension = self.green.minus_dimension
        self.greater_dimension = self.green.plus_dimension
        self.plus_eigenvalues = problem.hubbard_plus.eigenvalues
        self.minus_eigenvalues = problem.hubbard_minus.eigenvalues
        self.all_numerators = all_numerators
        self.all_denominators = all_denominators
        self.ground_state_energy = problem.hubbard_N.ground_state_energy
        self.U = problem.parameters.U
    def get_real_green(self,k_n,value):
        w = symbols('w',real=True)
        expression = self.green.real_part_in_k_space_up(k_n) 
        plt.figure(5)
        result = expression.subs(w,value)
        plt.scatter(value,re(result))
        return result
    def get_real_self_energy(self,n,value):
        w = symbols('w',real=True)
        expression = self.all_numerators[n]/self.all_denominators[n]
        result = expression.subs(w,value)
        return result
    def get_real_green_lesser(self,k_n,value):
        w = symbols('w',real=True)
        w_0 = self.green.ground_state_energy

        expression = 0
        for l in range(self.green.minus_dimension):
            expression = expression + self.green.lesser_weight_up_in_k_space(l, k_n).real / (w + self.green.minus_eigenvalues[l] - w_0)
        
        result = expression.subs(w,value)
        return result
    def get_real_green_greater(self,k_n,value):
        w = symbols('w',real=True)
        w_0 = self.green.ground_state_energy

        expression = 0
        for g in range(self.green.plus_dimension):
            expression = expression + self.green.greater_weight_up_in_k_space(g, k_n).real / (w - self.plus_eigenvalues[g] + w_0)

        
        result = expression.subs(w,value)
        return result
    def get_noninteracting(self,k_n,value):
        w = symbols('w',real=True)
        epsilon_of_k_n = self.green.hubbard_N.band_structure(k_n)
        G_0_of_k_n_sigma = 1/(w - epsilon_of_k_n)

        result = G_0_of_k_n_sigma.subs(w,value)

        return result
    def potential(self): # combine the weights for the green's function ; need nonzero weights for the green's function
        w_0 = self.ground_state_energy
        term_1 = 0
        for n,k_n in enumerate(self.allowed_wavevectors):
            roots_for_n = self.all_self_energy_roots[n]
            print("k_n",k_n,"poles of self energy",numpy.sort(roots_for_n))
            for r,c_r in enumerate(self.all_self_energy_weights[n]):
                if roots_for_n[r] < 0: # only looking at negative frequencies because of the fermi-dirac distribution
                    term_1 = term_1 - c_r * ( self.get_real_green(k_n, roots_for_n[r]) )
        print("term1",term_1)
        # term_2 = 0
        # for n,k_n in enumerate(self.allowed_wavevectors):
        #     for g in range(self.greater_dimension):
        #         w_g = self.plus_eigenvalues[g]
        #         term_2 = term_2 + self.get_real_self_energy(n,w_g - w_0) * self.green.greater_weight_up_in_k_space(g,k_n).real
        term_3 = 0
        for n,k_n in enumerate(self.allowed_wavevectors):
            for l in range(self.lesser_dimension):
                w_l = self.minus_eigenvalues[l]
                print("k_n",k_n,"pole of lesser green",w_0-w_l)
                term_3 = term_3 + self.get_real_self_energy(n,w_0-w_l) * self.green.lesser_weight_up_in_k_space(l,k_n).real
        print("term_3",term_3)
        print("potential",(term_1 + term_3) * pi)
        return (term_1 + term_3) * pi 
    def potential_guess(self):
        w_0 = self.ground_state_energy
        term_2 = 0
        for n,k_n in enumerate(self.allowed_wavevectors):
            for g in range(self.greater_dimension):
                w_g = self.plus_eigenvalues[g]
                term_2 = term_2 - ( self.green.greater_weight_up_in_k_space(g,k_n).real / self.get_noninteracting(k_n,w_g - w_0)  )
        term_3 = 0
        for n,k_n in enumerate(self.allowed_wavevectors):
            for l in range(self.lesser_dimension):
                w_l = self.minus_eigenvalues[l]
                term_3 = term_3 + ( self.green.lesser_weight_up_in_k_space(l,k_n).real / self.get_noninteracting(k_n,w_0-w_l) )
        print("potential guess",term_2,term_3,term_2 +term_3)
        return term_2 +term_3
    def potential_lesser(self):
        w_0 = self.ground_state_energy
        term_1 = 0
        for n,k_n in enumerate(self.allowed_wavevectors):
            roots_for_n = self.all_self_energy_roots[n]
            for r,c_r in enumerate(self.all_self_energy_weights[n]):
                term_1 = term_1 + c_r * self.get_real_green_lesser(k_n,roots_for_n[r]) 
        term_2 = 0
        # for n,k_n in enumerate(self.allowed_wavevectors):
        #     for g in range(self.greater_dimension):
        #         w_g = self.plus_eigenvalues[g]
        #         term_2 = term_2 - self.get_real_self_energy(n,w_g - w_0) * self.green.greater_weight_up_in_k_space(g,k_n)
        term_3 = 0
        # for n,k_n in enumerate(self.allowed_wavevectors):
        #     for l in range(self.lesser_dimension):
        #         w_l = self.minus_eigenvalues[l]
        #         term_3 = term_3 + self.get_real_self_energy(n,w_0-w_l) * self.green.lesser_weight_up_in_k_space(l,k_n)
        print("potential lesser",term_1)
        return term_1
    def potential_greater(self):
        w_0 = self.ground_state_energy
        term_1 = 0
        for n,k_n in enumerate(self.allowed_wavevectors):
            roots_for_n = self.all_self_energy_roots[n]
            for r,c_r in enumerate(self.all_self_energy_weights[n]):
                term_1 = term_1 + c_r * self.get_real_green_greater(k_n,roots_for_n[r]) 
        print("potential greater",term_1)
        return term_1
    def kinetic(self):
        greater_value = 0
        for n,k_n in enumerate(self.allowed_wavevectors):
            for g in range(self.greater_dimension):
                greater_value = greater_value - self.green.greater_weight_up_in_k_space(g,k_n) * self.hubbard_N.band_structure(k_n)
        lesser_value = 0
        for n,k_n in enumerate(self.allowed_wavevectors):
            for l in range(self.lesser_dimension):
                lesser_value = lesser_value + self.green.lesser_weight_up_in_k_space(l,k_n) * self.hubbard_N.band_structure(k_n)
        return (greater_value + lesser_value) * self.factor_for_spin * (1/2)




