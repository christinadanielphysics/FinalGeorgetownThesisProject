from sympy import *
from .tolerances import *

import numpy
import matplotlib.pyplot as plt
from cmath import pi

class Denominator:
    def __init__(self,green):
        self.green = green
        self.mu = green.problem.parameters.mu
        self.U = green.problem.parameters.U
        self.N_up = green.problem.N_up
        self.N_down = green.problem.N_down
        self.allowed_wavevectors = green.allowed_wavevectors
        self.ground_state_energy = green.ground_state_energy
        self.computational_zero = TOLERANCE_computational_zero_for_denominator_approximate # for weights close to zero
        # note: filter_lesser and filter_greater are for handling degenerate eigenvalues and their weights
    def filter_lesser(self,k_n):
        original_minus_dimension = self.green.minus_dimension
        original_minus_eigenvalues = self.green.minus_eigenvalues
        original_lesser_weights_for_k_n = []
        for l in range(original_minus_dimension):
            original_lesser_weights_for_k_n.append( self.green.lesser_weight_up_in_k_space(l,k_n) )

        new_minus_dimension = 0
        new_minus_eigenvalues = []
        new_lesser_weights_for_k_n = []
        unique_l = 0
        for l,value in enumerate(original_minus_eigenvalues):
            length_of_new_list = len(new_lesser_weights_for_k_n)
            if length_of_new_list != 0:
                if abs(value - new_minus_eigenvalues[length_of_new_list-1]) <= TOLERANCE_combine:
                    new_lesser_weights_for_k_n[length_of_new_list-1] = new_lesser_weights_for_k_n[length_of_new_list-1] + original_lesser_weights_for_k_n[l]
                else:
                    if abs(original_lesser_weights_for_k_n[l]) >= self.computational_zero: # only include significant weights and frequencies
                        new_minus_eigenvalues.append(value)
                        new_minus_dimension = new_minus_dimension + 1
                        new_lesser_weights_for_k_n.append( original_lesser_weights_for_k_n[l] )    
            else:
                if abs(original_lesser_weights_for_k_n[l]) >= self.computational_zero: # only include significant weights and frequencies
                    new_minus_eigenvalues.append(value)
                    new_minus_dimension = new_minus_dimension + 1
                    new_lesser_weights_for_k_n.append( original_lesser_weights_for_k_n[l] )

        return new_minus_dimension, new_minus_eigenvalues, new_lesser_weights_for_k_n
    def filter_greater(self,k_n):
        original_plus_dimension = self.green.plus_dimension
        original_plus_eigenvalues = self.green.plus_eigenvalues
        original_greater_weights_for_k_n = []
        for g in range(original_plus_dimension):
            original_greater_weights_for_k_n.append( self.green.greater_weight_up_in_k_space(g,k_n) )
        
        new_plus_dimension = 0
        new_plus_eigenvalues = []
        new_greater_weights_for_k_n = []
        unique_g = 0
        for g,value in enumerate(original_plus_eigenvalues):
            length_of_new_list = len(new_greater_weights_for_k_n)
            if length_of_new_list != 0:
                if abs(value - new_plus_eigenvalues[length_of_new_list-1]) <= TOLERANCE_combine:
                    new_greater_weights_for_k_n[length_of_new_list-1] = new_greater_weights_for_k_n[length_of_new_list-1] + original_greater_weights_for_k_n[g]
                else:
                    if abs(original_greater_weights_for_k_n[g]) >= self.computational_zero: # only include significant weights and frequencies
                        new_plus_eigenvalues.append(value)
                        new_plus_dimension = new_plus_dimension + 1
                        new_greater_weights_for_k_n.append( original_greater_weights_for_k_n[g] )
            else:
                if abs(original_greater_weights_for_k_n[g]) >= self.computational_zero: # only include significant weights and frequencies
                    new_plus_eigenvalues.append(value)
                    new_plus_dimension = new_plus_dimension + 1
                    new_greater_weights_for_k_n.append( original_greater_weights_for_k_n[g] )

        return new_plus_dimension, new_plus_eigenvalues, new_greater_weights_for_k_n
    def product_l(self,k_n):
        w = symbols('w',real=True)
        w_0 = self.ground_state_energy

        minus_dimension, minus_eigenvalues, lesser_weight_up_in_k_space = self.filter_lesser(k_n)


        product_l = 1
        for l in range(minus_dimension):
            W_l = minus_eigenvalues[l] - w_0
            product_l = product_l * (w + W_l)
        return product_l
    def product_m_not_l(self,l,k_n):
        w = symbols('w',real=True)
        w_0 = self.ground_state_energy

        minus_dimension, minus_eigenvalues, lesser_weight_up_in_k_space = self.filter_lesser(k_n)

        product_m_not_l = 1
        for m in range(minus_dimension):
            if m != l:
                W_m = minus_eigenvalues[m] - w_0
                product_m_not_l = product_m_not_l * (w + W_m)
        return product_m_not_l
    def product_g(self,k_n):
        w = symbols('w',real=True)
        w_0 = self.ground_state_energy

        plus_dimension, plus_eigenvalues, greater_weight_up_in_k_space = self.filter_greater(k_n)

        product_g = 1
        for g in range(plus_dimension):
            W_g = w_0 - plus_eigenvalues[g]
            product_g = product_g * (w + W_g)
        return product_g
    def product_h_not_g(self,g,k_n):
        w = symbols('w',real=True)
        w_0 = self.ground_state_energy

        plus_dimension, plus_eigenvalues, greater_weight_up_in_k_space = self.filter_greater(k_n)


        product_h_not_g = 1
        for h in range(plus_dimension):
            if h != g:
                W_h = w_0 - plus_eigenvalues[h]
                product_h_not_g = product_h_not_g * (w + W_h)
        return product_h_not_g
    def form_up(self,k_n):
        minus_dimension, minus_eigenvalues, lesser_weight_up_in_k_space = self.filter_lesser(k_n)
        plus_dimension, plus_eigenvalues, greater_weight_up_in_k_space = self.filter_greater(k_n)

        print("NUMBER OF LESSER POLES",len(minus_eigenvalues))
        print("NUMBER OF GREATER POLES",len(plus_eigenvalues))


        sum_g = 0
        for g in range(plus_dimension):
            if abs(greater_weight_up_in_k_space[g].real) >= self.computational_zero:
                sum_g = sum_g + greater_weight_up_in_k_space[g].real * self.product_h_not_g(g,k_n)
        term_1 = sum_g * self.product_l(k_n)

        sum_l = 0
        for l in range(minus_dimension):
            if abs(lesser_weight_up_in_k_space[l].real) >= self.computational_zero:
                sum_l = sum_l + lesser_weight_up_in_k_space[l].real * self.product_m_not_l(l,k_n)
        term_2 = sum_l * self.product_g(k_n)

        return term_1 + term_2
    def get_real_self_energy(self,n,value,all_numerators,all_denominators):
        w = symbols('w',real=True)
        expression = all_numerators[n]/all_denominators[n]
        result = expression.subs(w,value)
        return result
    def get_real_green(self,k_n,value):
        w = symbols('w',real=True)
        w_0 = self.ground_state_energy

        minus_dimension, minus_eigenvalues, lesser_weight_up_in_k_space = self.filter_lesser(k_n)
        plus_dimension, plus_eigenvalues, greater_weight_up_in_k_space = self.filter_greater(k_n)

        Re_G = 0
        for g in range(plus_dimension):
            if abs(greater_weight_up_in_k_space[g].real) >= self.computational_zero:
                Re_G = Re_G - greater_weight_up_in_k_space[g].real / (w - plus_eigenvalues[g] + w_0)
        for l in range(minus_dimension):
            if abs(lesser_weight_up_in_k_space[l].real) >= self.computational_zero:
                Re_G = Re_G + lesser_weight_up_in_k_space[l].real / (w + minus_eigenvalues[l] - w_0) 
    
        result = Re_G.subs(w,value) / (2 * pi)

        return result
    def potential(self,root_objects_for_all_n,weights_for_all_n,all_numerators,all_denominators):
        w_0 = self.ground_state_energy
        
        term_1 = 0
        for n,k_n in enumerate(self.allowed_wavevectors):
            roots_for_n = root_objects_for_all_n[n]
            for r,c_r in enumerate(weights_for_all_n[n]):
                if roots_for_n[r] <= 0:
                    if abs(roots_for_n[r]) <= TOLERANCE_computational_zero_for_denominator:
                        term_1 = term_1 + c_r * ( self.get_real_green(k_n, roots_for_n[r]) ) / 2
                    else:
                        term_1 = term_1 + c_r * ( self.get_real_green(k_n, roots_for_n[r]) )

            print("k_n =",k_n,"and term_1 =",term_1)


        term_3 = 0
        for n,k_n in enumerate(self.allowed_wavevectors):
            minus_dimension, minus_eigenvalues, lesser_weight_up_in_k_space = self.filter_lesser(k_n)
            for l in range(minus_dimension):
                w_l = minus_eigenvalues[l]
                if w_0 - w_l  <= 0:
                    if abs(lesser_weight_up_in_k_space[l].real) >= self.computational_zero:
                        if abs(w_0-w_l) <= TOLERANCE_computational_zero_for_denominator:                            
                            term_3 = term_3 + self.get_real_self_energy(n,w_0-w_l,all_numerators,all_denominators) * lesser_weight_up_in_k_space[l].real /2
                        else:
                            term_3 = term_3 + self.get_real_self_energy(n,w_0-w_l,all_numerators,all_denominators) * lesser_weight_up_in_k_space[l].real  
            print("k_n =",k_n,"and term_3 =",term_3)

        print("term 1",term_1)
        print("term 3",term_3)
        
        return (term_1 + term_3 * pi/2 ) 

