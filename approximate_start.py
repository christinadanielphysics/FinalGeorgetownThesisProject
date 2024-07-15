from Exact.Problem import *
from Exact.Parameters import *
from Approximate.TimeEvolution import *
from Approximate.Compressive_Sensing_Denominator import *
from Exact.Root_Finding.real_roots import get_real_roots_with_multiplicities
from Approximate.Compressive_Sensing_Green import *

from numpy import linspace,transpose,array
from cmath import exp
from sympy import symbols,nroots,apart,simplify,plot,cancel,diff
import matplotlib.pyplot as plt

from cosamp import cosamp
from numpy import floor,identity,pi,arange,conj
from numpy.random import rand
from scipy.fftpack import dct,idct
from numpy.fft import fft,fftfreq

import cvxpy as cvx

from scipy import sparse
from scipy.sparse.linalg import lsqr
from numpy.linalg import norm
import numpy

from scipy.fftpack import fft
from scipy.signal import find_peaks



class CompressiveSensingCosineTransform:
    def __init__(self,n,p,perm,number_of_peaks,sampling_f,V):
        self.V = V
        self.sampling_f = sampling_f
        self.number_of_peaks = number_of_peaks
        self.perm = perm
        # Psi = dct( identity(n) ) # build Psi
        # self.Theta = Psi[perm,:] # build Theta with rows of Psi
        self.n = n
        self.p = p
        self.L = int( floor(n/2) )
    def compute(self,x): # reconstruct the signal using compressive sensing
        # Number of sample points
        n = self.n
        # sample spacing
        T = 1.0 / self.sampling_f
        yf = fft(x)
        return (1/n) *  numpy.abs(yf[0:n//2])
    def plot(self,y,w):
        plt.plot(w, y)
        plt.grid()
        plt.show()
    def chop(self,all_w_values,all_weights): 
        # Find the peaks with a minimum prominence of 0.000001
        peaks = find_peaks( numpy.abs(all_weights), prominence=0.000001 )[0]
        return all_w_values[peaks], all_weights[peaks]

def constant_term(problem,k_n,greater_weights_up_in_k_space,lesser_weights_up_in_k_space,greater_w_differences,lesser_w_differences):
    band_structure = problem.hubbard_N.band_structure(k_n)
        
    greater_piece = 0
    for g,greater_weight_up_in_k_space in enumerate(greater_weights_up_in_k_space):
        greater_piece = greater_piece + greater_weight_up_in_k_space * ( greater_w_differences[g] )
        
    lesser_piece = 0
    for l,lesser_weight_up_in_k_space in enumerate(lesser_weights_up_in_k_space):
        lesser_piece = lesser_piece + lesser_weight_up_in_k_space * ( lesser_w_differences[l]  )
        
    return greater_piece + lesser_piece
def convert_weights_from_position_to_wavevector(input_weights,V,k_n):
    weights_for_k_n = []
    for index, input_weight in enumerate(input_weights):
        weight_for_k_n = 0
        for i in range(V):
            for j in range(V):
                weight_for_k_n = weight_for_k_n + (1/V) * exp(1j * k_n * (i - j)) * input_weight
        weights_for_k_n.append(weight_for_k_n.real)
    return weights_for_k_n
def get_denominator_roots(problem,numerator,denominator,k_n,integer_n):
    # Set the denominator and the numerator so that the root finder can access these polynomials quickly
    derivative = diff(denominator)
    problem.set_numerator_and_denominator_and_derivative(numerator,denominator,derivative)
    print("derivative",derivative)
    # ROOTS AND WEIGHTS for the real part of the self energy
    step = 0.1
    root_objects_denominator = []
    done = 0
    while done == 0:
        print("DENOMINATOR")
        print("k_n value is ", k_n)
        done = int(input("Enter 0 to continue, enter 1 to stop "))
        if done == 1:
            break
        plot(numerator/denominator)
        a_min = float(input("a_min"))
        b_max = a_min + 0.1
        # Find the roots of the denominator of the real part of the self energy
        roots_denominator = get_real_roots_with_multiplicities(problem,a_min,b_max,step,root_objects_denominator)
    return roots_denominator

def is_not_zero(array):
    for element in array:
        if element != 0:
            return True
    return False

def compute(U):
    # Problem Definition
    N_up = 2
    N_down = 2
    operator_spin = "up"
    t_0 = 0
    t_1 = 1
    t_2 = 0
    mu = 0
    V = 4
    parameters = Parameters(t_0,t_1,t_2,U,mu,V)
    problem = Problem(parameters,operator_spin,N_up,N_down)
    ground_state_energy = problem.hubbard_N.ground_state_energy




    # Adjustable compressive sensing parameters (there is also a parameter in the chop function)
    n = int(1000) # affects highest frequency magnitude that can obtain, also affects "angular frequency resolution" 
    highest_expected_angular_frequency = 9 # also affects "angular frequency resolution"




    p = 10
    number_of_expected_peaks = 5
    sampling_w = 2 * highest_expected_angular_frequency
    perm = floor( rand(p) * n ).astype(int) # random indices
    L = int( floor(n/2) )
    sampling_f = sampling_w / (2 * pi) # should be twice the highest expected frequency in Hz
    sampling_period = 1/sampling_f
    end_time = sampling_period * (n-1) 
    print("2 pi /end time",2 * pi/end_time)

    t = linspace(0,end_time,n)
    f = linspace(0, (n-1)*sampling_f/(2*n), n) # fft
    # f = linspace(0, (n-1)*sampling_f/(n), n)
    # calculate the frequency axis
    # f = fftfreq(n, d=1/sampling_f)
    f = linspace(0.0, 1.0/(2.0*sampling_period), n//2) # just fft
    delta_f = f[1] - f[0]
    print("delta_w",delta_f * 2 * pi)
    w = 2 * pi * f





    # all_w_values_greater = w[:L] # fft

    all_w_values_greater = w 
    all_w_values_lesser = w * (-1)
    
    
    time_evolution = TimeEvolution(problem,t)
    greater_weights_for_all_k_n = {}
    lesser_weights_for_all_k_n = {}
    w_values_greater_for_all_k_n = {}
    w_values_lesser_for_all_k_n = {}

    

    for integer_n,k_n in enumerate(problem.allowed_wavevectors):
        print("integer_n",integer_n)

        lesser_real, lesser_imag = time_evolution.time_evolve_lesser_in_k_space(k_n)
        greater_real, greater_imag = time_evolution.time_evolve_greater_in_k_space(k_n)
        # time_evolution.plot_lesser_up_in_k_space(k_n)
        # time_evolution.plot_greater_up_in_k_space(k_n)

        # # # Get angular frequency differences and weights from time-dependent data
        cosine_transform_object = CompressiveSensingCosineTransform(n,p,perm,number_of_expected_peaks,sampling_f,V)
        
        all_weights_lesser_even = []
        all_weights_greater_even = []
        if is_not_zero(lesser_imag) == True:
            all_weights_lesser_even = cosine_transform_object.compute(lesser_imag)     
            print("integer n",integer_n,"original sum",sum(all_weights_lesser_even))
        if is_not_zero(greater_imag) == True:
            all_weights_greater_even = cosine_transform_object.compute(-greater_imag) # just want the weight, not negative of the weight, which is apparent from the formula for the t-dependent greater green's function which has a negative sign out in front
            print("integer n",integer_n,"original sum",sum(all_weights_greater_even))
        
        # only save/use nonzero weights with corresponding frequency difference
        plt.figure(1)
        weights_greater_even = []
        w_values_greater_even = []
        weights_lesser_even = []
        w_values_lesser_even = []
        if all_weights_greater_even != []:
            w_values_greater_even, weights_greater_even = cosine_transform_object.chop(all_w_values_greater, all_weights_greater_even)
            #plt.scatter(w_values_greater_even,weights_greater_even,c="green",marker="o")
        if all_weights_lesser_even != []:
            w_values_lesser_even, weights_lesser_even = cosine_transform_object.chop(all_w_values_lesser, all_weights_lesser_even)
            #plt.scatter(w_values_lesser_even,weights_lesser_even,c="red",marker="x")

        # normalize last for each k-value, enforce the SUM of lesser and greater weights being 1. This is the sum rule for a spectral function that depends on one k-value. 
        sum_of_lesser_and_greater_weights = 0
        for weight in weights_lesser_even:
            sum_of_lesser_and_greater_weights = sum_of_lesser_and_greater_weights + weight
        for weight in weights_greater_even:
            sum_of_lesser_and_greater_weights = sum_of_lesser_and_greater_weights + weight
        for index,weight in enumerate(weights_lesser_even):
            weights_lesser_even[index] = weights_lesser_even[index]/sum_of_lesser_and_greater_weights
        for index,weight in enumerate(weights_greater_even):
            weights_greater_even[index] = weights_greater_even[index]/sum_of_lesser_and_greater_weights
        print("integer_n",integer_n,"sum of CS lesser weights",sum(weights_lesser_even),"sum of CS greater weights", sum(weights_greater_even))


        # save weights and frequencies for each k-value with dictionary
        greater_weights_for_all_k_n[integer_n] = weights_greater_even
        w_values_greater_for_all_k_n[integer_n] = w_values_greater_even
        lesser_weights_for_all_k_n[integer_n] = weights_lesser_even
        w_values_lesser_for_all_k_n[integer_n] = w_values_lesser_even
    
    
    #problem.plot_spectral_function_for_wavevectors()



    
    green_object = Compressive_Sensing_Green(problem, greater_weights_for_all_k_n, lesser_weights_for_all_k_n, w_values_greater_for_all_k_n, w_values_lesser_for_all_k_n)
    #green_object.plot_real_part_in_k_space()




    
    T_expectation = 0
    for integer_n,k_n in enumerate(problem.allowed_wavevectors):
        for g,weight in enumerate(greater_weights_for_all_k_n[integer_n]):
            T_expectation = T_expectation - weight * problem.hubbard_N.band_structure(k_n)
    for integer_n,k_n in enumerate(problem.allowed_wavevectors):
        for l,weight in enumerate(lesser_weights_for_all_k_n[integer_n]):
            T_expectation = T_expectation + weight * problem.hubbard_N.band_structure(k_n)
    
    print("T_expectation",T_expectation)
    


    self_energy_weights_for_all_k_n = {}
    denominator_roots_for_all_k_n = {}
    sigma_at_greater_poles_for_all_k_n = {}
    sigma_at_lesser_poles_for_all_k_n = {}
    for integer_n,k_n in enumerate(problem.allowed_wavevectors):

        print("\n\ninteger_n =",integer_n)
        
        # Compute the noninteracting Green's function based on the bandstructure
        G0_for_k_n = problem.green.noninteracting_from_bandstructure(k_n)
        
        # Compute the constant
        constant = constant_term(problem, k_n, greater_weights_for_all_k_n[integer_n], lesser_weights_for_all_k_n[integer_n], w_values_greater_for_all_k_n[integer_n], w_values_lesser_for_all_k_n[integer_n])
        # print("constant for k_n:", constant)

        
        # Use lesser and greater weights to form the denominator of the real part of the self energy
        object_for_denominator = Compressive_Sensing_Denominator(greater_weights_for_all_k_n[integer_n], lesser_weights_for_all_k_n[integer_n], w_values_greater_for_all_k_n[integer_n], w_values_lesser_for_all_k_n[integer_n])
        denominator = object_for_denominator.form()
        print("denominator",denominator)

        
        # Use lesser and greater weights to form the numerator of the real part of the self energy
        # Subtract and add the constant
        w = symbols('w',real=True)
        band_structure_for_k_n = problem.hubbard_N.band_structure(k_n)
        numerator = (w - band_structure_for_k_n - constant) * denominator - object_for_denominator.product_over_l() * object_for_denominator.product_over_g()
        print("numerator",numerator)

        if denominator == 0:
            print("denominator is zero")
            denominator = 1e-17*w
        denominator_roots = nroots(denominator,maxsteps=200000)
        denominator_roots_for_all_k_n[integer_n] = denominator_roots
        print("denominator_roots",denominator_roots)
        # denominator_roots_for_all_k_n[integer_n] = get_denominator_roots(problem,numerator,denominator,k_n,integer_n)



        # see what the real part of the self energy is at poles of the green's function
        sigma_at_lesser_poles = []
        sigma_at_greater_poles = []
        for w_g_minus_w_0 in w_values_greater_for_all_k_n[integer_n]:
            real_part_of_self_energy = numerator/denominator
            value = simplify(real_part_of_self_energy.subs(w,w_g_minus_w_0))
            sigma_at_greater_poles.append(value)
        for w_0_minus_w_l in w_values_lesser_for_all_k_n[integer_n]:
            real_part_of_self_energy = numerator/denominator
            value = simplify(real_part_of_self_energy.subs(w,w_0_minus_w_l))
            sigma_at_lesser_poles.append(value)
        sigma_at_lesser_poles_for_all_k_n[integer_n] = sigma_at_lesser_poles
        sigma_at_greater_poles_for_all_k_n[integer_n] = sigma_at_greater_poles


        # make sure none of the self energy poles have degeneracy!!



        denominator_weights = []
        numerator_weights = []
        weights = []

        # form a new denominator based on the roots of the denominator
        new_denominator = 1
        for root in denominator_roots:
            new_denominator = new_denominator * (w - root)
        print("new denominator",new_denominator)
        #plot(numerator/new_denominator) # ensure the poles are visible


        for root in denominator_roots:

            # cancel the factor with this denominator_root by forming a new denominator without that factor
            denominator_expression = 1
            for den_root in denominator_roots:
                if den_root != root:
                    denominator_expression = denominator_expression * (w - den_root)
           
            print("denominator expression",denominator_expression)
            denominator_weight = 1
            if denominator_expression != 1:
                denominator_weight = simplify(denominator_expression.subs(w,root))
            numerator_weight = simplify(numerator.subs(w,root))
            print("numerator_weight",numerator_weight,"denominator_weight",denominator_weight)
            numerator_weights.append(numerator_weight)
            weight = numerator_weight/denominator_weight
            weights.append(weight)
        self_energy_weights_for_all_k_n[integer_n] = weights

        # sum rule for each k-value
        local_moment = 0
        for weight in weights:
            local_moment = local_moment + weight
        print("The local moment for k_n=",k_n,"is",local_moment)


    # Compute Expectation Values
    term_1 = 0
    for integer_n,k_n in enumerate(problem.allowed_wavevectors):
        roots_for_one_k_value = denominator_roots_for_all_k_n[integer_n]
        for r,c_r in enumerate(self_energy_weights_for_all_k_n[integer_n]):
            if roots_for_one_k_value[r] <= TOLERANCE_computational_zero_for_denominator:
                if abs(roots_for_one_k_value[r]) <= TOLERANCE_computational_zero_for_denominator:
                    term_1 = term_1 + c_r * green_object.real_part_in_k_space(integer_n).subs(w,roots_for_one_k_value[r]) / 2
                else:
                    term_1 = term_1 + c_r * green_object.real_part_in_k_space(integer_n).subs(w,roots_for_one_k_value[r])
        print("for k_n =",k_n,"term_1 =",term_1)
    term_3 = 0
    for integer_n,k_n in enumerate(problem.allowed_wavevectors):
        w_0_minus_w_l_values = w_values_lesser_for_all_k_n[integer_n]
        for l,lesser_weight in enumerate(lesser_weights_for_all_k_n[integer_n]):
            if w_0_minus_w_l_values[l] <= TOLERANCE_computational_zero_for_denominator:
                if abs(lesser_weight) >= TOLERANCE_computational_zero_for_denominator:
                    if abs(w_0_minus_w_l_values[l]) <= TOLERANCE_computational_zero_for_denominator:
                        sigma_values = sigma_at_lesser_poles_for_all_k_n[integer_n]
                        term_3 = term_3 + lesser_weight * sigma_values[l] / 2
                    else:
                        sigma_values = sigma_at_lesser_poles_for_all_k_n[integer_n]
                        term_3 = term_3 + lesser_weight * sigma_values[l]
        print("for k_n =",k_n,"term_3 =",term_3)
                    

    print("term 1",term_1)
    print("term 3",term_3)

    U_expectation = (term_1 + term_3 * pi/2)
    N_expectation = mu *(N_up + N_down)
    H_expectation = T_expectation + U_expectation

    print("U summation",U_expectation)
    print("T_expectation =",T_expectation)
    print("H_expectation =",H_expectation)
    print("exact ground state energy",ground_state_energy)



    plt.figure(5)
    plt.scatter(U,U_expectation,color='red',label='calculated')
    plt.scatter(U,ground_state_energy - T_expectation,color='blue',marker='x',label='exact')
    plt.scatter(U,U_expectation/(ground_state_energy - T_expectation ),label='ratio',marker='.',color='green')


step = 2
U = 1
for index in range(10):
    compute(U+0.5*index)

plt.legend()
plt.show()


