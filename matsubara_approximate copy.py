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

w = symbols('w')

def define_problem():
    N_up = 3
    N_down = 3
    V = 6
    U = 1
    t_1 = 1

    operator_spin = "up"
    t_0 = 0
    t_2 = 0
    mu = U/2

    parameters = Parameters(t_0,t_1,t_2,U,mu,V)
    problem = Problem(parameters,operator_spin,N_up,N_down)

    return problem, U

def plot_time_evolution(problem): 

    # ADJUSTABLE TIME AND FREQUENCY PARAMETERS
    # There is also a parameter in the chop function
    # affects highest frequency magnitude that can obtain
    # also affects "angular frequency resolution" 
    # highest_expected_angular_frequency should be twice the highest expected angular frequency 
    
    n = 300
    highest_expected_angular_frequency = 50

    # Automated Conversions
    sampling_w = 2 * highest_expected_angular_frequency
    sampling_f = sampling_w / (2 * pi) 
    sampling_period = 1/sampling_f
    end_time = sampling_period * (n-1)
    t = linspace(0,end_time,n)

    sub_figure,sub_axes = plt.subplots( len(problem.allowed_wavevectors), sharex=True)
    sub_figure.tight_layout()
    for sub_axis in sub_axes:
        sub_axis.spines['left'].set_linewidth(2)
        sub_axis.spines['right'].set_linewidth(2)
        sub_axis.spines['top'].set_linewidth(2)
        sub_axis.spines['bottom'].set_linewidth(2)
        sub_axis.xaxis.set_tick_params(which='major', size=10, width=2, direction='in', top='on', labelsize=5)
        sub_axis.xaxis.set_tick_params(which='minor', size=7, width=2, direction='in', top='on', labelsize=5)
        sub_axis.yaxis.set_tick_params(which='major', size=10, width=2, direction='in', right='on', labelsize=5)
        sub_axis.yaxis.set_tick_params(which='minor', size=7, width=2, direction='in', right='on', labelsize=5)
        
    imaginary_greater_vs_time = {}
    imaginary_lesser_vs_time = {}
    time_evolution = TimeEvolution(problem,t)

    for integer_n,k_n in enumerate(problem.allowed_wavevectors):

        # Lesser
        g_real, g_imag = time_evolution.time_evolve_lesser_in_k_space(k_n)
        imaginary_lesser_vs_time[integer_n] = g_imag
        # sub_axes[integer_n].plot(time_evolution.t_values,g_real,color="darkred",linewidth=1)
        # sub_axes[integer_n].plot(time_evolution.t_values,g_imag,color="red",linewidth=1)
        # sub_axes[integer_n].set_title(r"$k_n =$"+str(round(k_n,2)),size=5)
        # sub_axes[integer_n].set_xlabel('t', size=10)

        # Greater
        g_real, g_imag = time_evolution.time_evolve_greater_in_k_space(k_n)
        imaginary_greater_vs_time[integer_n] = g_imag
        sub_axes[integer_n].plot(time_evolution.t_values,g_real,color="darkgreen",linewidth=1)
        sub_axes[integer_n].plot(time_evolution.t_values,g_imag,color="yellowgreen",linewidth=1)
        sub_axes[integer_n].set_title(r"$k_n =$"+str(round(k_n,2)),size=5)
        sub_axes[integer_n].set_xlabel('t', size=10)

    # Hide x labels and tick labels for all but bottom plot.
    for sub_axis in sub_axes:
        sub_axis.label_outer()

    plt.savefig('./Figures_for_Thesis/approximate_greater_vs_time_for_6_sites.png', dpi=800, transparent=False, bbox_inches='tight')

    return imaginary_greater_vs_time, imaginary_lesser_vs_time, n, sampling_f

def is_not_zero(array):
    for element in array:
        if element != 0:
            return True
    return False

def transform(n,sampling_f,signal):
    T = 1.0 / sampling_f
    yf = fft(signal)
    return (1/n) *  numpy.abs(yf[0:n//2])

def chop(all_w_values,all_weights): 
    # Find the peaks with a minimum prominence of 0.000001
    peaks = find_peaks( numpy.abs(all_weights), prominence=0.000001 )[0]
    return all_w_values[peaks], all_weights[peaks]

def lesser_and_greater_information(problem, imaginary_greater_vs_time, imaginary_lesser_vs_time, n, sampling_f, U):

    f = linspace(0, (n-1)*sampling_f/n, n) 
    w = 2 * pi * f
    all_w_values_greater = w 
    all_w_values_lesser = w * (-1)

    greater_weights_for_all_k_n = {}
    lesser_weights_for_all_k_n = {}
    w_values_greater_for_all_k_n = {}
    w_values_lesser_for_all_k_n = {}

    # Convert to the angular frequency domain
    for integer_n,k_n in enumerate(problem.allowed_wavevectors):
        all_weights_lesser_even = []
        all_weights_greater_even = []
        if is_not_zero(imaginary_lesser_vs_time[integer_n]) == True:
            all_weights_lesser_even = transform(n,sampling_f,imaginary_lesser_vs_time[integer_n])     
        if is_not_zero(imaginary_greater_vs_time[integer_n]) == True:
            all_weights_greater_even = transform(n,sampling_f,-imaginary_greater_vs_time[integer_n])

        # Only use nonzero weights with corresponding frequency difference
        weights_greater_even = []
        w_values_greater_even = []
        weights_lesser_even = []
        w_values_lesser_even = []
        if all_weights_greater_even != []:
            w_values_greater_even, weights_greater_even = chop(all_w_values_greater, all_weights_greater_even)
        if all_weights_lesser_even != []:
            w_values_lesser_even, weights_lesser_even = chop(all_w_values_lesser, all_weights_lesser_even)

        # Normalize last for each k-value, enforce the SUM of lesser and greater weights being 1. 
        # This is the sum rule for a spectral function that depends on one k-value. 
        sum_of_lesser_and_greater_weights = sum(weights_lesser_even) + sum(weights_greater_even)
        for index,weight in enumerate(weights_lesser_even):
            weights_lesser_even[index] = weights_lesser_even[index] / sum_of_lesser_and_greater_weights
        for index,weight in enumerate(weights_greater_even):
            weights_greater_even[index] = weights_greater_even[index] / sum_of_lesser_and_greater_weights

        # Save weights and frequencies for each k-value
        greater_weights_for_all_k_n[integer_n] = weights_greater_even
        w_values_greater_for_all_k_n[integer_n] = w_values_greater_even
        lesser_weights_for_all_k_n[integer_n] = weights_lesser_even
        w_values_lesser_for_all_k_n[integer_n] = w_values_lesser_even


    return greater_weights_for_all_k_n, lesser_weights_for_all_k_n, w_values_greater_for_all_k_n, w_values_lesser_for_all_k_n

def plot_lesser_and_greater_information(problem,greater_weights_for_all_k_n, lesser_weights_for_all_k_n, w_values_greater_for_all_k_n, w_values_lesser_for_all_k_n):

    sub_figure,sub_axes = plt.subplots( len(problem.allowed_wavevectors), sharex=True)
    sub_figure.tight_layout()
    
    for sub_axis in sub_axes:
        sub_axis.spines['left'].set_linewidth(1)
        sub_axis.spines['right'].set_linewidth(1)
        sub_axis.spines['top'].set_linewidth(1)
        sub_axis.spines['bottom'].set_linewidth(1)
        sub_axis.xaxis.set_tick_params(which='major', size=10, width=1, direction='in', top='on', labelsize=5)
        sub_axis.xaxis.set_tick_params(which='minor', size=7, width=1, direction='in', top='on', labelsize=5)
        sub_axis.yaxis.set_tick_params(which='major', size=10, width=1, direction='in', right='on', labelsize=5)
        sub_axis.yaxis.set_tick_params(which='minor', size=7, width=1, direction='in', right='on', labelsize=5)

    for integer_n,k_n in enumerate(problem.allowed_wavevectors):
        # Lesser
        sub_axes[integer_n].scatter(w_values_lesser_for_all_k_n[integer_n], lesser_weights_for_all_k_n[integer_n], color="darkorange", marker='o',s=50)
        sub_axes[integer_n].vlines(w_values_lesser_for_all_k_n[integer_n], zeros(len(w_values_lesser_for_all_k_n[integer_n])),lesser_weights_for_all_k_n[integer_n],color="lightgrey",linestyle='dashed',linewidth=1)
        # Greater
        sub_axes[integer_n].scatter(w_values_greater_for_all_k_n[integer_n], greater_weights_for_all_k_n[integer_n], color="darkblue", marker='o',s=50)
        sub_axes[integer_n].vlines(w_values_greater_for_all_k_n[integer_n], zeros(len(w_values_greater_for_all_k_n[integer_n])), greater_weights_for_all_k_n[integer_n],color="lightgrey",linestyle='dashed',linewidth=1)
        
        sub_axes[integer_n].set_title(r"$k_n =$"+str(round(k_n,2)),size=5)
        sub_axes[integer_n].set_xlabel(r'$\omega$', size=10)

    # Hide x labels and tick labels for all but bottom plot.
    for sub_axis in sub_axes:
        sub_axis.label_outer()

    plt.savefig('./Figures_for_Thesis/approximate_lesser_and_greater_information_for_6_sites.png', dpi=800, transparent=False, bbox_inches='tight')

def real_green(integer_n, problem, greater_weights_for_all_k_n, lesser_weights_for_all_k_n, w_values_greater_for_all_k_n, w_values_lesser_for_all_k_n):
    Re_G = 0
    for g,weight in enumerate(greater_weights_for_all_k_n[integer_n]):
        w_g_values = w_values_greater_for_all_k_n[integer_n]
        Re_G = Re_G + weight / (w - w_g_values[g]) 
    for l,weight in enumerate(lesser_weights_for_all_k_n[integer_n]):
        w_l_values = w_values_lesser_for_all_k_n[integer_n]
        Re_G = Re_G + weight / (w - w_l_values[l]) 
    return Re_G 
    
def plot_real_green(problem, greater_weights_for_all_k_n, lesser_weights_for_all_k_n, w_values_greater_for_all_k_n, w_values_lesser_for_all_k_n):

    fig = plt.figure(figsize=(5,5))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.xaxis.set_tick_params(which='major', size=10, width=2, direction='in', top='on', labelsize=20)
    ax.xaxis.set_tick_params(which='minor', size=7, width=2, direction='in', top='on', labelsize=20)
    ax.yaxis.set_tick_params(which='major', size=10, width=2, direction='in', right='on', labelsize=20)
    ax.yaxis.set_tick_params(which='minor', size=7, width=2, direction='in', right='on', labelsize=20)
    ax.spines['left'].set_linewidth(2)
    ax.spines['right'].set_linewidth(2)
    ax.spines['top'].set_linewidth(2)
    ax.spines['bottom'].set_linewidth(2)

    w_values = linspace(-5,5,500)
    for integer_n,k_n in enumerate(problem.allowed_wavevectors):
        Re_G = real_green(integer_n, problem, greater_weights_for_all_k_n, lesser_weights_for_all_k_n, w_values_greater_for_all_k_n, w_values_lesser_for_all_k_n)
        Re_G_lamdified = lambdify(w, Re_G, modules=['numpy'])
        Re_G_values = Re_G_lamdified(w_values)
        my_linestyle = "solid"
        if integer_n == 4 or integer_n == 5:
            my_linestyle = "dashed"
        plt.plot( w_values, Re_G_values, label="$k_n=$"+str(round(k_n,2)), linestyle = my_linestyle)
    
    ax.set_xlabel(r"$\omega$",labelpad=10,size=20)
    plt.legend()
    plt.savefig("Figures_for_Thesis/approximate_green_real_part_for_6_sites.png",dpi=800, transparent=False, bbox_inches='tight')
    plt.show()

def band_structure(problem,k_n):
    hubbard_N = problem.hubbard_N
    band_structure = hubbard_N.band_structure(k_n).real
    return band_structure

def noninteracting_real_from_band_structure(problem,w,k_n):
    g_0_real = 1/(w - band_structure(problem,k_n) ) 
    return g_0_real

def denominator(problem,greater_weights_for_all_k_n, lesser_weights_for_all_k_n, w_values_greater_for_all_k_n, w_values_lesser_for_all_k_n ):
    denominators = {}
    for integer_n,k_n in enumerate(problem.allowed_wavevectors):
        object_for_denominator = Compressive_Sensing_Denominator(greater_weights_for_all_k_n[integer_n], lesser_weights_for_all_k_n[integer_n], w_values_greater_for_all_k_n[integer_n], w_values_lesser_for_all_k_n[integer_n])
        denominator = object_for_denominator.form()
        if denominator == 0:
            denominator = 1e-17*w # may be needed for numerical stability
        denominators[integer_n] = denominator
    return denominators
    
def numerator(problem,greater_weights_for_all_k_n, lesser_weights_for_all_k_n, w_values_greater_for_all_k_n, w_values_lesser_for_all_k_n ):
    denominators = denominator(problem,greater_weights_for_all_k_n, lesser_weights_for_all_k_n, w_values_greater_for_all_k_n, w_values_lesser_for_all_k_n )
    numerators = {}
    for integer_n,k_n in enumerate(problem.allowed_wavevectors):
        object_for_denominator = Compressive_Sensing_Denominator(greater_weights_for_all_k_n[integer_n], lesser_weights_for_all_k_n[integer_n], w_values_greater_for_all_k_n[integer_n], w_values_lesser_for_all_k_n[integer_n])
        band_structure_for_k_n = band_structure(problem,k_n)
        numerator = (w - band_structure_for_k_n) * denominators[integer_n] - object_for_denominator.product_over_l() * object_for_denominator.product_over_g()
        numerators[integer_n] = numerator
    return numerators
    
def real_self_energy(problem,numerators,denominators):
    real_self_energies = {}
    for integer_n,k_n in enumerate(problem.allowed_wavevectors):
        real_self_energy = numerators[integer_n] / denominators[integer_n]
        real_self_energies[integer_n] = real_self_energy
    return real_self_energies

def plot_real_self_energy(problem,real_self_energies):
    fig = plt.figure(figsize=(5,5))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.xaxis.set_tick_params(which='major', size=10, width=2, direction='in', top='on', labelsize=20)
    ax.xaxis.set_tick_params(which='minor', size=7, width=2, direction='in', top='on', labelsize=20)
    ax.yaxis.set_tick_params(which='major', size=10, width=2, direction='in', right='on', labelsize=20)
    ax.yaxis.set_tick_params(which='minor', size=7, width=2, direction='in', right='on', labelsize=20)
    ax.spines['left'].set_linewidth(2)
    ax.spines['right'].set_linewidth(2)
    ax.spines['top'].set_linewidth(2)
    ax.spines['bottom'].set_linewidth(2)

    w_values = linspace(-10,10,500)
    for integer_n,k_n in enumerate(problem.allowed_wavevectors):

        sigma = real_self_energies[integer_n] # symbolic expression 
        lamdified_self_energy = lambdify(w, sigma, modules=['numpy'])
        sigma_values = lamdified_self_energy(w_values)

        my_linestyle = "solid"
        if integer_n == 4 or integer_n == 5:
            my_linestyle = "dashed"
        plt.plot(w_values, sigma_values, label="$k_n=$"+str(round(k_n,2)), linestyle = my_linestyle)
    
    ax.set_xlabel(r"$\omega$",labelpad=10,size=20)
    ax.set_ylim([-10,10])
    plt.legend()
    plt.savefig("Figures_for_Thesis/approximate_real_self_energy_for_6_sites.png",dpi=800, transparent=False, bbox_inches='tight')
    plt.show()

def self_energy_poles(problem,denominators):
    denominator_roots_for_all_k_n = {}
    for integer_n,k_n in enumerate(problem.allowed_wavevectors):
        denominator = denominators[integer_n]
        denominator_roots = nroots(denominator,maxsteps=200000)
        denominator_roots_for_all_k_n[integer_n] = denominator_roots
    return denominator_roots_for_all_k_n

def self_energy_weights(problem,denominator_roots_for_all_k_n,numerators):
    
    self_energy_weights_for_all_k_n = {}
    for integer_n,k_n in enumerate(problem.allowed_wavevectors):

        weights = []
        denominator_roots = denominator_roots_for_all_k_n[integer_n]
        numerator = numerators[integer_n]

        # Form a new denominator based on the roots of the denominator
        new_denominator = 1
        for root in denominator_roots:
            if im(root) != 0:
                print("complex root")
                return
            new_denominator = new_denominator * (w - root)
        
        # Cancel the factor with this denominator_root by forming a new denominator without that factor
        for root in denominator_roots:
            denominator_expression = 1
            for den_root in denominator_roots:
                if den_root != root:
                    denominator_expression = denominator_expression * (w - den_root)
           
            # Get the 'denominator weight'
            denominator_weight = 1
            if denominator_expression != 1: # may be needed for numerical stability
                denominator_weight = simplify( denominator_expression.subs(w,root) )
            
            # Get the 'numerator weight'
            numerator_weight = simplify( numerator.subs(w,root) )

            # Store the self energy weight for that root
            weights.append(numerator_weight/denominator_weight)
        
        # Store the self energy weights for that k_n value
        self_energy_weights_for_all_k_n[integer_n] = weights

    return self_energy_weights_for_all_k_n

def plot_imaginary_self_energy(problem, denominator_roots_for_all_k_n, self_energy_weights_for_all_k_n):
    fig = plt.figure(figsize=(5,5))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.xaxis.set_tick_params(which='major', size=10, width=2, direction='in', top='on', labelsize=20)
    ax.xaxis.set_tick_params(which='minor', size=7, width=2, direction='in', top='on', labelsize=20)
    ax.yaxis.set_tick_params(which='major', size=10, width=2, direction='in', right='on', labelsize=20)
    ax.yaxis.set_tick_params(which='minor', size=7, width=2, direction='in', right='on', labelsize=20)
    ax.spines['left'].set_linewidth(2)
    ax.spines['right'].set_linewidth(2)
    ax.spines['top'].set_linewidth(2)
    ax.spines['bottom'].set_linewidth(2)

    for integer_n,k_n in enumerate(problem.allowed_wavevectors):

        poles = denominator_roots_for_all_k_n[integer_n]
        weights = self_energy_weights_for_all_k_n[integer_n]

        plt.scatter(poles, weights , marker='x',s=30,color='firebrick')
        plt.vlines(poles,zeros(len(poles)),weights,color='lightgrey', linewidth=1,linestyle='dashed',label="$k_n=$"+str(round(k_n,2)))
    
    ax.set_xlabel(r"$\omega$",labelpad=10,size=20)
    plt.legend()
    plt.savefig("Figures_for_Thesis/approximate_imaginary_self_energy_for_6_sites.png",dpi=800, transparent=False, bbox_inches='tight')
    plt.show()



def main():
    problem, U = define_problem()       

    # 0) Time Evolution
    imaginary_greater_vs_time, imaginary_lesser_vs_time, n, sampling_f = plot_time_evolution(problem)

    greater_weights_for_all_k_n, lesser_weights_for_all_k_n, w_values_greater_for_all_k_n, w_values_lesser_for_all_k_n = lesser_and_greater_information(problem, imaginary_greater_vs_time, imaginary_lesser_vs_time, n, sampling_f, U)
    
    # 1) Lesser and Greater Information
    plot_lesser_and_greater_information(problem,greater_weights_for_all_k_n, lesser_weights_for_all_k_n, w_values_greater_for_all_k_n, w_values_lesser_for_all_k_n)
    
    # 2) Real Part of Green Function
    plot_real_green(problem, greater_weights_for_all_k_n, lesser_weights_for_all_k_n, w_values_greater_for_all_k_n, w_values_lesser_for_all_k_n)
    
    denominators = denominator(problem,greater_weights_for_all_k_n, lesser_weights_for_all_k_n, w_values_greater_for_all_k_n, w_values_lesser_for_all_k_n )

    numerators = numerator(problem,greater_weights_for_all_k_n, lesser_weights_for_all_k_n, w_values_greater_for_all_k_n, w_values_lesser_for_all_k_n )

    real_self_energies = real_self_energy(problem,numerators,denominators)

    denominator_roots_for_all_k_n = self_energy_poles(problem,denominators)

    self_energy_weights_for_all_k_n = self_energy_weights(problem,denominator_roots_for_all_k_n,numerators)
    
    # 3) Real Part of Self Energy
    plot_real_self_energy(problem,real_self_energies)
    
    # 4) Imaginary Part of Self Energy
    plot_imaginary_self_energy(problem, denominator_roots_for_all_k_n, self_energy_weights_for_all_k_n)



main()



