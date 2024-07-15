from Exact.Problem import *
from Exact.Parameters import *
from Exact.Root_Finding.real_roots import get_real_roots_with_multiplicities
from Exact.tolerances import *
from Exact.Root_Finding.Root import *
from Exact.weight import *
from Exact.ExpectationValue import *

from sympy import plot,degree,nroots,symbols
from cmath import pi
import pandas as pd
import itertools
from scipy.integrate import quad
import scipy 

w = symbols('w',real=True)
w_n = symbols('w_n',real=True)

# Eight Atoms
# markers = itertools.cycle(('o', 'o', 'o', 'o', 'o', '.', '.', '.')) 
# colors = itertools.cycle(('blue','navy','royalblue','deepskyblue','silver','lightblue','lightgrey','silver','blue','navy','royalblue','deepskyblue','silver','lightblue','lightgrey','silver','blue','navy','royalblue','deepskyblue','silver','lightblue','lightgrey','silver','blue','navy','royalblue','deepskyblue','silver','lightblue','lightgrey','silver','blue','navy','royalblue','deepskyblue','silver','lightblue','lightgrey','silver','blue','navy','royalblue','deepskyblue','silver','lightblue','lightgrey','silver','blue','navy','royalblue','deepskyblue','silver','lightblue','lightgrey','silver'))
# sizes = itertools.cycle((80,80,80,80,80,40,40,40))
# w_values = linspace(-100,100,100)

# Six Atoms
# markers = itertools.cycle(('o', 'o', 'o', 'o', '.', '.', '.', '.')) 
# colors = itertools.cycle(('red','grey','salmon','firebrick','silver','coral','red','grey','salmon','firebrick','silver','coral','red','grey','salmon','firebrick','silver','coral','red','grey','salmon','firebrick','silver','coral','red','grey','salmon','firebrick','silver','coral','red','grey','salmon','firebrick','silver','coral','red','grey','salmon','firebrick','silver','coral','red','grey','salmon','firebrick','silver','coral','red','grey','salmon','firebrick','silver','coral','red','grey','salmon','firebrick','silver','coral'))
# sizes = itertools.cycle((80,80,80,80,40,40,40,40))
# w_values = linspace(-12,12,1000)

# Four Atoms
markers = itertools.cycle(('o', 'o', 'o', '.', '.', '.', '.', '.')) 
colors = itertools.cycle(('forestgreen','grey','yellowgreen','silver','forestgreen','grey','yellowgreen','silver','forestgreen','grey','yellowgreen','silver','forestgreen','grey','yellowgreen','silver','forestgreen','grey','yellowgreen','silver','forestgreen','grey','yellowgreen','silver','forestgreen','grey','yellowgreen','silver'))
sizes = itertools.cycle((80,80,80,40,40,40,40,40))
w_values = arange(-5,5,0.001)



def band_structure(problem,k_n):
    hubbard_N = problem.hubbard_N
    band_structure = hubbard_N.band_structure(k_n).real
    return band_structure

def plot_eigenvalues(eigenvalues):
    fig = plt.figure(figsize=(5,5))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.xaxis.set_tick_params(which='major', size=10, width=0, direction='in', top='off', labelsize=0)
    ax.xaxis.set_tick_params(which='minor', size=7, width=0, direction='in', top='off', labelsize=0)
    ax.yaxis.set_tick_params(which='major', size=10, width=2, direction='in', right='off', labelsize=20)
    ax.yaxis.set_tick_params(which='minor', size=7, width=2, direction='in', right='off', labelsize=20)
    ax.spines['left'].set_linewidth(2)
    ax.spines['right'].set_linewidth(0)
    ax.spines['top'].set_linewidth(0)
    ax.spines['bottom'].set_linewidth(0)
    
    visible_ticks = { "right": False}

    array = []
    for index,value in enumerate(eigenvalues):
        array.append(1)

    plt.tick_params(**visible_ticks)
    
    plt.hlines(xmin=zeros(len(eigenvalues)),xmax=array,y=eigenvalues,color='blue',linewidth=1)

    ax.set_ylabel(r'$E_m$', labelpad=10, size=20)
    ax.set_xlabel("")
    plt.savefig('./Figures_for_Thesis/eigenvalues_four.png', dpi=800, transparent=False, bbox_inches='tight')
    plt.show()

def plot_weights_for_each_k_n(problem,all_minus_eigenvalues,all_minus_weights,all_plus_eigenvalues,all_plus_weights):


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

    for n,k_n in enumerate(problem.allowed_wavevectors):
        col = next(colors)
        mark = next(markers)
        size = next(sizes)

        plt.vlines(all_minus_eigenvalues[n],zeros(len(all_minus_eigenvalues[n])),all_minus_weights[n],linewidth=1,color='lightgrey',linestyle='dashed')
        plt.vlines(all_plus_eigenvalues[n],zeros(len(all_plus_eigenvalues[n])),all_plus_weights[n],linewidth=1,color='lightgrey',linestyle='dashed')
        plt.scatter(all_minus_eigenvalues[n],all_minus_weights[n],marker=mark,s=size,label=r"$k_n =$"+str(round(k_n,2)),color=col)
        plt.scatter(all_plus_eigenvalues[n],all_plus_weights[n],marker=mark,s=size,label=r"$k_n =$"+str(round(k_n,2)),color=col)


    ax.set_xlabel(r'$\omega$', labelpad=10, size=20)
    plt.legend()
    plt.savefig('./Figures_for_Thesis/weights_for_each_k_n_for_4_sites_with_U=5_and_t=1.png', dpi=800, transparent=False, bbox_inches='tight')
    plt.show()   

def plot_lesser_and_greater(problem,all_minus_eigenvalues,all_minus_weights,all_plus_eigenvalues,all_plus_weights):


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

    for n,k_n in enumerate(problem.allowed_wavevectors):

        col = "orange"
        size = 40
        mark = 'o'
        plt.vlines(all_minus_eigenvalues[n],zeros(len(all_minus_eigenvalues[n])),all_minus_weights[n],linewidth=1,color=col,linestyle='dashed')
        plt.scatter(all_minus_eigenvalues[n],all_minus_weights[n],marker=mark,s=size,color=col)
        col = "blue"
        size = 40
        marker = 'o'
        plt.scatter(all_plus_eigenvalues[n],all_plus_weights[n],marker=mark,s=size,color=col)
        plt.vlines(all_plus_eigenvalues[n],zeros(len(all_plus_eigenvalues[n])),all_plus_weights[n],linewidth=1,color=col,linestyle='dashed')


    ax.set_xlabel(r'$\omega$', labelpad=10, size=20)
    plt.legend()
    plt.savefig('./Figures_for_Thesis/lesser_and_greater_for_4_sites_with_U=5_and_t=1.png', dpi=800, transparent=False, bbox_inches='tight')
    plt.show()   
       
def noninteracting_real_from_band_structure(problem,k_n):
    g_0_real = 1/(1j*w_n - band_structure(problem,k_n) ) 
    return g_0_real

def plot_real_noninteracting(problem):
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

    for n,k_n in enumerate(problem.allowed_wavevectors):
        col = next(colors)
        linestyle = 'solid'
        if n==4 or n==5:
            linestyle ='dashed'
        Re_G0 = noninteracting_real_from_band_structure(problem,k_n)
        Re_G0_lamdified = lambdify(w, Re_G0, modules=['numpy'])
        Re_G0_values = Re_G0_lamdified(w_values)
   
        plt.plot(w_values,Re_G0_values,label="$k_n =$"+str(round(k_n,2)),color=col,linestyle=linestyle,linewidth=2)

    ax.set_xlabel(r'$\omega$', labelpad=10, size=20)
    ax.set_ylim([-6,6])
    plt.legend()
    plt.savefig('./Figures_for_Thesis/real_noninteracting_for_4_sites_with_t=1.png', dpi=800, transparent=False, bbox_inches='tight')
    plt.show()  

def real_green(problem,minus_dimension, minus_eigenvalues, lesser_weights, plus_dimension, plus_eigenvalues, greater_weights):
    Re_G = 0
    for g in range(plus_dimension):
        Re_G = Re_G + greater_weights[g].real / (1j*w_n - plus_eigenvalues[g]) 
    for l in range(minus_dimension):
        Re_G = Re_G + lesser_weights[l].real / (1j*w_n - minus_eigenvalues[l]) 
    return Re_G

def plot_real_interacting(problem, minus_dimension, minus_eigenvalues, lesser_weights,plus_dimension, plus_eigenvalues, greater_weights):
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

    for n,k_n in enumerate(problem.allowed_wavevectors):
        col = next(colors)
        linestyle = 'solid'
        if n==4 or n==5:
            linestyle ='dashed' 
            
        Re_G = real_green(problem, minus_dimension[n], minus_eigenvalues[n], lesser_weights[n], plus_dimension[n], plus_eigenvalues[n], greater_weights[n]) 
        Re_G_lamdified = lambdify(w, Re_G, modules=['numpy'])
        Re_G_values = Re_G_lamdified(w_values)
        
        plt.plot(w_values,Re_G_values,label="$k_n =$"+str(round(k_n,2)),color=col,linestyle=linestyle,linewidth=2)


    ax.set_xlabel(r'$\omega$', labelpad=10, size=20)
    ax.set_ylim([-2,2])
    plt.legend()
    plt.savefig('./Figures_for_Thesis/real_interacting_for_4_sites_with_t=1_and_U=5.png', dpi=800, transparent=False, bbox_inches='tight')
    plt.show()  

    return Re_G_values

def plot_real_self_energy(problem,numerators,denominators):
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

    for n,k_n in enumerate(problem.allowed_wavevectors):
        col = next(colors)
        linestyle = 'solid'
        if n==3:
            linestyle ='dashed' 

        sigma_lamdified = lambdify(w, numerators[n]/denominators[n], modules=['numpy'])
        sigma_values = sigma_lamdified(w_values)
        plt.plot(w_values, sigma_values, label="$k_n =$"+str(round(k_n,2)),c=col,linestyle=linestyle,linewidth=2)

    ax.set_xlabel(r'$\omega$', labelpad=10, size=20)
    ax.set_ylim([-1,1])
    plt.legend()
    plt.savefig('./Figures_for_Thesis/real_self_energy_for_4_sites_with_t=1_and_U=5.png', dpi=800, transparent=False, bbox_inches='tight')
    plt.show()  
    return sigma_values

def plot_imaginary_self_energy(problem,numerators,denominators):
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


    all_weights = {}
    all_denominator_roots = {}

    for n,k_n in enumerate(problem.allowed_wavevectors):
        denominator_roots = nroots(denominators[n],maxsteps=200000) 
        number_of_self_energy_poles = len(denominator_roots)
        print("NUMBER OF SELF ENERGY POLES:",number_of_self_energy_poles)
        w = symbols('w')
        weights = []
        for denominator_root in denominator_roots:
            if im(denominator_root) != 0:
                print("complex root!")
                return
            numerator_lamdified = lambdify(w,numerators[n],'numpy')
            numerator_value = numerator_lamdified(denominator_root)
            expression = 1
            for den_root in denominator_roots:
                if den_root != denominator_root:
                    expression = expression * (w - den_root) # cancel the factor with this denominator_root by forming a new denominator without that factor
            denominator_value = 1
            if expression != 1:
                denominator_value = simplify(expression.subs(w,denominator_root))
            weight = numerator_value/denominator_value
            weights.append(weight)
            plt.scatter(denominator_root,weight,s=30,color='blue',marker='x')
        all_weights[n] = weights
        all_denominator_roots[n] = denominator_roots
        plt.vlines(denominator_roots,zeros(number_of_self_energy_poles),weights,color='lightgrey',linestyles='dashed',linewidth=1,label="$k_n =$"+str(round(k_n,2)))


    ax.set_xlabel(r'$\omega$', labelpad=10, size=20)
    plt.legend()
    plt.savefig('./Figures_for_Thesis/imaginary_self_energy_for_4_sites_with_t=1_and_U=5.png', dpi=800, transparent=False, bbox_inches='tight')
    plt.show()  

    return all_denominator_roots, all_weights



def exact(U,t_1):
    mu = U/2
    N_up = 2
    N_down = 2
    operator_spin = "up"
    t_0 = 0
    t_2 = 0
    V = 4


    parameters = Parameters(t_0,t_1,t_2,U,mu,V)
    problem = Problem(parameters,operator_spin,N_up,N_down)
    w_0 = problem.hubbard_N.ground_state_energy

    all_numerators = {}
    all_denominators = {}
    all_minus_eigenvalues = {}
    all_plus_eigenvalues = {}
    all_minus_weights = {}
    all_plus_weights = {}
    all_minus_dimensions = {}
    all_plus_dimensions = {}
    for n,k_n in enumerate(problem.allowed_wavevectors):
        numerator, denominator, derivative = problem.form_numerator_and_denominator_and_derivative(k_n)
        problem.set_numerator_and_denominator_and_derivative(numerator,denominator,derivative)
        all_numerators[n] = numerator
        all_denominators[n] = denominator

        new_minus_dimension, new_minus_eigenvalues, new_lesser_weights_for_k_n = problem.denominator.filter_lesser(k_n)
        new_plus_dimension, new_plus_eigenvalues, new_greater_weights_for_k_n = problem.denominator.filter_greater(k_n)
        for index,value in enumerate(new_minus_eigenvalues): 
            new_minus_eigenvalues[index] = w_0 - new_minus_eigenvalues[index] 
        for index,value in enumerate(new_plus_eigenvalues):
            new_plus_eigenvalues[index] = new_plus_eigenvalues[index] - w_0
        all_minus_eigenvalues[n] = new_minus_eigenvalues
        all_plus_eigenvalues[n] = new_plus_eigenvalues
        all_minus_weights[n] = new_lesser_weights_for_k_n
        all_plus_weights[n] = new_greater_weights_for_k_n
        all_minus_dimensions[n] = new_minus_dimension
        all_plus_dimensions[n] = new_plus_dimension

    # # 0) Eigenvalues
    #plot_eigenvalues(problem.hubbard_N.eigenvalues)

    # # 1) Lesser and Greater Weights for different k_n values
    #plot_weights_for_each_k_n(problem,all_minus_eigenvalues,all_minus_weights,all_plus_eigenvalues,all_plus_weights)

    # # 2) Real Part of Noninteracting Green's Function
    #plot_real_noninteracting(problem)

    # # 3) Real Part of Interacting Green's function
    #real_g = plot_real_interacting(problem, all_minus_dimensions, all_minus_eigenvalues, all_minus_weights, all_plus_dimensions, all_plus_eigenvalues, all_plus_weights)

    # # 4) Plot real part of self energy
    #real_sigma = plot_real_self_energy(problem,all_numerators,all_denominators)

    # # 5) Compute poles and weights of Self Energy
    # all_denominator_roots, all_weights = plot_imaginary_self_energy(problem,all_numerators,all_denominators)

    

    
    U_expectation = 0
    for n,k_n in enumerate(problem.allowed_wavevectors):

        G = real_green(problem, all_minus_dimensions[n], all_minus_eigenvalues[n], all_minus_weights[n], all_plus_dimensions[n], all_plus_eigenvalues[n], all_plus_weights[n]) 
        G0 = noninteracting_real_from_band_structure(problem,k_n)
        sigma = (1/G0 - 1/G)
        

        expression_1 = re(G)* re(sigma) * 1/(2*pi)
        expression_2 = -im(G) * im(sigma) * 1/(2*pi)
        expression_3 = im(G) * re(sigma) * 1/(2*pi)
        expression_4 = re(G) * im(sigma) * 1/(2*pi)

        expression = expression_1 + expression_2 + expression_3 + expression_4


        expression_lamdified = lambdify(w_n,expression,modules='scipy')

        result = quad(lambda omega_n: expression_lamdified(omega_n),-numpy.inf,numpy.inf)

        U_expectation = U_expectation + result[0] 








    # # 9) Kinetic Energy
    ke_sum_1 = 0
    for n,k_n in enumerate(problem.allowed_wavevectors):
        g_weights = all_plus_weights[n]
        for g_weight in g_weights:
            ke_sum_1 = ke_sum_1 + g_weight * band_structure(problem,k_n)
    
    ke_sum_2 = 0
    for n,k_n in enumerate(problem.allowed_wavevectors):
        l_weights = all_minus_weights[n]
        for l_weight in l_weights:
            ke_sum_2 = ke_sum_2 + l_weight * band_structure(problem,k_n)


    print("U_expectation",2*U_expectation) 
    print("computed ground state energy", mu*N_down-mu*(N_up+N_down)-(-(re( 2*ke_sum_2 + U_expectation)  ) ))
    print("exact ground state energy",w_0)
    print("-mu * N = ",-mu*(N_up + N_down))

    # vary U --- plt.scatter(U,-(-(re( 2*ke_sum_2 + U_expectation*U) + U*N_down ) - 2*(N_up+N_down)),c='red',marker ='x')
    plt.scatter(U,mu*N_down-mu*(N_up+N_down)-(-(re( 2*ke_sum_2 + U_expectation)  ) ),c='red',marker ='x',s=20)

            




  
        

        
