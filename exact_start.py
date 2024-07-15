from Exact.Problem import *
from Exact.Parameters import *
from Exact.Root_Finding.real_roots import get_real_roots_with_multiplicities
from Exact.tolerances import *
from Exact.Root_Finding.Root import *
from Exact.weight import *
from Exact.ExpectationValue import *

from sympy import plot,degree,nroots,symbols
from cmath import pi

def get_user_input(numerator,denominator):
    print("Time to look for roots! Pay attention to lower and upper bounds of a new root. The more precise the better.")
    plot(numerator/denominator)
    done = int(input("Enter 0 if there is a new root to identify. Otherwise, enter 1\n"))
    if done == 1:
        return None,None,1
    else:
        print("Testing your memory... Where was the root?")
        a_min = float(input("Enter a lower bound fairly close to the root "))
        b_max = float(input("Enter an upper bound fairly close to the root "))
        
        return a_min,b_max,done

def display_root_and_weight_information(root_objects,weight):
    print("root_value:",root_objects[-1].root_value)
    print("multiplicity:",root_objects[-1].multiplicity)
    print("weight:",weight)

def exact(U):
    N_up = 2
    N_down = 2
    operator_spin = "up"
    t_0 = 0
    t_1 = 1
    t_2 = 0
    mu = U/2 
    V = 4
    parameters = Parameters(t_0,t_1,t_2,U,mu,V)
    problem = Problem(parameters,operator_spin,N_up,N_down)

    i = 0
    j = 0

    # problem.plot_spectral_function(i,j)
    # problem.hubbard_N.plot_band_structure()    
    # problem.plot_noninteracting_based_on_band_structure()
    # problem.plot_green_real_part_in_k_space()
    # problem.print_constant()
    
    #problem.plot_spectral_function_for_wavevectors()

    
    all_numerators = {}
    all_denominators = {}
    root_objects_for_all_n = {}
    weights_for_all_n = {}
    for n,k_n in enumerate(problem.allowed_wavevectors):

        print("")

        print("WAVEVECTOR =",k_n)

        numerator, denominator, derivative = problem.form_numerator_and_denominator_and_derivative(k_n)
        problem.set_numerator_and_denominator_and_derivative(numerator,denominator,derivative)
        all_numerators[n] = numerator
        all_denominators[n] = denominator

        denominator_roots = nroots(denominator,maxsteps=200000) 
        print("NUMBER OF SELF ENERGY POLES:",len(denominator_roots))

        w = symbols('w')

        sum_of_weights = 0
        weights = []

        for denominator_root in denominator_roots:
            numerator_lamdified = lambdify(w,numerator,'numpy')
            numerator_value = numerator_lamdified(denominator_root)
            
            # cancel the factor with this denominator_root by forming a new denominator without that factor
            expression = 1
            for den_root in denominator_roots:
                if den_root != denominator_root:
                    expression = expression * (w - den_root)

            denominator_value = 1
            if expression != 1:
                denominator_value = simplify(expression.subs(w,denominator_root))


            weight = numerator_value/denominator_value
            weights.append(weight)
            sum_of_weights = sum_of_weights + weight
        # print("SUM OF WEIGHTS",sum_of_weights)
        root_objects_for_all_n[n] = denominator_roots
        weights_for_all_n[n] = weights

        print("")
    
    
    denominator_object_for_potential = Denominator(problem.green)
    U_expectation_value = denominator_object_for_potential.potential(root_objects_for_all_n,weights_for_all_n,all_numerators,all_denominators)


    expectation_value = ExpectationValue(root_objects_for_all_n,weights_for_all_n,problem,all_numerators,all_denominators)
    T_expectation_value = re(expectation_value.kinetic())
    
    N_expectation_value = mu * (N_up + N_down)


    print("\n RESULTS: \n")
    print("<T> =",T_expectation_value)
    print("<U> =",U_expectation_value)
    print("<H> =",T_expectation_value + U_expectation_value - N_expectation_value)
    print("ground state energy from ED =", problem.hubbard_N.ground_state_energy )


    plt.figure(1)
    plt.scatter(U,U_expectation_value,color='red',label='calculated')
    plt.scatter(U,problem.hubbard_N.ground_state_energy - T_expectation_value,color='blue',marker='x',label='exact')
    plt.scatter(U,U_expectation_value/(problem.hubbard_N.ground_state_energy - T_expectation_value ),label='ratio',marker='.',color='green')


step = 3
U = 1
for index in range(10):
    exact(U+0.5*index)

plt.legend()
plt.show()
