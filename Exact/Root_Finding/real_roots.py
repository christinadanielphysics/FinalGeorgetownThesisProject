from Exact.Root_Finding.Polynomial_in_Chebyshev import *
from Exact.Root_Finding.Multiplicity_Computation import *
from Exact.Root_Finding.Root import *
from Exact.tolerances import *

from sympy import symbols,degree
from scipy.optimize import brenth


def search_interval(a,b,pbar_in_t,pbar_in_t_prime,is_derivative,polynomial_in_chebyshev,derivative_polynomial_in_chebyshev):
    
    if is_derivative == True:
        computed_root_value_in_t = brenth(derivative_polynomial_in_chebyshev.evaluate, -1, 1, maxiter=MAXITER_brent)
        if abs( polynomial_in_chebyshev.evaluate(computed_root_value_in_t) ) < TOLERANCE_brent_for_derivative:
            multiplicity_object = Multiplicity_Computation(pbar_in_t_prime,computed_root_value_in_t)
            computed_multiplicity = multiplicity_object.get_multiplicity()
            if computed_multiplicity != 2:
                print("Expected multiplicity is 2 but computed multiplicity is ",computed_multiplicity)
            else:
                print("Expected multiplicity is 2 and computed multiplicity is ",computed_multiplicity)
            t_root_with_multiplicity = Root(computed_root_value_in_t,computed_multiplicity,a,b)
            w_root_with_multiplicity = Root(t_root_with_multiplicity.convert_root_from_t_to_w(),computed_multiplicity,a,b)
            return w_root_with_multiplicity
        else:
            false_alarm = Root(None,0,a,b)
            return false_alarm
    else:
        computed_root_value_in_t = brenth(polynomial_in_chebyshev.evaluate, -1, 1, maxiter=MAXITER_brent)
        if abs( polynomial_in_chebyshev.evaluate(computed_root_value_in_t) ) < TOLERANCE_brent_for_polynomial:
            multiplicity_object = Multiplicity_Computation(pbar_in_t,computed_root_value_in_t)
            computed_multiplicity = multiplicity_object.get_multiplicity()
            if computed_multiplicity != 1:
                print("Expected multiplicity is 1 but computed multiplicity is ",computed_multiplicity)
            else:
                print("Expected multiplicity is 1 and computed multiplicity is ",computed_multiplicity)
            t_root_with_multiplicity = Root(computed_root_value_in_t,computed_multiplicity,a,b)
            w_root_with_multiplicity = Root(t_root_with_multiplicity.convert_root_from_t_to_w(),computed_multiplicity,a,b)
            return w_root_with_multiplicity
        else:
            false_alarm = Root(None,0,a,b)
            return false_alarm

def get_real_roots_with_multiplicities(problem,a_min,b_max,step,root_objects):

    print("looking for a root in the specified interval...")
    t = symbols('t')
    a = a_min
    b = a_min + step

    root_found = False
    while b <= b_max:
        print("a = ",a," and b =",b)
        w_in_terms_of_t = a + ((t+1)*(b-a)/2)
        pbar_in_t = problem.get_denominator(w_in_terms_of_t)
        pbar_in_t_prime = problem.get_derivative_of_denominator(w_in_terms_of_t)
        the_polynomial_in_chebyshev = Polynomial_in_Chebyshev(pbar_in_t)
        derivative_polynomial_in_chebyshev = Polynomial_in_Chebyshev(pbar_in_t_prime)
        if the_polynomial_in_chebyshev.evaluate(-1) * the_polynomial_in_chebyshev.evaluate(1) <= 0:
            w_root_with_multiplicity = search_interval(a,b,pbar_in_t,pbar_in_t_prime,False,the_polynomial_in_chebyshev,derivative_polynomial_in_chebyshev) 
            if w_root_with_multiplicity.root_value != None:
                root_objects.append(w_root_with_multiplicity)
                print("root_value:",w_root_with_multiplicity.root_value,"multiplicity",w_root_with_multiplicity.multiplicity)
                root_found = True
                break
        if derivative_polynomial_in_chebyshev.evaluate(-1) * derivative_polynomial_in_chebyshev.evaluate(1) <= 0:
            w_root_with_multiplicity = search_interval(a,b,pbar_in_t,pbar_in_t_prime,True,the_polynomial_in_chebyshev,derivative_polynomial_in_chebyshev)
            if w_root_with_multiplicity.root_value != None:
                root_objects.append(w_root_with_multiplicity)
                print("root_value:",w_root_with_multiplicity.root_value,"multiplicity",w_root_with_multiplicity.multiplicity)
                root_found = True
                break
        a = a + step
        b = b + step
    
    if root_found == True:
        print("Success!")
    else:
        print("Altering the interval slightly.")
        a_min = a_min - 0.005
        b_max = b_max + 0.005
        get_real_roots_with_multiplicities(problem,a_min,b_max,step,root_objects)
    return root_objects