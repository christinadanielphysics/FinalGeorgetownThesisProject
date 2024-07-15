from Exact.tolerances import *

import matplotlib.pyplot as plt
from numpy import linspace
from sympy import symbols,simplify


# used
def compute_weight_symbolically(root_object,problem):
    w = symbols('w',real=True)
    numerical_root = root_object.root_value
    multiplicity = root_object.multiplicity

    expression = ( (w - numerical_root)**multiplicity ) * problem.get_numerator(w) / problem.get_denominator(w)
    simplified_expression = simplify(expression)
    weight = simplified_expression.subs(w,numerical_root)
    return weight


