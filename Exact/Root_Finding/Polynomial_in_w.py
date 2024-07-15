from sympy import Poly,degree,symbols
from numpy import flip

class Polynomial_in_w:
    def __init__(self,expression):
        self.expression = expression
        self.degree = degree(expression)
        self.coefficients = flip(Poly(expression).all_coeffs())
    def compute_symbolic_derivative_in_w(self):
        w = symbols('w',real=True)
        a = self.coefficients
        derivative_expression = 0
        for m in range(self.degree):
            derivative_expression = derivative_expression + a[m+1] * (m+1) * (w**m)
        return derivative_expression