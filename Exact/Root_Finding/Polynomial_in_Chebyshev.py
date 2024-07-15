from sympy import degree, lambdify, symbols
from numpy import cos, arccos, linspace, absolute, zeros, allclose, dot
from numpy.linalg import solve

class Polynomial_in_Chebyshev:
    def __init__(self,expression_in_t):
        self.expression_in_t = expression_in_t
        self.b_vector = self.calculate_b_vector()
        self.degree = degree(expression_in_t)
    def calculate_b_vector(self):
        t = symbols('t')
        p_bar = self.expression_in_t
        n = 0
        if p_bar == 0:
            pass
        else:
            n = degree(p_bar)
        p_bar_numerical_evaluation = lambdify(t,p_bar)
        t_i_values = linspace(-1,1,num=(n+1),endpoint=True)
        absolute_values = absolute(p_bar_numerical_evaluation(t_i_values))
        max_value = 0
        if n==0:
            if p_bar == 0:
                max_value = 1
            else:
                max_value = float(self.expression_in_t)
        else:
            max_value = max(absolute_values)
        y_bar_i_values = [] # vector for scaled y_bar_i values
        tau = zeros(  (n+1,n+1)  ) # matrix for Chebyshev polynomials
        for i,t_i in enumerate(t_i_values):
            y_bar_i_values.append(p_bar_numerical_evaluation(t_i)/max_value)
            for j in range(n+1):
                tau[i][j] = cos( j * arccos(t_i) )
        calculated_b_vector = solve(tau,y_bar_i_values) # solve tau b = y
        if allclose(dot(tau,calculated_b_vector), y_bar_i_values) == False:
            raise Exception("There is an issue with the linear algebra.")
        return calculated_b_vector
    def evaluate(self,t_value):
        sum = 0
        for j,b_j in enumerate(self.b_vector):
            sum = sum + b_j * cos( j * arccos( float(t_value) ) )
        return sum