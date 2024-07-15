from sympy import lambdify,symbols,simplify,div,rem,diff

class Multiplicity_Computation:
    def __init__(self,symbolic_f,root_value):
        self.symbolic_f = symbolic_f
        self.root_value = root_value
    def get_multiplicity(self):
        t = symbols('t')
        divisor = (t - self.root_value)
        expression = self.symbolic_f 

        multiplicity = 0
        while True:
            quotient,remainder = div(expression,divisor)
            if diff(remainder) != 0:
                break
            expression = quotient + remainder/divisor
            multiplicity = multiplicity + 1

        return multiplicity

