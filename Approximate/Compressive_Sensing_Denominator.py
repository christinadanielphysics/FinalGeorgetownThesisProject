from sympy import symbols

class Compressive_Sensing_Denominator:
    def __init__(self,greater_weights_for_k_n,lesser_weights_for_k_n,w_values_greater,w_values_lesser):
        self.greater_weights_for_k_n = greater_weights_for_k_n
        self.lesser_weights_for_k_n = lesser_weights_for_k_n
        self.w_values_greater = w_values_greater # w_g - w_0 
        self.w_values_lesser = w_values_lesser # w_0 - w_l 
    def product_over_l(self):
        w = symbols('w',real=True)
        product_over_l = 1
        for value in self.w_values_lesser:
            W_l = -value
            product_over_l = product_over_l * (w + W_l)
        return product_over_l
    def product_over_h_not_g(self,g):
        w = symbols('w',real=True)
        product_over_h_not_g = 1
        for h,value in enumerate(self.w_values_greater):
            if h != g:
                W_h = -value
                product_over_h_not_g = product_over_h_not_g * (w + W_h)
        return product_over_h_not_g
    def sum_over_g(self):
        sum_over_g = 0
        for g,c_g in enumerate(self.greater_weights_for_k_n):
            if c_g != 0:
                sum_over_g = sum_over_g + c_g * self.product_over_h_not_g(g)
        return sum_over_g
    def product_over_g(self):
        w = symbols('w',real=True)
        product_over_g = 1
        for value in self.w_values_greater:
            W_g = -value
            product_over_g = product_over_g * (w + W_g)
        return product_over_g
    def product_over_m_not_l(self,l):
        w = symbols('w',real=True)
        product_over_m_not_l = 1
        for m,value in enumerate(self.w_values_lesser):
            if m != l:
                W_m = -value
                product_over_m_not_l = product_over_m_not_l * (w + W_m)
        return product_over_m_not_l
    def sum_over_l(self):
        sum_over_l = 0
        for l,c_l in enumerate(self.lesser_weights_for_k_n):
            if c_l != 0:
                sum_over_l = sum_over_l + c_l * self.product_over_m_not_l(l)
        return sum_over_l
    def form(self):
        term_1 = self.product_over_l() * self.sum_over_g()
        term_2 = self.product_over_g() * self.sum_over_l()
        denominator = term_1 + term_2
        return denominator
