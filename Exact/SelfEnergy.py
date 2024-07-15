from .Denominator import *

class SelfEnergy:
    def __init__(self,green):
        self.green = green
    def constant_term_up(self,k_n):
        w_0 = self.green.ground_state_energy
        band_structure = self.green.hubbard_N.band_structure(k_n)
        
        greater_piece = 0
        for g in range(self.green.plus_dimension):
            greater_piece = greater_piece + self.green.greater_weight_up_in_k_space(g, k_n).real * ( self.green.plus_eigenvalues[g] - w_0 - band_structure )
        
        lesser_piece = 0
        for l in range(self.green.minus_dimension):
            lesser_piece = lesser_piece + self.green.lesser_weight_up_in_k_space(l, k_n).real * ( w_0 - self.green.minus_eigenvalues[l] - band_structure )
        
        return greater_piece + lesser_piece
    def constant_term_down(self,k_n):
        w_0 = self.green.ground_state_energy
        band_structure = self.green.hubbard_N.band_structure(k_n)

        greater_piece = 0
        for g in range(self.green.plus_dimension):
            greater_piece = greater_piece + self.green.greater_weight_down_in_k_space(g, k_n).real * ( self.green.plus_eigenvalues[g] - w_0 - band_structure )
        
        lesser_piece = 0
        for l in range(self.green.minus_dimension):
            lesser_piece = lesser_piece + self.green.lesser_weight_down_in_k_space(l, k_n).real * ( w_0 - self.green.minus_eigenvalues[l] - band_structure )
        
        return greater_piece + lesser_piece

    
    

