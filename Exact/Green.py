from sympy import symbols
from numpy import conjugate
from cmath import exp

# The column v[:, i] is the normalized eigenvector corresponding to the eigenvalue w[i]

class Green:
    def __init__(self,problem):
        self.problem = problem
        self.storage = problem.storage

        self.V = problem.parameters.V
        self.allowed_wavevectors = problem.allowed_wavevectors
        self.hubbard_N = problem.hubbard_N
        self.ground_state_energy = problem.hubbard_N.ground_state_energy

        self.mu = problem.parameters.U / 2 # new: using canonical ensemble




        self.plus_eigenvectors = problem.hubbard_plus.eigenvectors
        self.minus_eigenvectors = problem.hubbard_minus.eigenvectors

        self.plus_eigenvalues = problem.hubbard_plus.eigenvalues
        self.minus_eigenvalues = problem.hubbard_minus.eigenvalues

        self.minus_dimension = problem.hubbard_minus.dimension
        self.plus_dimension = problem.hubbard_plus.dimension





    def greater_weight_up(self,i,j,g):
        eigenvector_plus = self.plus_eigenvectors[:,g]
        matrix_element_1 = self.storage.plus_c_dagger_up_ground(eigenvector_plus,j)
        matrix_element_2 = conjugate( self.storage.plus_c_dagger_up_ground(eigenvector_plus,i) )
        return matrix_element_1 * matrix_element_2
    def greater_weight_down(self,i,j,g):
        eigenvector_plus = self.plus_eigenvectors[:,g]
        matrix_element_1 = self.storage.plus_c_dagger_down_ground(eigenvector_plus,j)
        matrix_element_2 = conjugate( self.storage.plus_c_dagger_down_ground(eigenvector_plus,i) )
        return matrix_element_1 * matrix_element_2
    def lesser_weight_up(self,i,j,l):
        eigenvector_minus = self.minus_eigenvectors[:,l]
        matrix_element_1 = self.storage.minus_c_up_ground(eigenvector_minus,i)
        matrix_element_2 = conjugate( self.storage.minus_c_up_ground(eigenvector_minus,j) )
        return matrix_element_1 * matrix_element_2
    def lesser_weight_down(self,i,j,l):
        eigenvector_minus = self.minus_eigenvectors[:,l]
        matrix_element_1 = self.storage.minus_c_down_ground(eigenvector_minus,i)
        matrix_element_2 = conjugate( self.storage.minus_c_down_ground(eigenvector_minus,j) )
        return matrix_element_1 * matrix_element_2
    def spectral_function_up(self,i,j):
        greater_angular_frequency_differences = []
        greater_weights = []
        lesser_angular_frequency_differences = []
        lesser_weights = []
        for g in range(self.plus_dimension):
            greater_angular_frequency_differences.append(self.plus_eigenvalues[g] - self.ground_state_energy)
            greater_weights.append( self.greater_weight_up(i,j,g) )
        for l in range(self.minus_dimension):
            lesser_angular_frequency_differences.append(self.ground_state_energy - self.minus_eigenvalues[l])
            lesser_weights.append( self.lesser_weight_up(i,j,l) )
        return greater_angular_frequency_differences, greater_weights, lesser_angular_frequency_differences, lesser_weights
    def spectral_function_down(self,i,j):
        greater_angular_frequency_differences = []
        greater_weights = []
        lesser_angular_frequency_differences = []
        lesser_weights = []
        for g in range(self.plus_dimension):
            greater_angular_frequency_differences.append(self.plus_eigenvalues[g] - self.ground_state_energy)
            greater_weights.append( self.greater_weight_down(i,j,g) )    
        for l in range(self.minus_dimension):
            lesser_angular_frequency_differences.append(self.ground_state_energy - self.minus_eigenvalues[l])
            lesser_weights.append( self.lesser_weight_down(i,j,l) )  
        return greater_angular_frequency_differences, greater_weights, lesser_angular_frequency_differences, lesser_weights
    def greater_weight_up_in_k_space(self,g,k_n):
        weight = 0
        for i in range(self.V):
            for j in range(self.V):
                weight = weight + (1/self.V) * exp(1j * k_n * (i - j) ) * self.greater_weight_up(i,j,g)
        return weight
    def greater_weight_down_in_k_space(self,g,k_n):
        weight = 0
        for i in range(self.V):
            for j in range(self.V):
                weight = weight + (1/self.V) * exp(1j * k_n * (i - j) ) * self.greater_weight_down(i,j,g)
        return weight
    def lesser_weight_up_in_k_space(self,l,k_n):
        weight = 0
        for i in range(self.V):
            for j in range(self.V):
                weight = weight + (1/self.V) * exp(1j * k_n * (i - j) ) * self.lesser_weight_up(i,j,l)
        return weight
    def lesser_weight_down_in_k_space(self,l,k_n):
        weight = 0
        for i in range(self.V):
            for j in range(self.V):
                weight = weight + (1/self.V) * exp(1j * k_n * (i - j) ) * self.lesser_weight_down(i,j,l)
        return weight
    def noninteracting_from_bandstructure(self,k_n):
        w = symbols('w',real=True)
        epsilon_of_k_n = self.hubbard_N.band_structure(k_n)
        G_0_of_k_n_sigma = 1/(w - epsilon_of_k_n)
        return G_0_of_k_n_sigma
    def real_part_in_k_space_up(self,k_n):
        w = symbols('w',real=True)
        w_0 = self.ground_state_energy
        Re_G = 0
        for g in range(self.plus_dimension):
            Re_G = Re_G + self.greater_weight_up_in_k_space(g, k_n).real / (w - self.plus_eigenvalues[g] + w_0)
        for l in range(self.minus_dimension):
            Re_G = Re_G + self.lesser_weight_up_in_k_space(l, k_n).real / (w + self.minus_eigenvalues[l] - w_0)
        return Re_G
    def real_part_in_k_space_down(self, k_n):
        w = symbols('w',real=True)
        w_0 = self.ground_state_energy
        real_part = 0
        for g in range(self.plus_dimension):
            Re_G = Re_G + self.greater_weight_down_in_k_space(g, k_n).real / (w - self.plus_eigenvalues[g] + w_0)
        for l in range(self.minus_dimension):
            Re_G = Re_G + self.lesser_weight_down_in_k_space(l, k_n).real / (w + self.minus_eigenvalues[l] - w_0)
        return Re_G
    def spectral_function_up_for_one_wavevector(self, k_n):
        greater_angular_frequency_differences = []
        greater_weights = []
        lesser_angular_frequency_differences = []
        lesser_weights = []
        for g in range(self.plus_dimension):
            greater_angular_frequency_differences.append(self.plus_eigenvalues[g] - self.ground_state_energy)
            greater_weights.append( self.greater_weight_up_in_k_space(g,k_n) )
        for l in range(self.minus_dimension):
            lesser_angular_frequency_differences.append(self.ground_state_energy - self.minus_eigenvalues[l])
            lesser_weights.append( self.lesser_weight_up_in_k_space(l,k_n) )
        return greater_angular_frequency_differences, greater_weights, lesser_angular_frequency_differences, lesser_weights
    def spectral_function_down_for_one_wavevector(self, k_n):
        greater_angular_frequency_differences = []
        greater_weights = []
        lesser_angular_frequency_differences = []
        lesser_weights = []
        for g in range(self.plus_dimension):
            greater_angular_frequency_differences.append(self.plus_eigenvalues[g] - self.ground_state_energy)
            greater_weights.append( self.greater_weight_down_in_k_space(g,k_n) )
        for l in range(self.minus_dimension):
            lesser_angular_frequency_differences.append(self.ground_state_energy - self.minus_eigenvalues[l])
            lesser_weights.append( self.lesser_weight_down_in_k_space(l,k_n) )
        return greater_angular_frequency_differences, greater_weights, lesser_angular_frequency_differences, lesser_weights

