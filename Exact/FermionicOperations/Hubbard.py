from .Basis import *

from scipy import linalg as LA
from cmath import pi,cos
import matplotlib.pyplot as plt
from numpy import arange

class Hubbard:
    def __init__(self,N_up,N_down,parameters):
        """
        Completed on 12/22/22
        """

        self.N_up = N_up
        self.N_down = N_down
        self.V = parameters.V
        self.parameters = parameters

        self.ground_state_energy = None
        self.ground_state_vector = None
        self.eigenvalues = None
        self.eigenvectors = None

        self.basis_object = Basis(self.V,N_up,N_down)
        self.basis_object.form()
        self.basis = self.basis_object.basis
        self.dimension = len(self.basis)

        self.hopping_matrix = zeros((self.V,self.V))
        self.kinetic_matrix_up = zeros((self.dimension, self.dimension))
        self.kinetic_matrix_down = zeros((self.dimension, self.dimension))
        self.kinetic_matrix = zeros((self.dimension, self.dimension))
        self.interaction_matrix = zeros((self.dimension, self.dimension))
        self.number_matrix = zeros((self.dimension, self.dimension))
        self.hamiltonian_matrix = zeros((self.dimension, self.dimension))
    def form_hopping_matrix(self):
        """
        Completed on 12/22/22
        """
        for i in range(self.V):
            for j in range(self.V):
                if abs(i-j) == 0:
                    self.hopping_matrix[i][j] = self.parameters.t_0
                elif abs(i-j) == 1 or abs(i-j) == self.V-1:
                    self.hopping_matrix[i][j] = self.parameters.t_1
                elif abs(i-j) == 2 or abs(i-j) == self.V-2:
                    self.hopping_matrix[i][j] = self.parameters.t_2
                else:
                    self.hopping_matrix[i][j] = 0
    def form_up_kinetic_matrix(self):
        """
        Completed on 12/22/22
        """
        for col_index,ket in enumerate(self.basis):
            for j in range(self.V):
                result_1 = ket.apply_up_annihilation(j)
                if result_1.coefficient != 0:
                    for i in range(self.V):
                        result_2 = result_1.apply_up_creation(i)
                        if result_2.coefficient != 0:
                            resulting_index = self.basis_object.get_index(result_2)
                            new_coefficient = result_2.coefficient * self.hopping_matrix[i][j]
                            self.kinetic_matrix_up[resulting_index][col_index] = self.kinetic_matrix_up[resulting_index][col_index] + new_coefficient
    def form_down_kinetic_matrix(self):
        """
        Completed on 12/22/22
        """
        for col_index,ket in enumerate(self.basis):
            for j in range(self.V):
                result_1 = ket.apply_down_annihilation(j)
                if result_1.coefficient != 0:
                    for i in range(self.V):
                        result_2 = result_1.apply_down_creation(i)
                        if result_2.coefficient != 0:
                            resulting_index = self.basis_object.get_index(result_2)
                            new_coefficient = result_2.coefficient * self.hopping_matrix[i][j]
                            self.kinetic_matrix_down[resulting_index][col_index] = self.kinetic_matrix_down[resulting_index][col_index] + new_coefficient
    def form_kinetic_matrix(self):
        """
        Completed on 12/22/22
        """
        self.form_hopping_matrix()
        self.form_up_kinetic_matrix()   
        self.form_down_kinetic_matrix()
        self.kinetic_matrix = self.kinetic_matrix_up + self.kinetic_matrix_down
    def form_interaction_matrix(self):
        """
        Completed on 12/22/22
        """
        for col_index,ket in enumerate(self.basis):
            for l in range(self.V):
                result_1 = ket.apply_down_annihilation(l)
                if result_1.coefficient != 0:
                    result_2 = result_1.apply_down_creation(l)
                    if result_2.coefficient != 0:
                        result_3 = result_2.apply_up_annihilation(l)
                        if result_3.coefficient != 0:
                            result_4 = result_3.apply_up_creation(l)
                            if result_4.coefficient != 0:
                                resulting_index = self.basis_object.get_index(result_4)
                                new_coefficient = result_4.coefficient * self.parameters.U
                                self.interaction_matrix[resulting_index][col_index] = self.interaction_matrix[resulting_index][col_index] + new_coefficient
    def form_number_matrix(self):
        """
        Completed on 12/22/22
        """
        for row_index in range(self.dimension):
            self.number_matrix[row_index][row_index] = (-1) * self.parameters.mu * (self.N_up + self.N_down)
    def form_hamiltonian_matrix(self):
        """
        Completed on 12/22/22
        """
        self.form_kinetic_matrix()
        self.form_interaction_matrix()
        self.form_number_matrix()
        self.hamiltonian_matrix = self.kinetic_matrix + self.interaction_matrix + self.number_matrix
    def compute_eigenvalues_and_eigenvectors(self):
        """
        Completed on 12/22/22
        """
        w,v = LA.eigh(self.hamiltonian_matrix)
        self.eigenvalues = w
        self.eigenvectors = v
        self.ground_state_energy = w[0]
        self.ground_state_vector = v[:, 0] # The column v[:, i] is the normalized eigenvector corresponding to the eigenvalue w[i]
    def allowed_wavevectors(self):
        allowed_wavevectors = []
        for n in range(self.V):
            k_n = 2 * pi * n / self.V
            allowed_wavevectors.append(k_n)
        return allowed_wavevectors
    def band_structure(self,k_n):
        band_structure = self.parameters.t_1 * 2 * cos( k_n )
        return band_structure
    def plot_band_structure(self):
        
        allowed_wavevectors = self.allowed_wavevectors()
        band_structure = []
        for k_n in allowed_wavevectors:
            band_structure.append( self.band_structure(k_n) )

        fine_wavevectors = arange(0,2*pi,0.01)
        fine_points = []
        for k_n in fine_wavevectors:
            fine_points.append( self.band_structure(k_n) )

        plt.figure(1)
        plt.scatter(allowed_wavevectors, band_structure ,color='black')
        plt.plot(fine_wavevectors,fine_points,color='black')
        plt.xlim([0,2 * pi])
        plt.xlabel(r"$k_n$")
        plt.ylabel(r"$\varepsilon(k_n)$")
        plt.savefig("./Exact/Figures/band_structure.png",dpi=800)
        plt.show()