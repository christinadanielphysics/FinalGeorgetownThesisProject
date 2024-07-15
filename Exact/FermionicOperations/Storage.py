from .c_dagger_ground import *
from .c_ground import *

from numpy import dot

class Storage:
    def __init__(self,problem):
        self.problem = problem
        self.ground_state_vector = problem.hubbard_N.ground_state_vector
        self.V = problem.parameters.V
        self.c_dagger_up_ground_eigenstates = []
        self.c_dagger_down_ground_eigenstates = []
        self.c_up_ground_eigenstates = []
        self.c_down_ground_eigenstates = []
    def form_c_dagger_up_ground_eigenstates(self):
        for spatial_index in range(self.V):
            vector = c_dagger_up_ground(spatial_index, self.problem)
            self.c_dagger_up_ground_eigenstates.append(vector)
    def form_c_dagger_down_ground_eigenstates(self):
        for spatial_index in range(self.V):
            vector = c_dagger_down_ground(spatial_index, self.problem)
            self.c_dagger_down_ground_eigenstates.append(vector)
    def form_c_up_ground_eigenstates(self):
        for spatial_index in range(self.V):
            vector = c_up_ground(spatial_index, self.problem)
            self.c_up_ground_eigenstates.append(vector)
    def form_c_down_ground_eigenstates(self):
        for spatial_index in range(self.V):
            vector = c_down_ground(spatial_index, self.problem)
            self.c_down_ground_eigenstates.append(vector)
    def plus_c_dagger_up_ground(self,eigenvector_plus,spatial_index):
        return dot(self.c_dagger_up_ground_eigenstates[spatial_index],eigenvector_plus)
    def plus_c_dagger_down_ground(self,eigenvector_plus,spatial_index):
        return dot(self.c_dagger_down_ground_eigenstates[spatial_index],eigenvector_plus)
    def minus_c_up_ground(self,eigenvector_minus,spatial_index):
        return dot(self.c_up_ground_eigenstates[spatial_index],eigenvector_minus)
    def minus_c_down_ground(self,eigenvector_minus,spatial_index):
        return dot(self.c_down_ground_eigenstates[spatial_index],eigenvector_minus)
    


