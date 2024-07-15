from .OccupationState import *

from itertools import combinations
from numpy import zeros,array_equal

class Basis:
    def __init__(self,V,N_up,N_down):
        """
        Completed 12/22/22
        """
        self.basis = []
        self.spatial_indices = list(range(V))
        self.up_spin_tuples = combinations(self.spatial_indices,N_up)
        self.down_spin_tuples = combinations(self.spatial_indices,N_down)
        self.V = V
    def convert_tuples_to_bit_lists(self,combinations):
        """
        Completed 12/22/22
        """
        bit_lists = []
        for my_tuple in list(combinations):
            bit_list = zeros(self.V)
            for spatial_index in my_tuple:
                bit_list[spatial_index] = 1
            bit_lists.append(bit_list)
        return bit_lists
    def form(self):
        """
        Completed 12/22/22
        """
        up_lists = self.convert_tuples_to_bit_lists(self.up_spin_tuples)
        down_lists = self.convert_tuples_to_bit_lists(self.down_spin_tuples)
        for up_list in up_lists:
            for down_list in down_lists:
                self.basis.append( OccupationState(1, up_list, down_list) )
    def display(self):
        """
        Completed 12/22/22
        """
        for index,state in enumerate(self.basis):
            print("index "+str(index),end=" | state: ")
            state.display()
            print("\n")
    def get_index(self,trial_state):
        """
        Completed 12/22/22
        """
        trial_up_list = trial_state.up_spin_list
        trial_down_list = trial_state.down_spin_list
        for basis_index,basis_state in enumerate(self.basis):
            if array_equal(basis_state.up_spin_list,trial_up_list):
                if array_equal(basis_state.down_spin_list,trial_down_list):
                    return basis_index
        print("\nBasis index not found for: ",end="")
        trial_state.display()
        







