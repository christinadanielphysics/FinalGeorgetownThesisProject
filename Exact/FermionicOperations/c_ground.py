from .OccupationState import *

from numpy import zeros 

def c_up_ground(operator_index, problem):
    ground_state_vector = problem.hubbard_N.ground_state_vector
    vector_minus = zeros(problem.hubbard_minus.dimension)
    for basis_index,basis_state in enumerate(problem.hubbard_N.basis):
        initial_state = OccupationState(ground_state_vector[basis_index],basis_state.up_spin_list,basis_state.down_spin_list)
        end_state = initial_state.apply_up_annihilation(operator_index)
        if end_state.coefficient != 0:
            end_index = problem.hubbard_minus.basis_object.get_index(end_state)
            vector_minus[end_index] = vector_minus[end_index] + end_state.coefficient


    return vector_minus


def c_down_ground(operator_index, problem):
    ground_state_vector = problem.hubbard_N.ground_state_vector
    vector_minus = zeros(problem.hubbard_minus.dimension)
    for basis_index,basis_state in enumerate(problem.hubbard_N.basis):
        initial_state = OccupationState(ground_state_vector[basis_index],basis_state.up_spin_list,basis_state.down_spin_list)
        end_state = initial_state.apply_down_annihilation(operator_index)
        if end_state.coefficient != 0:
            end_index = problem.hubbard_minus.basis_object.get_index(end_state)
            vector_minus[end_index] = vector_minus[end_index] + end_state.coefficient
    
    
    return vector_minus