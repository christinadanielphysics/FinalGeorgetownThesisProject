from .OccupationState import *

from numpy import zeros

def c_dagger_up_ground(operator_index, problem):
    ground_state_vector = problem.hubbard_N.ground_state_vector
    vector_plus = zeros(problem.hubbard_plus.dimension)
    for basis_index,basis_state in enumerate(problem.hubbard_N.basis):
        initial_state = OccupationState(ground_state_vector[basis_index],basis_state.up_spin_list,basis_state.down_spin_list)
        end_state = initial_state.apply_up_creation(operator_index)
        if end_state.coefficient != 0:
            end_index = problem.hubbard_plus.basis_object.get_index(end_state)
            vector_plus[end_index] = vector_plus[end_index] + end_state.coefficient
    

    return vector_plus

def c_dagger_down_ground(operator_index, problem):
    ground_state_vector = problem.hubbard_N.ground_state_vector
    vector_plus = zeros(problem.hubbard_plus.dimension)
    for basis_index,basis_state in enumerate(problem.hubbard_N.basis):
        initial_state = OccupationState(ground_state_vector[basis_index],basis_state.up_spin_list,basis_state.down_spin_list)
        end_state = initial_state.apply_down_creation(operator_index)
        if end_state.coefficient != 0:
            end_index = problem.hubbard_plus.basis_object.get_index(end_state)
            vector_plus[end_index] = vector_plus[end_index] + end_state.coefficient


    return vector_plus
