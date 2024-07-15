from .OccupationState import *
from .Basis import *
from ..Parameters import *
from .Hubbard import *
from .Storage import *
from Problem import *

import matplotlib.pyplot as plt
from numpy import conjugate
from cmath import exp, pi

print("\n",end="")
print("TEST 1: Print an OccupationState")
coefficient = -2
up_spin_list = [0,1,0]
down_spin_list = [1,0,1]
occupation_state = OccupationState(coefficient, up_spin_list, down_spin_list)
occupation_state.display()
print("\n")

print("TEST 2: copy_spin_list")
new_up_spin_list = occupation_state.copy_spin_list(up_spin_list)
new_up_spin_list.append(1)
print("new:",new_up_spin_list)
print("old:",up_spin_list)
print("\n",end="")

print("TEST 3: number_of_passes_for_up")
passes = occupation_state.number_of_passes_for_up(1)
print(passes)
print("\n",end="")

print("TEST 4: number_of_passes_for_down")
passes = occupation_state.number_of_passes_for_down(1)
print(passes)
print("\n",end="")

print("TEST 5: apply_up_creation")
print("(c†2↑)",end="")
occupation_state.display()
print(" = ",end="")
occupation_state.apply_up_creation(2).display()
print("\n")

print("TEST 6: apply_down_creation")
print("(c†1↓)",end="")
occupation_state.display()
print(" = ",end="")
occupation_state.apply_down_creation(1).display()
print("\n")

print("TEST 7: apply_up_annihilation")
print("(c1↑)",end="")
occupation_state.display()
print(" = ",end="")
occupation_state.apply_up_annihilation(1).display()
print("\n")

print("TEST 8: apply_down_annihilation")
print("(c0↓)",end="")
occupation_state.display()
print(" = ",end="")
occupation_state.apply_down_annihilation(0).display()
print("\n")

print("TEST 9: form and display")
V = 3
N_up = 1
N_down = 2
basis_object = Basis(V,N_up,N_down)
basis_object.form()
basis_object.display()

print("TEST 10: get_index")
print("trial state: ", end="")
occupation_state.display()
print("\n",end="")
index = basis_object.get_index(occupation_state)
print("index:",index)
print("state using index: ", end="")
basis_object.basis[index].display()
print("\n")

# test cases 1-10 were completed on 12/22/22








