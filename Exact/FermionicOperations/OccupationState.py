from termcolor import colored

class OccupationState:
    def __init__(self, coefficient, up_spin_list, down_spin_list):
        """
        Completed 12/22/22
        """
        self.coefficient = coefficient
        self.up_spin_list = up_spin_list
        self.down_spin_list = down_spin_list
    def copy_spin_list(self,spin_list):
        """
        Completed 12/22/22
        """
        return spin_list.copy()
    def number_of_passes_for_up(self,operator_index):
        """
        Completed 12/22/22
        """
        passes = 0
        for index,value in enumerate(self.up_spin_list):
            if operator_index == index:
                break
            else:
                if value == 1:
                    passes = passes + 1
        return passes
    def number_of_passes_for_down(self,operator_index):
        """
        Completed 12/22/22
        """
        passes = sum(self.up_spin_list)
        for index,value in enumerate(self.down_spin_list):
            if operator_index == index:
                break
            else:
                if value == 1:
                    passes = passes + 1
        return passes
    def apply_up_creation(self,operator_index):
        """
        Completed 12/22/22
        """
        if self.up_spin_list[operator_index] == 1:
            return OccupationState(0,[],[])
        else:
            new_up_spin_list = self.copy_spin_list(self.up_spin_list)
            new_up_spin_list[operator_index] = 1
            new_coefficient = self.coefficient * (-1)**self.number_of_passes_for_up(operator_index)
            return OccupationState(new_coefficient, new_up_spin_list, self.down_spin_list)
    def apply_down_creation(self,operator_index):
        """
        Completed 12/22/22
        """
        if self.down_spin_list[operator_index] == 1:
            return OccupationState(0,[],[])
        else:
            new_down_spin_list = self.copy_spin_list(self.down_spin_list)
            new_down_spin_list[operator_index] = 1
            new_coefficient = self.coefficient * (-1)**self.number_of_passes_for_down(operator_index)
            return OccupationState(new_coefficient, self.up_spin_list, new_down_spin_list)
    def apply_up_annihilation(self,operator_index):
        """
        Completed 12/22/22
        """
        if self.up_spin_list[operator_index] == 1:
            new_up_spin_list = self.copy_spin_list(self.up_spin_list)
            new_up_spin_list[operator_index] = 0
            new_coefficient = self.coefficient * (-1)**self.number_of_passes_for_up(operator_index)
            return OccupationState(new_coefficient,new_up_spin_list,self.down_spin_list)
        else:
            return OccupationState(0,[],[])
    def apply_down_annihilation(self,operator_index):
        """
        Completed 12/22/22
        """
        if self.down_spin_list[operator_index] == 1:
            new_down_spin_list = self.copy_spin_list(self.down_spin_list)
            new_down_spin_list[operator_index] = 0
            new_coefficient = self.coefficient * (-1)**self.number_of_passes_for_down(operator_index)
            return OccupationState(new_coefficient,self.up_spin_list,new_down_spin_list)
        else:
            return OccupationState(0,[],[])
    def display(self):
        """
        Completed 12/22/22
        """
        print(colored('('+str(self.coefficient)+')','green',attrs=['bold']), end="")
        for index,value in enumerate(self.up_spin_list):
            if value == 1:
                print(colored("(c†"+str(index)+"↑)",'blue',attrs=['bold']), end="")
        for index,value in enumerate(self.down_spin_list):
            if value == 1:
                print(colored("(c†"+str(index)+"↓)",'red',attrs=['bold']), end="")
        print(colored("|0>",attrs=['bold']), end="")
    
