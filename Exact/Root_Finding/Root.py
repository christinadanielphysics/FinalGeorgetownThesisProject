class Root:
    def __init__(self,root_value,multiplicity,a,b):
        self.root_value = root_value
        self.multiplicity = multiplicity
        self.a = a
        self.b = b
    def convert_root_from_t_to_w(self):
        return self.a+((self.root_value+1)*(self.b-self.a)/2)

