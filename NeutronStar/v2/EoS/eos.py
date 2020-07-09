import numpy as np
from Params.params import params as prm

class eqn_state(prm):

    def __init__(self):
        super().__init__()

    def eos(self,p_bar):
        params = self.get_params()
        A_r =params["A_r"]
        A_nr = params["A_nr"]
        e_bar = A_r*p_bar**self.gamma_r + A_nr*p_bar**(1.0/self.gamma_nr)
        #A_r*p_bar**self.gamma_r +
        return e_bar
    
    
