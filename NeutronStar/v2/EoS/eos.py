import numpy as np
from Params.params import params as prm

class eqn_state(prm):
    def __init__(self):
        super().__init__()
    def eos(self,p_bar = 0.0):
        params = self.get_params()
        A_r =params["A_r"]
        A_nr = params["A_nr"]
        print(A_r, A_nr)
    
    
