try:
    from StarParams.params import Params as Prm
except ImportError:
    print("ImportErrot-Params")

#-----------------------------------------------------------------------------------

class Eos(Prm):

    def __init__(self):
        super().__init__()

    #-----------------------------------------------------------------------------------

    def get_EoS_params(self):
        eos_params = self.get_params()
        params = {
            "e0":eos_params["e0_rl"]
        }
        return params
    #-----------------------------------------------------------------------------------

    def eos(self, barP):

        params = self.get_params()
        A_rl = params["A_rl"]
        A_nrl = params["A_nrl"]
        gamma_nrl = params["gamma_nrl"]

        barE = A_rl*barP + A_nrl*(barP**(1./gamma_nrl))

        # A0 = 0.862
        # barE = A0*barP**(-1./2.)

        return barE

    #-----------------------------------------------------------------------------------
    
    def barE(self, barP):

        params = self.get_params()
        A_rl = params["A_rl"]
        A_nrl = params["A_nrl"]
        gamma_nrl = params["gamma_nrl"]

        barE = A_rl*barP + A_nrl*(barP**(1./gamma_nrl))

        # A0 = 0.862
        # barE = A0*barP**(-1./2.)

        return barE

    #-----------------------------------------------------------------------------------
