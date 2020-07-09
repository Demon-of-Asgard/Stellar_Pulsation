try:
    from StarParams.params import Params as prm
except ImportError:
    print("ImportErrot-Params")

#-----------------------------------------------------------------------------------

class Eos(prm):

    def __init__(self):
        super().__init__()

    def eos(self, barP):

        params = self.get_params()
        Krl = params["Krl"]
        e0_rl = params["e0_rl"] 
        gamma_rl = params["gamma_rl"]

        barK = Krl*e0_rl**(gamma_rl-1) 
        bare = (barP/barK)**(1./gamma_rl) 

        return bare

#-----------------------------------------------------------------------------------