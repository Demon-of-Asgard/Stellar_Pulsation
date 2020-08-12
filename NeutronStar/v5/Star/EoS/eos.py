try:
    from StarParams.params import Params as Prm
except ImportError:
    print("ImportErrot-Params")

#-----------------------------------------------------------------------------------

class Eos(Prm):

    def __init__(self):
        super().__init__()

    def eos(self, barP):

        params = self.get_params()
        A_rl = params["A_rl"]
        A_nrl = params["A_nrl"]
        gamma_nrl = params["gamma_nrl"]

        barE = A_rl*barP + A_nrl*(barP**(1./gamma_nrl))

        return barE

#-----------------------------------------------------------------------------------