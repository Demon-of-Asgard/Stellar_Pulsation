from StarParams.params import Parameters as Prm

#-----------------------------------------------------------------------------------

class EoSPolytrope():

    def __init__(self):
        self.__gamma_rl = 4./3.
        self.__gamma_nrl = 5./3.
        self.__A_rl = 2.891 #  --->neutrons only
        self.__A_nrl = 2.572 #2.4216

    def print_EoS_params(self):
        ''' Only to check the values. '''

        left_width = 15
        rigt_width = 15
        EoS_params = self.get_EoS_params()

        print(" ","EOS PARAMETERS".center(left_width+rigt_width, '-'))
        print(" ","="*(left_width+rigt_width))

        for item, value in  EoS_params.items():
            print(" ",item.ljust(left_width, '.')+" | "+'{:8.6e}'.format(value).ljust(rigt_width))
        print("\n")
    #-----------------------------------------------------------------------------------

    def get_EoS_params(self):

        paramObj = Prm()
        params = paramObj.get_constants()

        mN = params["mN"]
        c  = params["c"]
        pi = params["pi"]
        hbar  = params["hbar"]
        e0_rl = mN**4*c**5/(3.*pi**2*hbar**3)
        
        EoS_params = {
            "gamma_rl": self.__gamma_rl,
            "gamma_nrl": self.__gamma_nrl,
            "A_rl": self.__A_rl,
            "A_nrl": self.__A_nrl,
            "e0":e0_rl,
        }
        return EoS_params

    #-----------------------------------------------------------------------------------
    
    def barE(self, barP):
        barE = self.__A_rl*barP + self.__A_nrl*(barP**(1./self.__gamma_nrl))
        return barE

#-----------------------------------------------------------------------------------