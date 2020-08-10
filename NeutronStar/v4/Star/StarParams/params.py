
try:
    import numpy as np
except ImportError:
    print("ImportError-Numpy")

class Params:
    def __init__(self):
        
        # Constant values. 
        self.__pi = np.pi
        self.__c = 3.0e10  # [cm/s]
        self.__Ms = 1.989e33  # [gm]
        self.__h = 6.62607004e-27  # [erg s]
        self.__hbar = self.__h/(2.0*self.__pi)  # [erg s]
        self.__mn = 1.6749e-24  # [gm]
        self.__me = 9.10e-28  # [gm]
        self.__R0 = 1.476  # [km]

        #Model dependent choices.
        self.__A = 2.15
        self.__Z = 1.0
        self.__alpha_rl = self.__R0
        self.__gamma_rl = 4./3.
        self.__gamma_nrl = 5./3.

        # EoS PARAMS WITH NEUTRONS ONLY.
        # self.__A_rl = 2.4216  
        # self.__A_nrl = 2.8663

        # EoS PARAMS WITH NEUTRONS, PROTONS AND ELECTRONS.
        self.__A_rl = 2.891 #  --->neutrons only
        self.__A_nrl = 2.572 #2.4216

        #Initial values.
        self.__barM0 = 0.0
        self.__barP0 = 0.01
        self.__dr = 1.0e-4
        self.__r0 = self.__dr

    #-------------------------------------------------------------

    def get_params(self):

        hbar = self.__hbar
        c = self.__c
        pi = self.__pi
        mn = self.__mn
        Ms = self.__Ms
        R0 = self.__R0
        alpha_rl = self.__alpha_rl
        gamma_rl = self.__gamma_rl
        gamma_nrl = self.__gamma_nrl
        Z = self.__Z
        A = self.__A
        A_rl = self.__A_rl
        A_nrl = self.__A_nrl

        e0_rl = mn**4*c**5/(3.*pi**2*hbar**3)
     
        dict_params = {
            "c":self.__c,
            "pi":self.__pi,
            "R0":self.__R0,
            "Ms":self.__Ms,
            "e0_rl" : e0_rl,
            "alpha_rl": alpha_rl,
            "gamma_rl" : gamma_rl,
            "gamma_nrl":gamma_nrl,
            "A_rl": A_rl,
            "A_nrl": A_nrl,
        }

        return dict_params
    
    def init_params(self):
        r0 = self.__r0
        dr = self.__dr
        barP0 = self.__barP0
        barM0 = self.__barM0

        init_params = {

            "r0":r0,
            "dr":dr,
            "barP0":barP0,
            "barM0":barM0,
        }
        return init_params

    #-------------------------------------------------------------

if __name__ == "__main__":
    obj = Params()
    print(obj.get_params())
#-------------------------------------------------------------
