
try:
    import numpy as np
except ImportError:
    print("ImportError-Numpy")

class Params:
    def __init__(self):

        self.__pi = np.pi
        self.__c = 3.0e10  # [cm/s]
        self.__Ms = 1.989e33  # [gm]
        self.__h = 6.62607004e-27  # [erg s]
        self.__hbar = self.__h/(2.0*self.__pi)  # [erg s]
        self.__mn = 1.6749e-24  # [gm]
        self.__me = 9.10e-28  # [gm]

        self.__A = 2.15
        self.__Z = 1.0
        self.__R0 = 1.476  # [km]
        self.__alpha_rl = self.__R0
        self.__gamma_rl = 4./3.
        self.__gamma_nrl = 5./3.
        self.__A_rl = 2.8663
        self.__A_nrl = 2.4216

        self.__barM0 = 0.0
        self.__barP0 = 0.01
        self.__dr = 0.01
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


        #Krl = (hbar*c/(12.*pi**2))*(3.0*pi**2*Z/(A*mn*c**2))**gamma_rl
        #e0_rl = ((R0/alpha_rl)**gamma_rl/Krl)**(1/(gamma_rl-1))
        e0_rl = mn**4*c**5/(3.*pi**2*hbar**3)
        #print(e0_rl)
        # barKrl = Krl*e0_rl**(gamma_rl-1)
        # beta_rl = 1.0e15*(4.0*pi*e0_rl)/(Ms*c**2*barKrl**(1/gamma_rl))

        dict_params = {
            "c":self.__c,
            "pi":self.__pi,
            "R0":self.__R0,
            "Ms":self.__Ms,
            #"Krl": Krl,
            "e0_rl" : e0_rl,
            #"barKrl": barKrl,
            "alpha_rl": alpha_rl,
            #"beta_rl" : beta_rl,
            "gamma_rl" : gamma_rl,
            "gamma_nrl":gamma_nrl,
            "A_rl": A_rl,
            "A_nrl": A_nrl,
        }

        return dict_params
    
    def init_params(self):
        r0 = self.__r0
        barP0 = self.__barP0
        barM0 = self.__barM0

        init_params = {

            "r0":r0,
            "barP0":barP0,
            "barM0":barM0,
        }

        return init_params

    #-------------------------------------------------------------

if __name__ == "__main__":
    obj = Params()
    print(obj.get_params())
#-------------------------------------------------------------
