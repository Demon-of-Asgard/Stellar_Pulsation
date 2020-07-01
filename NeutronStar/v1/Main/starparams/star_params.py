if __name__ == "__main__":
    try:
        import numpy as np
    except:
        print("ImportError-Numpy")

class params:
    try:
        import numpy as np
    except:
        print("ImportError-Numpy")

    pi = np.pi
    c = 3.0e10  # [cm/s]
    h = 6.62607004e-27  # [erg s]
    hbar = h/(2.0*pi)  # [erg s]

    mn = 1.6749e-24  # [gm]
    me = 9.10e-28  # [gm]
    A = 1.0
    Z = 1.0

    R0 = 1.473  # [cm]
    Ms = 1.989e33  # [gm]

    alpha_rl = 1.473
    alpha_nrl = 1.0  # [cm]
    gamma_rl = 4./3.
    gamma_nrl = 5./3.
    barP0rl = 1.0e-14 #Dimensionless central preasure.
    barP0nrl = 1.0e-6#""

    #---------------------------------------------------------------------


    def get_rlparams(self):

        c = self.c
        pi = self.pi
        hbar = self.hbar

        Ms = self.Ms
        R0 = self.R0
        A = self.A
        Z = self.Z
        mn = self.mn

        gamma_rl = self.gamma_rl
        alpha_rl = self.alpha_rl
        barP0_rl = self.barP0rl

        Krl = (hbar*c/(12.*pi**2))*(3.0*pi**2*Z/(A*mn*c**2))**gamma_rl
        e0_rl = ((R0/alpha_rl)**gamma_rl/Krl)**(1/(gamma_rl-1))
        barKrl = Krl*e0_rl**(gamma_rl-1)
        beta_rl = 1.0e15*(4.0*pi*e0_rl)/(Ms*c**2*barKrl**(1/gamma_rl))

        params = {
            "c": c,
            "R0": R0,
            "Ms": Ms,
            "alpha": alpha_rl,
            "beta": beta_rl,
            "gamma": gamma_rl,
            "eps0": e0_rl,
            "barP0": barP0_rl,
            "K": Krl,
            "barK": barKrl,
        }

        return params

#---------------------------------------------------------------------


    def get_nrlparams(self):

        c = self.c
        pi = self.pi
        hbar = self.hbar

        Ms = self.Ms
        R0 = self.R0
        A = self.A
        Z = self.Z
        mn = self.mn
        me = self.me
        
        gamma_nrl = self.gamma_nrl
        alpha_nrl = self.alpha_nrl
        barP0_nrl = self.barP0nrl

        Knrl = (hbar**2/(15.0*pi**2*mn))*(3.0*pi**2*Z/(A*mn*c**2))**gamma_nrl
        e0_nrl = ((R0/alpha_nrl)**gamma_nrl/Knrl)**(1/(gamma_nrl-1))
        barKnrl = Knrl*e0_nrl**(gamma_nrl-1)
        beta_nrl = 1.0e15*(4.*pi*e0_nrl)/(Ms*c**2*barKnrl**(1/gamma_nrl))
        
        params = {
            "c" : c,
            "R0": R0,
            "Ms": Ms,
            "alpha": alpha_nrl,
            "beta": beta_nrl,
            "gamma": gamma_nrl,
            "eps0": e0_nrl,
            "barP0": barP0_nrl,
            "K": Knrl,
            "barK": barKnrl,
        }

        return params

    #---------------------------------------------------------------------


if __name__ == "__main__":
    pm = params()
    
