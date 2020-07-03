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
    A = 2.1452
    Z = 1

    R0 = 1.476  # [cm]
    Ms = 1.989e33  # [gm]
    Anr = 2.4216
    Ar  = 2.8663
    e0 = 5.346e36 #erg/cm^3 

    alpha_ = R0
    gamma_r = 1.0
    gamma_nr = 5./3.

    def __init__(self):
        pass

    def get_params(self):
        
        params = {
            "pi" : self.pi,
            "c" : self.c,
            "hbar" :  self.hbar,
            "mn" : self.mn,
            "me" :  self.me,
            "R0" : self.R0,
            "Ms" : self.Ms,

            "A" : self.A,
            "Z" : self.Z,

            "A_nr" : self.Anr,
            "A_r" : self.Ar,
            "e0" : self.mn**4*self.c**5/(3.*self.pi**2*self.hbar**3),

            "alpha" : self.R0,
            "gamma_r" : self.gamma_r,
            "gamma_nr" : self.gamma_nr
        }

        return params
