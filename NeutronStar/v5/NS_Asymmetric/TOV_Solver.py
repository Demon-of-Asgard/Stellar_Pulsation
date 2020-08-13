'''
TOV SOLVER
'''
#-------------------------------------------------------------------------------------------------------------------------------------------
try:
    import matplotlib.pyplot as plt 
    pyplot_import_error = False
except ImportError:
    pyplot_import_error = True
    print("Unable to import matplotlib.pyplot @ TOV_Solver. Make sure that the \
     package is installed to the current env.")

import numpy as np 

#-------------------------------------------------------------------------------------------------------------------------------------------

try:
    from EoS.eos import Eos_Asym as Eos_Asym
except ImportError:
    print("Unable to import Eos from EoS.eos @ TOV_Solver.")
    exit(1)

try:
    from StarParams.params import Params as Prm
except:
    print("Unable to import Params from StarParams.params @ TOV_Solver.")
    exit(1)
#-------------------------------------------------------------------------------------------------------------------------------------------

class Star():

    u = []
    barM = []
    barP = []
    r = []
    r0 = 0.0
    dr = 0.0

#-------------------------------------------------------------------------------------------------------------------------------------------

    def __init__(self, u0):

        self.u0 = u0
        cent_values = self.get_star_central_values() #Central values of the star.

        self.u.append(cent_values["u_init"])
        self.barM.append(cent_values["barM_init"])
        self.barP.append(cent_values["barP_init"])
        self.r0 = cent_values["r0"]
        self.dr = cent_values["dr"]

#-------------------------------------------------------------------------------------------------------------------------------------------

    def __del__(self):
            objname = self.__class__.__name__
            print(" ", (objname +" Removed.").center(30,"-"))

#-------------------------------------------------------------------------------------------------------------------------------------------
    
    def get_star_central_values(self):

        _prm = Prm()
        init_values = _prm.get_star_init_values()
        star_params = _prm.get_star_params()

        Ms = star_params["Ms"]
        c = star_params["c"]
        pi = star_params["pi"]

        r0 = init_values ["r0"]
        dr = init_values ["dr"]

        barP0 = 0.0
        barM0 = 0.0 #(4.*pi*e0_rl/(3.*Ms*c**2))*r0**3.*self.eos(barP0) #Value of barM0 for given barP0 according to the EoS

        star_central_vals = {
            "r0":r0,
            "dr":dr,
            "barP_init":barP0,
            "barM_init":barM0,
            "u_init":self.u0,
        }

        return star_central_vals

#-------------------------------------------------------------------------------------------------------------------------------------------

    def print_central_values(self):

        left_width = 15
        rigt_width = 15
        star_central_values = self.get_star_central_values()

        print(" ","STAR CENTRAL VALUES".center(left_width+rigt_width, '-'))
        print(" ","="*(left_width+rigt_width))

        for item, value in  star_central_values.items():
            print(" ",item.ljust(left_width, '.')+" | "+'{:8.6e}'.format(value).ljust(rigt_width))
            #print(" "+"_"*(left_width+rigt_width+4))
        print("\n")
#-------------------------------------------------------------------------------------------------------------------------------------------

def main():
    boundary = "="*100
    print(" > ", boundary, "\n")
    NS = Star(0.0)
    
    star_params = Prm()
    star_params.print_star_params()
    star_params.print_conv_factors()
    NS.print_central_values()
    del NS

    print("\n > ",boundary)

#-------------------------------------------------------------------------------------------------------------------------------------------

if __name__ == "__main__":
    main()

#-------------------------------------------------------------------------------------------------------------------------------------------
