'''
TOV SOLVER
'''
#------------------------------------------------------------------------------------------------------------
import numpy as np 

try:
    import matplotlib.pyplot as plt 
    pyplot_import_error = False

except ImportError:
    pyplot_import_error = True
    print("Unable to import matplotlib.pyplot @ TOV_Solver. Make sure that the \
     package is installed to the current env.")

#------------------------------------------------------------------------------------------------------------

#from EoS.EoS_Asym import EoS_Asym as EoS
#from EoS_parametric.EoS_asym_parametric import EoSAsymParametric as EoS
from EoS_Polytrope.EoS_polytrope import EoSPolytrope as EoS
from StarParams.params import Parameters as Prm 
from Root_finding.Rootfinding import RootFinder as Root
from Plot.plot import Plot as Plot 

#------------------------------------------------------------------------------------------------------------

class Star():

    barM = []
    barP0 = 0.0
    barP = []
    r = []
    r0 = 0.0
    dr = 0.0

#------------------------------------------------------------------------------------------------------------


    def __init__(self, barP0):

        self.barP0 = barP0
        cent_values = self.get_star_central_values() #Central values of the star.
        self.barM.append(cent_values["barM_init"])
        self.barP.append(barP0)
        self.r.append(cent_values["r0"])
        self.dr = cent_values["dr"]

        self.print_central_values()

#------------------------------------------------------------------------------------------------------------


    def __del__(self):
            objname = self.__class__.__name__
            print(" ", (objname +" Removed.").center(30,"-"))

#------------------------------------------------------------------------------------------------------------
    
    def get_star_central_values(self):

        prmObj = Prm()
        init_values = prmObj.get_star_init_values()
        star_params = prmObj.get_constants()
        eosObj = EoS()


        Ms = star_params["Ms"]
        c  = star_params["c"]
        pi = star_params["pi"]

        r0 = init_values ["r0"]
        dr = init_values ["dr"]

        EoS_params = eosObj.get_EoS_params()

        #Value of barM0 for given barP0 according to the EoS
        barM0 = (4.*pi*EoS_params["e0"]/(3.*Ms*c**2))*r0**3.*eosObj.barE(self.barP0) 
        
        star_central_vals = {
            "r0":r0,
            "dr":dr,
            "barM_init":barM0,
        }

        return star_central_vals

#------------------------------------------------------------------------------------------------------------

    def print_central_values(self):

        left_width = 15
        rigt_width = 15
        star_central_values = self.get_star_central_values()

        print(" ","STAR CENTRAL VALUES".center(left_width+rigt_width, '-'))
        print(" ","="*(left_width+rigt_width))

        for item, value in  star_central_values.items():
            print(" ",item.ljust(left_width, '.')+" | "+'{:8.6e}'.format(value).ljust(rigt_width))
        print("\n")

#------------------------------------------------------------------------------------------------------------

    def dbarM_dr(self, r, barP, barM):
        '''
        RHS of the mass balance equation in the TOV eqns.
        '''
        
        eosObj = EoS()
        prmObj = Prm()
        params = prmObj.get_constants()

        Ms = params["Ms"]
        c  = params["c"]
        pi = params["pi"]

        eosparams = eosObj.get_EoS_params()
        e0 = eosparams["e0"]
        
        barE = eosObj.barE(barP)

        dbarMdr = (4.0*pi*e0/(Ms*c**2))*(r**2.0*barE)

        return dbarMdr

#------------------------------------------------------------------------------------------------------------

    def dbarP_dr(self, r, barM, barP):
        '''
        RHS of the force balance eqn in the TOV eqns.
        '''
        eosObj  = EoS()
        prmObj  = Prm()

        params  = prmObj.get_constants()

        R0 = params["R0"]
        pi = params["pi"]
        c  = params["c"]
        Ms = params["Ms"]

        eosObj = EoS()
        eosparams = eosObj.get_EoS_params()
        e0 = eosparams["e0"]
        barE   = eosObj.barE(barP)

        f1 = -R0*barM*barE/r**2
        f2 = 1.0 +(barP/(barE))
        f3 = 1.0 + (4.0*pi*e0/(Ms*c**2))*(r**3)*(barP/barM)
        f4 = 1.0 - 2.0*R0*(barM/r)

        dbarP_dr = f1*f2*f3/f4

        return dbarP_dr


    #------------------------------------------------------------------------------------------------------------

    def TOV_solver(self):

        '''
        Estimate the next value of bar P, barM and r by solving TOV eqns
        using simple Euler method. The iteration of the while-loop terminates when barP
        hit negative value.
        '''

        outfname = "./Output/barP0.dat"

        break_loop = False
        i = 0
        dr = self.dr

        while(True):

            i += 1

            r_lst = self.r[-1]
            barM_lst = self.barM[-1]
            barP_lst = self.barP[-1]

            '''
            Runge-Kutta implementation.

            The differential eqns. are coupled. Runge-Kutta is implemented
            here only for one variable at a time per iteration. That is in each iteration, value of next
            barP is calculated assuming barM is a constant during that iteration. Similarely for other
            parameters.
            '''

            
            barPK1 = dr * self.dbarP_dr(r_lst, barM_lst, barP_lst)

            if barP_lst+(0.5*barPK1) >= 0.0:
                break_loop = False
                barPK2 = dr * self.dbarP_dr(r_lst+(0.5*dr), barM_lst, barP_lst+(0.5*barPK1))
            else:
                break_loop = True
                break

            if barP_lst+(0.5*barPK2) >= 0.0:
                break_loop = False
                barPK3 = dr * self.dbarP_dr(r_lst+(0.5*dr), barM_lst, barP_lst+(0.5*barPK2))
            else:
                break_loop = True
                break
             
            if barP_lst+barPK3 >= 0.0:
                break_loop = False
                barPK4 = dr * self.dbarP_dr(r_lst+dr, barM_lst, barP_lst+barPK3)
            else:
                break_loop = True
                break
                
            barMK1 = dr * self.dbarM_dr(r_lst, barP_lst, barM_lst)
            barMK2 = dr * self.dbarM_dr(r_lst+(0.5*dr), barP_lst, barM_lst+(0.5*barMK1))
            barMK3 = dr * self.dbarM_dr(r_lst+(0.5*dr), barP_lst, barM_lst +(0.5*barMK2))
            barMK4 = dr * self.dbarM_dr(r_lst+dr, barP_lst, barM_lst+barMK3)

            barP_nxt = barP_lst + (1.0/6.0)*(barPK1+(2.0*barPK2)+(2.0*barPK3)+barPK4)
            barM_nxt = barM_lst + (1.0/6.0)*(barMK1+(2.0*barMK2)+(2.0*barMK3)+barMK4)

            r_nxt = r_lst + dr

            if (barP_nxt >= 0):
                break_loop = False
                self.barP.append(barP_nxt)
                self.barM.append(barM_nxt)
                self.r.append(r_nxt)
                
            else:
                break_loop = True
                break
        
        if break_loop == True:
                print(" > r={:5.4e} barM={:5.4e} barP={:5.4e}".format(self.r[-1],self.barM[-1], self.barP[-1]))
                outdata = [self.r[-1], self.barM[-1], outfname]
                output = np.array([self.r,self.barM, self.barP])
                outdata = [self.r[-1], self.barM[-1], outfname]
                np.savetxt(outfname, output.T, fmt="%6.5e", delimiter="\t")
                return outdata

#------------------------------------------------------------------------------------------------------------


def main():

    boundary = "="*50

    print("\n", " > ", boundary, "\n")

    eosObj = EoS()
    prmObj = Prm()

    eosObj.print_EoS_params()
    prmObj.print_constants()
    prmObj.print_conv_factors()

    print("\n > ",boundary)
    
    #------------------------------------------------------------------------------------------------------------

    us = []
    u0 = 20.0
    scale = 5.0
    i = 0

    for i in range(30):
        us.append(u0*np.exp(-i/scale))

    us = us[::-1]
    #barP0s = np.array([eosObj.barP(u) for u in us ])
    barP0s = [
        5.0e-5, 8.0e-5, 1.0e-4, 5.0e-4, 0.001,0.002, 0.004, 0.008, 0.01, 0.02, 0.04,0.08,0.1,0.50,
        1.0, 8.0, 16.0, 32.0, 64.0, 100.0, 200.0, 400.0, 800.0, 1.0e3, 2.0e3,
        ]

    R = []
    barM = []
    Nil = []
    
    for barP0 in barP0s:
        NS = Star(barP0) # Instance of class NS
        print("\n > barP0: {:8.7e}".format(barP0))
        print(" --------------------------------------------------\n")
        outdata = NS.TOV_solver() # Build star

        R.append(outdata[0])
        barM.append(outdata[1])
        Nil.append(float("nan"))
        print(" --------------------------------------------------\n")
        del NS
        print("\n\n")

    RVsbarM = np.array([R, barM, Nil])
    RVsbarM_fname = "./Output/RVsbarM.dat"
    np.savetxt(RVsbarM_fname, RVsbarM.T, fmt="%6.5e", delimiter="\t")
    plotObj = Plot()  # Used for plotting
    plotObj.single_plot(RVsbarM_fname)

    print("\n > ",boundary)

#------------------------------------------------------------------------------------------------------------

if __name__ == "__main__":
    main()

#------------------------------------------------------------------------------------------------------------

