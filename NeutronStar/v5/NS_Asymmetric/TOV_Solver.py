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

from EoS.EoS_Asym import EoS_Asym as EoS_Asym
from StarParams.params import Parameters as Prm 
from Root_finding.Rootfinding import RootFinder as Root
from Plot.plot import Plot as Plot 

#------------------------------------------------------------------------------------------------------------


class Star():

    u = []
    barM = []
    barP = []
    r = []
    r0 = 0.0
    dr = 0.0

#------------------------------------------------------------------------------------------------------------


    def __init__(self, u_init):

        self.u_lower = 0.0
        self.u_upper = 100.0
        self.u_init = u_init

        cent_values = self.get_star_central_values() #Central values of the star.

        self.u.append(cent_values["u_init"])
        self.barM.append(cent_values["barM_init"])
        self.barP.append(cent_values["barP_init"])
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
        eosObj = EoS_Asym()


        Ms = star_params["Ms"]
        c  = star_params["c"]
        pi = star_params["pi"]

        r0 = init_values ["r0"]
        dr = init_values ["dr"]

        barP0 = eosObj.barP(self.u_init)
        EoS_params = eosObj.get_EoS_params()

        #Value of barM0 for given barP0 according to the EoS
        barM0 = (4.*pi*EoS_params["e0"]/(3.*Ms*c**2))*r0**3.*eosObj.barE(self.u_init) 
        
        star_central_vals = {
            "r0":r0,
            "dr":dr,
            "barP_init":barP0,
            "barM_init":barM0,
            "u_init":self.u_init,
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
            #print(" "+"_"*(left_width+rigt_width+4))
        print("\n")

#------------------------------------------------------------------------------------------------------------

    def dbarM_dr(self, r, u, barM):
        '''
        RHS of the mass balance equation in the TOV eqns.
        '''
        
        eosObj = EoS_Asym()
        prmObj = Prm()
        params = prmObj.get_constants()

        Ms = params["Ms"]
        c  = params["c"]
        pi = params["pi"]

        eosparams = eosObj.get_EoS_params()
        e0 = eosparams["e0"]
        
        barE = eosObj.barE(u)

        dbarMdr = (4.0*pi*e0/(Ms*c**2))*(r**2.*barE)

        return dbarMdr

#------------------------------------------------------------------------------------------------------------

    def dbarP_dr(self, r, barM, barP):
        '''
        RHS of the force balance eqn in the TOV eqns.
        '''
        eosObj  = EoS_Asym()
        prmObj  = Prm()
        rootObj = Root()

        params  = prmObj.get_constants()

        u = rootObj.Bisection(eosObj.barP, barP, self.u_lower, self.u_upper)
        
        R0 = params["R0"]
        pi = params["pi"]
        c  = params["c"]
        Ms = params["Ms"]

        eosObj = EoS_Asym()
        eosparams = eosObj.get_EoS_params()
        e0 = eosparams["e0"]
        barE   = eosObj.barE(u)

        f1 = -R0*barM*barE/r**2
        f2 = 1. +(barP/barE)
        f3 = 1. + (4.*pi*e0/(Ms*c**2))*(barM/r)
        f4 = 1. - 2.*R0*(barM/r)

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

        eosObj = EoS_Asym()
        rootObj = Root()

        i = 0
        dr = self.dr

        while(True):

            i += 1

            r_lst = self.r[-1]
            u_lst = self.u[-1]
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

            if barP_lst+(0.5*barPK1) >=0.0:
                barPK2 = dr * self.dbarP_dr(r_lst+(0.5*dr), barM_lst, barP_lst+(0.5*barPK1))
            else:
                output = np.array([self.r,self.barM, self.barP])
                np.savetxt(outfname, output.T, fmt="%6.5e", delimiter="\t")
                print(" > r={:5.4e} barM={:5.4e} barP={:5.4e}".format(self.r[-1],self.barM[-1], self.barP[-1]))
                outdata = [self.r[-1], self.barM[-1], outfname]
                return outdata

            if barP_lst+(0.5*barPK2) >=0.0:
                barPK3 = dr * self.dbarP_dr(r_lst+(0.5*dr), barM_lst, barP_lst+(0.5*barPK2))
            else:
                output = np.array([self.r,self.barM, self.barP])
                np.savetxt(outfname, output.T, fmt="%6.5e", delimiter="\t")
                print(" > r={:5.4e} barM={:5.4e} barP={:5.4e}".format(self.r[-1],self.barM[-1], self.barP[-1]))
                outdata = [self.r[-1], self.barM[-1], outfname]
                return outdata
            
            if barP_lst+barPK3 >= 0.0:
                barPK4 = dr * self.dbarP_dr(r_lst+dr, barM_lst, barP_lst+barPK3)
            else:
                output = np.array([self.r,self.barM, self.barP])
                np.savetxt(outfname, output.T, fmt="%6.5e", delimiter="\t")
                print(" > r={:5.4e} barM={:5.4e} barP={:5.4e}".format(self.r[-1],self.barM[-1], self.barP[-1]))
                outdata = [self.r[-1], self.barM[-1], outfname]
                return outdata 


            barMK1 = dr * self.dbarM_dr(r_lst, u_lst, barM_lst)
            barMK2 = dr * self.dbarM_dr(r_lst+(0.5*dr), u_lst, barM_lst+(0.5*barMK1))
            barMK3 = dr * self.dbarM_dr(r_lst+(0.5*dr), u_lst, barM_lst +(0.5*barMK2))
            barMK4 = dr * self.dbarM_dr(r_lst+dr, u_lst, barM_lst+barMK3)

            barP_nxt = barP_lst + (1.0/6.0)*(barPK1+(2.0*barPK2)+(2.0*barPK3)+barPK4)
            barM_nxt = barM_lst + (1.0/6.0)*(barMK1+(2.0*barMK2)+(2.0*barMK3)+barMK4)
            u_nxt = rootObj.Bisection(eosObj.barP, barP_nxt, self.u_lower, self.u_upper)
            #print("u_nxt: {} \t barP_nxt: {:5.4e}".format(u_nxt, barP_nxt))
            r_nxt = r_lst + dr

            #print("i: {} \t r: {} \t barP: {}".format(i,r_nxt, barP_nxt))
            
            if (barP_nxt >= 0) and (i <= 100000000):
                self.u.append(u_nxt)
                self.barP.append(barP_nxt)
                self.barM.append(barM_nxt)
                self.r.append(r_nxt)
                
            else:
                print(" > r={:5.4e} barM={:5.4e} barP={:5.4e}".format(self.r[-1],self.barM[-1], self.barP[-1]))
                outdata = [self.r[-1], self.barM[-1], outfname]
                output = np.array([self.r,self.barM, self.barP])
                print(" > r={:5.4e} barM={:5.4e} barP={:5.4e}".format(self.r[-1],self.barM[-1], self.barP[-1]))
                outdata = [self.r[-1], self.barM[-1], outfname]
                np.savetxt(outfname, output.T, fmt="%6.5e", delimiter="\t")
                return outdata

#------------------------------------------------------------------------------------------------------------


def main():

    boundary = "="*100

    eos_obj = EoS_Asym()

    print("\n", " > ", boundary, "\n")
  
    prmObj = Prm()
    prmObj.print_constants()
    prmObj.print_conv_factors()

    print("\n > ",boundary)
    
    #------------------------------------------------------------------------------------------------------------

    us = np.arange(0.05, 20.5, 0.05)

    R = []
    barM = []
    Nil = []
    
    for u in us:
        barP0 = barP0 = eos_obj.barP(u)
        NS = Star(barP0) # Instance of class NS
        #NS.print_central_values()
        print("\n > u: {}".format(u))
        #NS.print_params() # cross check the params.
        print(" --------------------------------------------------\n")
        outdata = NS.TOV_solver() # Build star

        R.append(outdata[0])
        barM.append(outdata[1])
        Nil.append(float("nan"))
        print(" --------------------------------------------------\n")
        # plot_obj = Plot() # Used for plotting
        # plot_obj.plot(outfname)
        del NS
        print("\n\n")

    RVsbarM = np.array([R, barM, Nil])
    RVsbarM_fname = "./Output/RVsbarM.dat"
    np.savetxt(RVsbarM_fname, RVsbarM.T, fmt="%6.5e", delimiter="\t")
    plotObj = Plot()  # Used for plotting
    plotObj.single_plot(RVsbarM_fname)
    #plotObj.double_plot(outdata[2])

    print("\n > ",boundary)

#------------------------------------------------------------------------------------------------------------

if __name__ == "__main__":
    main()

#------------------------------------------------------------------------------------------------------------

