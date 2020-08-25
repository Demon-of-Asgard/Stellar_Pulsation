#!/home/demon/anaconda3/bin/python3
'''
This code constructs the mass and pressure distribution of a neutron star using TOV equations.
    Refs:
        1)  Neutronstar for undergraduate.
            Authors: Richard R. Silbar, Sanjay Reddy.
            Arxiv: 0309041
            Url: https://arxiv.org/abs/nucl-th/0309041

        2)  Compact stars for undergraduates.
            Authors: I. Segert, M. Hampel, C. Greiner, J. S. Beilich.
            Arxiv: 0506417v1
            Url: https://arxiv.org/abs/astro-ph/0506417
'''

#-----------------------------------------------------------------------------------

import numpy as np
from EoS.eos import Eos as Eos
from Plot.plot import Plot as Plot

#-----------------------------------------------------------------------------------

class Star(Eos):
    '''
    Calss representing a Neutron star. Inherits properties of 
    class Eos, which inherits properties from class Params.
    '''

    barP = []
    barM = []
    r = []
    dr = 0.0

    #-----------------------------------------------------------------------------------

    def __init__(self, barP0):

        super().__init__()

        self.barP = []
        self.barM = []
        self.r = []
        self.dr = 0.0

        init_vals = self.get_init_vals(barP0)

        self.barP.append(init_vals["barP0"])
        self.barM.append(init_vals["barM0"])
        self.r.append(init_vals["r0"])
        self.dr = init_vals["dr"]

    #-----------------------------------------------------------------------------------

    def __del__(self):
        objname = self.__class__.__name__
        print(" > ", objname, "Deleted.")

    #-----------------------------------------------------------------------------------

    def get_init_vals(self, barP0):

        init_params = self.init_params()

        params = self.get_params()

        Ms = params["Ms"]
        c = params["c"]
        pi = params["pi"]
        e0_rl = params["e0_rl"]

        #barP0 = init_params["barP0"]
        r0 = init_params["r0"]
        dr = init_params["dr"]

        barM0 = (4.*pi*e0_rl/(3.*Ms*c**2))*r0**3.*self.eos(barP0)

        init_vals = {
            "r0":r0,
            "dr":dr,
            "barP0":barP0,
            "barM0":barM0,
        }
        print(" > Initial values: {}".format(init_vals))
        return init_vals


    #-----------------------------------------------------------------------------------

    def print_params(self):
        '''
        To cross check values of the star-parameters.
        '''

        params = self.get_params()
        for index in params:
            print(" > [{}:={:5.4e}]\n".format(index, params[index]))

    #-----------------------------------------------------------------------------------

    def dbarM_dr(self, r, barM, barP, params):

        '''
        RHS of the mass balance equation in the TOV.
        '''

        Ms = params["Ms"]
        c = params["c"]
        pi = params["pi"]
        e0 = params["e0_rl"]

        barE = self.eos(barP)

        dbarMdr = (1.0e15*4.*pi*e0/(Ms*c**2))*(r**2.*barE)

        return dbarMdr

    #-----------------------------------------------------------------------------------

    def dbarP_dr(self, r, barM, barP, params):

        '''
        RHS of the force balance eqn in the TOV.
        '''

        R0 = params["R0"]
        pi = params["pi"]
        c  = params["c"]
        e0_rl = params["e0_rl"]
        Ms = params["Ms"]

        barM = self.barM[-1]
        barP = self.barP[-1]
        r = self.r[-1]

        barE = self.eos(barP)

        f1 = -R0*barM*barE/r**2
        f2 = 1. +(barP/barE)
        f3 = 1. + (4.*pi*e0_rl/(Ms*c**2))*(r**3.0)*(barP/barM)
        f4 = 1. - 2.*R0*(barM/r)

        dbarPdr = f1*f2*f3/f4

        return dbarPdr

    #-----------------------------------------------------------------------------------


    def build(self):

        '''
        Estimate the next value of bar P, barM and r by solving 
        TOV eqnsusing Runge-Kutta method. The iteration of the 
        while-loop terminates when barP hit negative value.
        '''

        params = self.get_params()
        i = 0
        dr = self.dr
        break_loop = False

        while True:

            i += 1

            r_lst = self.r[-1]
            barM_lst = self.barM[-1]
            barP_lst = self.barP[-1]

            '''
            Runge-Kutta implementation.

            The differential eqns. are coupled. Runge-Kutta is 
            implemented here only for one variable at a time 
            per iteration. That is in each iteration, value of next
            barP is calculated assuming barM is a constant 
            during that iteration. Similarely for otherparameters.
            '''
            barPK1 = dr * self.dbarP_dr(r_lst, barM_lst, 
            barP_lst, params)

            if barP_lst+(0.5*barPK1) >= 0.0:
                break_loop = False
                barPK2 = dr * self.dbarP_dr(r_lst+(0.5*dr), 
                barM_lst, barP_lst+(0.5*barPK1), params)
            else:
                break_loop = True
                break

            if barP_lst+(0.5*barPK2) >= 0.0:
                break_loop = False
                barPK3 = dr * self.dbarP_dr(r_lst+(0.5*dr), 
                barM_lst, barP_lst+(0.5*barPK2), params)
            else:
                break_loop = True
                break

            if barP_lst+barPK3 >= 0.0:
                break_loop = False
                barPK4 = dr * self.dbarP_dr(r_lst+dr, 
                barM_lst, barP_lst+barPK3, params)
            else:
                break_loop = True
                break

            # barPK1 = dr * self.dbarP_dr(r_lst, barM_lst, barP_lst, params)
            # barPK2 = dr * self.dbarP_dr(r_lst+(0.5*dr), barM_lst, barP_lst+(0.5*barPK1), params)
            # barPK3 = dr * self.dbarP_dr(r_lst+(0.5*dr), barM_lst, barP_lst+(0.5*barPK2), params)
            # barPK4 = dr * self.dbarP_dr(r_lst+dr, barM_lst, barP_lst+barPK3, params)

            barMK1 = dr * self.dbarM_dr(r_lst, barM_lst,
            barP_lst, params)

            barMK2 = dr * self.dbarM_dr(r_lst+(0.5*dr), 
            barM_lst+(0.5*barMK1), barP_lst, params)

            barMK3 = dr * self.dbarM_dr(r_lst+(0.5*dr), 
            barM_lst +(0.5*barMK2), barP_lst, params)

            barMK4 = dr * self.dbarM_dr(r_lst+dr, 
            barM_lst+barMK3, barP_lst, params)

            barP_nxt = barP_lst + (1.0/6.0)*(barPK1+(2.0*barPK2)
            +(2.0*barPK3)+barPK4)

            barM_nxt = barM_lst + (1.0/6.0)*(barMK1+(2.0*barMK2)
            +(2.0*barMK3)+barMK4)

            r_nxt = r_lst + dr

            if (barP_nxt >= 0.0):
                self.barP.append(barP_nxt)
                self.barM.append(barM_nxt)
                self.r.append(r_nxt)

            else:
                break_loop = True
                break

        #-----------------------------------------------------------------------------------
        if break_loop == True:
            print(" > r={:5.4e} barM={:5.4e} barP={:5.4e}".format(self.r[-1],self.barM[-1], self.barP[-1]))

            outfname = "./Output/.barP0_"+str(self.barP[0])+".dat"
            output = np.array([self.r,self.barM, self.barP])
            np.savetxt(outfname, output.T, fmt="%6.5e", delimiter="\t")

            outdata = [self.r[-1], self.barM[-1], outfname]

            return outdata

#-----------------------------------------------------------------------------------

def main():

    barP0s = [
        5.0e-5, 8.0e-5, 1.0e-4, 5.0e-4, 0.001,0.002, 0.004, 0.008, 0.01, 0.02, 0.04,0.08,0.1,0.50,
        1.0, 8.0, 16.0, 32.0, 64.0, 100.0, 200.0, 400.0, 800.0, 1.0e3,2.0e3,
        ]

    R = []
    barM = []
    Nil = []

    for barP0 in barP0s:
        boundary = "="*100
        NS = Star(barP0) # Instance of class NS
        #NS.print_params() # cross check the params.
        print(boundary,"\n")
        outdata = NS.build() # Build star

        R.append(outdata[0])
        barM.append(outdata[1])
        Nil.append(float("nan"))
        print(boundary,"\n")

        # plot_obj = Plot() # Used for plotting
        # plot_obj.plot(outfname)
        del NS

    RVsbarM = np.array([R, barM, Nil])
    RVsbarM_fname = "./Output/RVsbarM.dat"
    np.savetxt(RVsbarM_fname, RVsbarM.T, fmt="%6.5e", delimiter="\t")
    plot_obj = Plot()  # Used for plotting
    plot_obj.single_plot(RVsbarM_fname)

#-----------------------------------------------------------------------------------

if __name__ == "__main__":
    main()

#-----------------------------------------------------------------------------------
