#!/home/demon/anaconda3/bin/python3

#-----------------------------------------------------------------------------------

import numpy as np
from EoS.eos import Eos as Eos
from Plot.plot import Plot as Plot

#-----------------------------------------------------------------------------------

class NS(Eos):

    barP = []
    barM = []
    r = []

    #-----------------------------------------------------------------------------------

    def __init__(self):
        super().__init__()
        init_params = self.get_init_params()
        self.barP.append(init_params["barP0"])
        self.barM.append(init_params["barM0"])
        self.r.append(init_params["r0"])
        
    #-----------------------------------------------------------------------------------

    def print_params(self):
        params = self.get_params()
        for index in params:
            print(" > [{}:={:5.4e}]\n".format(index, params[index]))
    #-----------------------------------------------------------------------------------

    def dbarM_dr(self):
        params = self.get_params()
        Ms = params["Ms"]
        c = params["c"]
        pi = params["pi"]
        eps0 = params["e0_rl"]

        barP = self.barP[-1]
        r = self.r[-1]

        dbarMdr = (1.e15*4.*pi*eps0/(Ms*c**2))*(r**2.*self.eos(barP))
        return dbarMdr
        
    #-----------------------------------------------------------------------------------

    def dbarP_dr(self):

        params = self.get_params()
        R0 = params["R0"]
        barM = self.barM[-1]
        barP = self.barP[-1]
        r = self.r[-1]

        dbarPdr = -R0*barM*self.eos(barP)/r**2
        return dbarPdr

    #-----------------------------------------------------------------------------------

    def build(self):
        i = 0
        dr = 0.1
        while(True and i<= 10000000):
            i += 1
            barP_lst = self.barP[-1]
            barP_nxt = barP_lst + self.dbarP_dr()*dr
            barM_lst = self.barM[-1]
            barM_nxt = barM_lst + self.dbarM_dr()* dr 
            r_lst = self.r[-1]
            r_nxt = r_lst + dr

            if barP_nxt >= 0.0:
                self.barP.append(barP_nxt)
                self.barM.append(barM_nxt)
                self.r.append(r_nxt)

            else:
                break

        #-----------------------------------------------------------------------------------

        print(" > r={} barM={:5.4e} barP={:5.4e}".format(i,self.barM[-1], self.barP[-1]))

        outfname = "output.dat"
        output = np.array([self.r,self.barM, self.barP])
        np.savetxt(outfname, output.T, fmt="%6.5e", delimiter="\t")

        return outfname

#-----------------------------------------------------------------------------------

def main():
    star_obj = NS()
    star_obj.print_params()
    print("--------------------------------------------------\n")
    outfanme = star_obj.build()
    print("--------------------------------------------------\n")
    plot_obj = Plot()
    plot_obj.plot(outfanme)


#-----------------------------------------------------------------------------------

if __name__ == "__main__":
    main()

#-----------------------------------------------------------------------------------
