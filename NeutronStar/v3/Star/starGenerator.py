#!/home/demon/anaconda3/bin/python3

from StarParams.params import Params as Prm
import numpy as np

class NS(Prm):

    barP = []
    barM = []
    r = []

    def __init__(self):
        super().__init__()
        init_params = self.get_init_params()
        self.barP.append(init_params["barP0"])
        self.barM.append(init_params["barM0"])
        self.r.append(init_params["r0"])
        

    def print_params(self):
        params = self.get_params()
        for index in params:
            print(" > [{}:={:5.4e}]\n".format(index, params[index]))
    
    def dbarM_dr(self):
        params = self.get_params()
        

    def dbarP_dr(self):
        pass

    def build(self):
        i = 0
        dr = 1.0
        while(True and i<= 10000):
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
    
    outfname = "output.dat"
    output = np.array([r,barM, barP]).T
    np.savetxt(outfname,output, fmt='%.18e', delimiter="\t")


def main():
    star_obj = NS()
    star_obj.build()

if __name__ == "__main__":
    main()
