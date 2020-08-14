#-------------------------------------------------------------------------------------------------

import numpy as np 
import matplotlib.pyplot as plt 

#-------------------------------------------------------------------------------------------------

from EoS.Eos import EoS_Asym as EoS

#-------------------------------------------------------------------------------------------------

class Plot:

    def __init__(self):
        pass

    #-------------------------------------------------------------------------------------------------


    def plot(self, u_list, barE, barP):

        plt.rc('text', usetex=True)
        plt.style.use('seaborn-whitegrid')

        lw = 1.5
        ls = '-'

        color = {
            "indigo": "#4B0082",
            "royal_blue": "#4169E1",
            "crimson": "#B22222",
            "salmon": "#FA8072",
        }

        plt.plot(
            u_list, barE, color = color["royal_blue"], linestyle = ls, 
            linewidth = lw, label = "E"
            )

        plt.plot(
            u_list, barP, color = color["salmon"], linestyle = ls, 
            linewidth = lw, label = "P"
            )
        plt.xlabel("u")
        plt.ylabel("[E/$\\varepsilon_0$]" +" and "+"[P/$\\varepsilon_0$]")

        plt.legend()
        plt.show()

#-------------------------------------------------------------------------------------------------

def main():
    
    u_list = np.arange(0.010,10.001, 0.001)
    _EoS_obj = EoS()

    barE = []
    barP = []

    for u in u_list:

        try:
            barE.append(_EoS_obj.E(u))
            barP.append(_EoS_obj.P(u))

        except ZeroDivisionError:
            print("'ZeroDivisionError' occured for u: {}".format(u))

    _plot_obj = Plot() 
    _plot_obj.plot(u_list,barE,barP)

#-------------------------------------------------------------------------------------------------

if __name__ == "__main__":
    main()

#-------------------------------------------------------------------------------------------------
