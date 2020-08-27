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


    def plot(self, u_list, E, P):

        plt.rc('text', usetex=True)
        plt.style.use('seaborn-whitegrid')

        lw = 2.0
        ls = '-'

        color = {
            "indigo": "#4B0082",
            "royal_blue": "#4169E1",
            "crimson": "#B22222",
            "salmon": "#FA8072",
        }

        # plt.plot(
        #     u_list, E, color = color["royal_blue"], linestyle = ls,
        #     linewidth = lw, label = "E"
        #     )

        plt.plot(
            u_list, P, color = color["salmon"], linestyle = ls,
            linewidth = lw, label = "P"
            )
        plt.xlabel("u")
        plt.ylabel("E, P$~[MeV/fm^3]$")

        plt.legend()
        plt.show()

        #-------------------------------------------------------------------------------------------------

        plt.plot(
            E, P, color = color["salmon"], linestyle = ls,
            linewidth = lw, label = "P"
            )
        plt.xlabel("E$~[MeV/fm^3]$")
        plt.ylabel("P$~[MeV/fm^3]$")

        plt.legend()
        plt.show()
#-------------------------------------------------------------------------------------------------

def main():

    u_list = np.arange(0.010,2.001, 0.0001)
    _EoS_obj = EoS()

    E = []
    P = []

    for u in u_list:

        try:
            E.append(_EoS_obj.E(u))
            P.append(_EoS_obj.P(u))

        except ZeroDivisionError:
            print("'ZeroDivisionError' occured for u: {}".format(u))

    _plot_obj = Plot()
    _plot_obj.plot(u_list, E, P)

#-------------------------------------------------------------------------------------------------

if __name__ == "__main__":
    main()

#-------------------------------------------------------------------------------------------------
