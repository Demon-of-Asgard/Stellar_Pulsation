import numpy as np 
import matplotlib.pyplot as plt 
from EoS.EoS_Asym import EoS_Asym as EoS
from Root_finding.Rootfinding import RootFinder as Root 


def main():

    u = np.arange(0.0, 10.0, 0.50)
    
    eosObj =EoS()
    eosObj.print_EoS_params()
    barE = []
    barP = []
    for ui in u:
        barP.append(eosObj.barP(ui,0.0))
        barE.append(eosObj.barE(barP[-1]))
    
    # plt.xscale("log")
    # plt.yscale("log")
    plt.plot(u,barE, "-")
    plt.plot(u, barP, "-")

    plt.show()

if __name__ == "__main__":
    main()