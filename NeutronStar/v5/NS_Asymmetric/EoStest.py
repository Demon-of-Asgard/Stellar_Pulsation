import numpy as np 
import matplotlib.pyplot as plt 
from EoS.EoS_Asym import EoS_Asym as EoS
from Root_finding.Rootfinding import RootFinder as Root 


def main():

    u = np.arange(0.0, 2.0, 1.0e-6)
    
    eosObj =EoS()
    eosObj.print_EoS_params()
    barE = []
    barP = []
    for ui in u:
        barP.append(eosObj.barP(ui,0.0))
        barE.append(eosObj.barE(ui))
    
    plt.xscale("log")
    plt.yscale("log")
    plt.plot(barP,barE, "-")
    #plt.plot(u, barE, "-")

    plt.show()

if __name__ == "__main__":
    main()