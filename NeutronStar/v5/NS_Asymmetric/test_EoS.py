#-------------------------------------------------------------------------------------------------

import numpy as np 
import matplotlib.pyplot as plt 

#-------------------------------------------------------------------------------------------------

from EoS.eos import Eos_Asym as EoS
from StarParams.params import Params as Prm
from StarParams.Unit_conversions import UnitConverter as Uconv 

#-------------------------------------------------------------------------------------------------

def plot(u_list, barE, barP):

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
        linewidth = lw, label = "E/$\\varepsilon_0$"
        )

    plt.plot(
        u_list, barP, color = color["salmon"], linestyle = ls, 
        linewidth = lw, label = "P/$\\varepsilon_0$"
        )
    plt.xlabel("u")
    plt.ylabel("[E/$\\varepsilon_0$]" +" and "+"[P/$\\varepsilon_0$]")

    plt.legend()
    plt.show()

#-------------------------------------------------------------------------------------------------

def main():
    
    u_list = np.arange(0.0,10.0, 0.001)
    _EoS_obj = EoS()
    _prm = Prm()
    
    star_params = _prm.get_star_params()
    e0 = star_params["e0"]

    _converter = Uconv()
    dict_conv_factors = _converter.get_conv_factors()
    ERG_TO_MEV = dict_conv_factors["ERG_TO_MEV"]
    e0MEV4 = e0/(ERG_TO_MEV)**4

    barE = []
    barP = []

    for u in u_list:
        barE.append(_EoS_obj.barE(u))
        barP.append(_EoS_obj.barP(u))
        
    plot(u_list,barE,barP)

#-------------------------------------------------------------------------------------------------

if __name__ == "__main__":
    main()

#-------------------------------------------------------------------------------------------------
