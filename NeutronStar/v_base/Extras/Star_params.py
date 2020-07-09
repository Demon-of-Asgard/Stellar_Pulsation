#!/home/demon/anaconda3/bin/python3
#---------------------------------------------------------------------
import numpy as np
#---------------------------------------------------------------------
def get_params(type):
    pi= np.pi
    c = 3.0e10 #[cm/s]
    Ms = 1.989e33 #[gm]
    h = 6.62607004e-27 #[erg s]
    hbar = h/(2.0*pi) #[erg s]
    mn = 1.6749e-24 #[gm]
    me = 9.10e-28 # [gm]
    A =2.1452
    Z = 1
    R0 = 1.476 #[km]
    alpha_rl = R0
    alpha_nrl = 0.05e3 #[cm]
    gamma_rl = 4./3.
    gamma_nrl = 5./3.
    #---------------------------------------------------------------------
    if type == "rl":
        Krl = (hbar*c/(12.*pi**2))*(3.0*pi**2*Z/(A*mn*c**2))**gamma_rl
        e0_rl = ((R0/alpha_rl)**gamma_rl/Krl)**(1/(gamma_rl-1))
        barKrl = Krl*e0_rl**(gamma_rl-1)
        beta_rl = 1.0e15*(4.0*pi*e0_rl)/(Ms*c**2*barKrl**(1/gamma_rl))
        params = [e0_rl, beta_rl]
        return params
        #---------------------------------------------------------------------

    if type == "nrl":
        Knrl = (hbar**2/(15.0*pi**2*me))*(3.0*pi**2*Z/(A*mn*c**2))**gamma_nrl
        e0_nrl = ((R0/alpha_nrl)**gamma_nrl/Knrl)**(1/(gamma_nrl-1))
        barKnrl = Knrl*e0_nrl**(gamma_nrl-1)
        beta_nrl = 1.0e15*(4.*pi*e0_nrl)/(Ms*c**2*barKnrl**(1/gamma_nrl))
        params = [e0_nrl, beta_nrl]
        return params
        #---------------------------------------------------------------------

def main():
    type = int(input("Enter the type. \n [1] Relativistic star \n [2] Non-relativistic star \n >"))
    if type == 1:
        params = get_params("rl")
        print("[e0, beta] = ", params)
    elif type == 2:
        params = get_params("nrl")
        print("[e0, beta] = ", params)
    else:
        print("Invalid input")

if __name__ == "__main__":
    main()
