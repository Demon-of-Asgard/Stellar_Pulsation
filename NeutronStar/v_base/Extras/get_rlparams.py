#---------------------------------------------------------------------

def get_rlparams(sp):

    c = sp.c
    pi = sp.pi
    hbar = sp.hbar
    Ms = sp.Ms
    R0 = sp.R0
    A = sp.A
    Z = sp.Z
    mn = sp.mn
    gamma_rl = sp.gamma_rl
    alpha_rl = sp.alpha_rl

    Krl = (hbar*c/(12.*pi**2))*(3.0*pi**2*Z/(A*mn*c**2))**gamma_rl
    e0_rl = ((R0/alpha_rl)**gamma_rl/Krl)**(1/(gamma_rl-1))
    barKrl = Krl*e0_rl**(gamma_rl-1)
    beta_rl = 1.0e15*(4.0*pi*e0_rl)/(Ms*c**2*barKrl**(1/gamma_rl))

    params = {
        "alpha":alpha_rl,
        "ep0":e0_rl, 
        "beta":beta_rl,
        "K":Krl,
        "barK":barKrl,
        "gamma":gamma_rl
        }
    print(params)
    return params

#---------------------------------------------------------------------
