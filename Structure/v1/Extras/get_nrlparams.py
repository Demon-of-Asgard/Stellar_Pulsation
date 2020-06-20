#---------------------------------------------------------------------

def get_nrlparams(sp):

    c = sp.c
    pi = sp.pi
    hbar = sp.hbar
    Ms = sp.Ms
    R0 = sp.R0
    A = sp.A
    Z = sp.Z
    mn = sp.mn
    me = sp.me
    gamma_nrl = sp.gamma_nrl
    alpha_nrl = sp.alpha_nrl

    Knrl = (hbar**2/(15.0*pi**2*me))*(3.0*pi**2*Z/(A*mn*c**2))**gamma_nrl
    e0_nrl = ((R0/alpha_nrl)**gamma_nrl/Knrl)**(1/(gamma_nrl-1))
    barKnrl = Knrl*e0_nrl**(gamma_nrl-1)
    beta_nrl = 1.0e15*(4.*pi*e0_nrl)/(Ms*c**2*barKnrl**(1/gamma_nrl))

    params = {
        "alpha": alpha_nrl,
        "ep0": e0_nrl,
        "beta": beta_nrl,
        "K": Knrl,
        "barK": barKnrl,
        "gamma": gamma_nrl
    }

    return params

#---------------------------------------------------------------------
