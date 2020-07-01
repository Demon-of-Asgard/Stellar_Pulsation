try:
    import numpy as np
except:
    print("ImportError-Numpy in 'build_star'")
    exit(1)
    
#----------------------------------------------------------------------------

def P_func(barP,barM,r,params):

    alpha = params["alpha"]
    gamma = params["gamma"]
    barK = params["barK"]
    Ms = params["Ms"]
    eps0 = params["eps0"]
    c = params["c"]
    pi = np.pi
    R0 = params["R0"]

    f1 = -alpha*barP**(1./gamma)*barM/r**2
    f2 = 1.0 + barK**(1.0/gamma)*barP**(1.0-1.0/gamma)
    f3 = 1.0 + (4.0*pi*eps0/(Ms*c**2))*(barP*r**3/(barM+1.0e-50))
    f4 = 1 - (2.0*R0*barM/r)

    delta_barPnxt = f1*f2*f3/f4

    return delta_barPnxt

#--------------------------------------------------------------------------------

def M_func(barP, barM, r, params):

    beta = params["beta"]
    gamma = params["gamma"]

    delta_barMnxt = beta*r**2*barP**(1/gamma)

    return delta_barMnxt
    
#--------------------------------------------------------------------------------

def build_star(params, outfname): 
    
    print("\n PARAMS:")
    print(" ---------------------------------\n")
    for entity in params:
        print(" > [",entity, ":", params[entity],"]\n")
    print("\n ---------------------------------\n")

    barP0 = params["barP0"]

    deltar = 1.0e-4
    #initial values.
    barP = [barP0]
    barM = [0.0]
    r = [0.0001]
    i = 0
    shoud_iter = True 

    while(shoud_iter == True):
        i += 1
        barP_i = barP[i-1] + deltar*P_func(barP[i-1], barM[i-1], r[i-1], params)
        barM_i = barM[i-1] + deltar*M_func(barP[i-1], barM[i-1], r[i-1], params)
        r_i = r[i-1] + deltar

        if(barP_i >= 0.):
            barP.append(barP_i)
            barM.append(barM_i)
            r.append(r_i)
        else:
            i -= 1
            shoud_iter = False 
    
    print()
    print(" Last entry: <{:3.0f}> {:5.4e}\t{:5.4e}\t{:5.4e}".format(i, r[i], barM[i], barP[i]))
    #barP = np.array(barP)*eps0
    output = np.array([r,barM,barP])
    np.savetxt(outfname,output.T,fmt="%6.5e",delimiter="\t\t")

    return outfname

#--------------------------------------------------------------------------------
