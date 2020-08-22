#!/home/demon/anaconda3/bin/python3

import numpy as np
import matplotlib.pyplot as plt
import star_params as sp


#Constants.
pi = np.pi
alpha = 0.05
c = 3.e10
Ms = 1.989e33 #[gm]
eps0 = 7.463e39
beta = 0.005924
gamma = 5.0/3.0
delta_r = 1.0e0
centP = 1.0e-4
R0 = 1.47



def fP(M_i,P_i,r_i):

    barK, eps0, beta = sp.get_params("rl", alpha)
    f1 = -alpha*P_i**(1/gamma)*M_i/r_i**2
    f2 = 1.0 + barK**(1/gamma)*P_i**(1-(1/gamma))
    f3 = 1+(4.*pi*eps0/(Ms*c**2))*(P_i*r_i**3/(M_i+0.00001))
    f4 = 1-R0*M_i/r_i**2
    f = f1*f2*f3/f4
    return f

def fM(M_i,P_i,r_i):

    f = beta*r_i**2*P_i**(1/gamma)
    return f

def make_star():

    i = 0
    barM = [0.0]
    barP = [centP]
    r = [0.000001]
    
    while True:
        
        k1M = delta_r*fM(barM[i],barP[i],r[i])
        k2M = delta_r*fM(barM[i]+k1M/2.0,barP[i],r[i]+delta_r/2.0)
        k3M = delta_r*fM(barM[i]+k2M/2.0,barP[i],r[i]+delta_r/2.0)
        k4M = delta_r*fM(barM[i]+k3M/2.0,barP[i],r[i]+delta_r/2.0)


        if barP[i]>0.0:
            k1P = delta_r*fP(barM[i],barP[i],r[i])
        else:
            break
        if barP[i]+k1P/2.0 >0:
            k2P = delta_r*fP(barM[i],barP[i]+k1P/2.0,r[i]+delta_r/2.0)
        else:
            break
        if barP[i]+k2P/2.0:
            k3P = delta_r*fP(barM[i],barP[i]+k2P/2.0,r[i]+delta_r/2.0)
        else:
            break
        if barP[i]+k3P/2.0:
            k4P = delta_r*fP(barM[i],barP[i]+k3P/2.0,r[i]+delta_r/2.0)
        else:
            break


        barM_nxt = barM[i] + (1.0/6.0)*(k1M+2.0*k2M+2.0*k3M+k4M)
        barP_nxt = barP[i] + (1.0/6.0)*(k1P+2.0*k2P+2.0*k3P+k4P)

        #print("r: {0:8.5e} Mr: {1:8.5e} Pr: {2:8.5e}".format(r[i],barM_nxt,barP_nxt))

        if barP_nxt > 0.0:
            barM.append(barM_nxt)
            barP.append(barP_nxt)
            r.append(r[i]+delta_r)
            i += 1
        
        else:
            break

        if i>10000000:
            break
        else:
            pass
    
    print("Iterations = {}".format(i))
    print("r:",r[-1]," ","M:",barM[-1])
    outfile = "output.dat"
    out = np.array([r,barM,barP]).T
    np.savetxt(outfile,out, fmt = '%.8e', newline = '\n', delimiter = '\t')
    return(outfile)
    



def plot(outfile):

    linestyle = {"densly_dashdotted":  (0, (3, 1, 1, 1)),}
    dataE = np.genfromtxt(outfile)
    fig,ax = plt.subplots(1,2)
    ax[0].plot(dataE[:,0]/1000,dataE[:,1],color="red", linestyle="-")
    ax[0].set_xlabel("r [$10^3$ km]")
    ax[1].plot(dataE[:,0]/1000,dataE[:,2], color = "red", linestyle = linestyle["densly_dashdotted"])
    ax[1].set_xlabel("r [$10^3$ km]")
    plt.show()

def main():

    outfile = make_star()
    plot(outfile)

if __name__ == "__main__":
    main()
