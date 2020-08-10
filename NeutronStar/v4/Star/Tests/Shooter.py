import numpy as np 
'''
Program to find the value ou for a given pressure. 
'''

def BE_fun(u):
    A =-122.2
    B = 65.39
    sigma = 2.112
    mN = 939.56
    E0 = 22.1

    e_by_n = mN + E0*u**(2.0/3.0) + (A/2.0)*u + (B/(1.0+sigma))*u**sigma 
    BE = e_by_n - mN

    return BE


def P_fun(u):
    A = -118.2
    B = 65.39
    sigma = 2.112
    E0 = 22.1

    P =  (2.0/3.0)*E0*u**(5.0/3.0) + (A/2.0)*(u**2) + (B*sigma/(1.0+sigma))*(u**(1+sigma))

    return P



def get_u(be_n):

    u = 0
    should_loop = True
    i = 0
    while should_loop:

        i = i+1
        diff = (be_n - BE_fun(u))/(1.0e-100+be_n)
        
        if  abs(diff) <= 1.0e-10:
            print ("#Iterations: ", i, u,be_n, diff)
            return u
        else:
            u = u+0.001
            #print("#Iterations: ", i, u, be_n, diff)





def Plot():
    import matplotlib.pyplot as plt 
    plt.style.use('seaborn-whitegrid')
    color = {
        "indigo": "#4B0082",
        "royal_blue": "#4169E1",
        "crimson": "#B22222",
        "salmon": "#FA8072",
    }
    u = np.arange(0.0,2, 0.01)
    BE = np.array([BE_fun(uu) for uu in u])
    P = np.array([P_fun(uu) for uu in u])
    P_prime = np.array([P_fun(get_u(be_n)) for be_n in BE])
    y0 = [0 for uu in u]
    print("Plotting")
    plt.plot(u, BE, label = "BE")
    plt.plot(u, y0,"k-.", lw = 1.0)
    plt.show()
    plt.plot(u, P, label = "P", zorder = 1)
    plt.scatter(u, P_prime, marker=".", color = color ["crimson"], linewidth = 0.025, label="P'", zorder = 2)
    plt.plot(u, y0, "k-.", lw= 1.0 )
    plt.legend(loc = 0)
    plt.show()
    

if __name__ == "__main__":
    Plot()
