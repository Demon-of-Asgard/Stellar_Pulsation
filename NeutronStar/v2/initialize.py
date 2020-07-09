import numpy as np

from Params.params import params as prms 

obj = prms()
params = obj.get_params()
P0 = params["e0"]

def P(x):
    e0 = params["e0"]
    p = e0*((2.*x**3-3.0*x)*(1.0+x**2)**(1.0/2.0)+3.0*np.arcsinh(x))
    return p


def init_kF():
    x = 0.0
    
    while(P0-P(x)>=0.0):
        print(" > x:{:5.3e} \t P:{:5.3e}".format(x,P(x)))
        x = x+0.001
        
    print("\n > x: {:5.3e}".format(x))


if __name__ == "__main__":
    init_kF()
    