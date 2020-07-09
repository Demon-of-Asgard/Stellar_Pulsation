
#----------------------------------------------------------------------------------------

try:
    import numpy as np 
except ImportError:
    print("Numpy-ImportError. Please install numpy to the current venv. Thank you!")

from EoS.eos import eqn_state # Reminder: The class eqn_state inherits the class params

#----------------------------------------------------------------------------------------

class NS(eqn_state):

    def __init__(self):
        super().__init__()
    
    #----------------------------------------------------------------------------------------

    def Pfunc(self, P_bar, M_bar, r):
        params = self.get_params()

        R0 = params["R0"]

        e_bar = self.eos(P_bar)

        dP_dr = -R0*M_bar*e_bar/(r**2.0)

        return dP_dr

    #----------------------------------------------------------------------------------------

    def Mfunc(self, P_bar, M_bar, r):

        params = self.get_params()

        pi = params["pi"]
        Ms = params["Ms"]
        c  = params["c"]
        e0 = params["e0"]

        beta  = (4.0*pi*e0)/(Ms*c**2)

        e_bar = self.eos(P_bar)

        dM_dr = beta*e_bar/(r**2.0)

        return dM_dr

    #----------------------------------------------------------------------------------------

    def construct(self):

        loop_again = True
        params = self.get_params()
        P_bar = [params["P0"]]
        M_bar = [params["M0"]]
        dr = 1.0
        R = [0.01] 
        i = 0

        while(loop_again):
            
            r = R[i]
            p_bar =  P_bar[i] + self.Pfunc(P_bar[i], M_bar[i], r)*dr
            m_bar = M_bar[i] + self.Mfunc(P_bar[i], M_bar[i], r)*dr
            r = r+dr
            #print(" > i : {} P_bar: {}".format(i, p_bar))
            if (p_bar>=0):
                P_bar.append(p_bar)
                M_bar.append(m_bar)
                R.append(r)
                print(" >", i, r, p_bar)
                i = i+1
            else:
                loop_again = False

        output = np.array([R, M_bar, P_bar]).T
        np.savetxt("output.dat",output)
        return "Success"

    #----------------------------------------------------------------------------------------

def main():
    Star = NS()
    Star.construct()
            
if __name__ == "__main__":
    main()

#----------------------------------------------------------------------------------------
