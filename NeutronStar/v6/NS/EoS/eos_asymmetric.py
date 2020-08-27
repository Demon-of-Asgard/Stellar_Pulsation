from StarParams.params import Parameters as Prm
from Root.root import RootFinder as Root
#-----------------------------------------------------------------------------------

class EoSAsym:
    def __init__(self):

        prmObj = Prm()
        conv_factors = prmObj.get_conv_factors()
        constants = prmObj.get_constants()

        '''Global constants.'''
        self.__mE = constants["mE"]
        self.__mP = constants["mP"] 
        self.__mN = constants["mN"] #[gm]
        self.__c  = constants["c"]  #[cm/s]
        self.__pi = constants["pi"]
        self.__hbar = constants["hbar"] # [erg . s]
        self.__eN0 = self.__mN*self.__c**2 # ([erg])Rest mass energy of neutron.


        '''EoS specific constants.'''
        self.__A = -118.2 * conv_factors["MEV_TO_ERG"] # A[MeV] --> A[erg]
        self.__B = 65.39 * conv_factors["MEV_TO_ERG"] # B[MeV] --> B[erg]
        self.__sigma = 2.112 # dimension less
        self.__S0 = 30.0* conv_factors["MEV_TO_ERG"] #[erg]
        self.__n0 = 0.16/(10**(-13))**3 # [#/(fm^3)] --> [#/(cm^3)]
        self.__k0 = (3.0*self.__pi**2.0*self.__hbar**3.0*self.__n0)**(1.0/3.0)
        self.__static_alpha = 1.0   # Dimension less
        self.__is_dynamic_alpha = False
        self.u_lower = 0.0 # u = (n/n0)
        self.u_upper = 100.0
        self.k_lower = -(1.0/10.0)*self.__k0*self.u_upper**(1.0/3.0)
        self.k_upper = self.__k0*self.u_upper**(1.0/3.0)


        self.__e0 = self.__mN**4*self.__c**5/(3.*self.__pi**2*self.__hbar**3) # [erg/cm^-3] # This quantity is arbitrary
        self.__EF0 = (
            (3.0/5.0)*(1.0/(2.0*self.__mN))
            *(3.0*self.__pi**2*self.__hbar**3*self.__n0/2.0)**(2.0/3.0)
            ) # [erg] # <KE> of a nucleon corresponding to k_F(n) = k_F(n0)

    #-----------------------------------------------------------------------------------

    def get_EoS_params(self):

        EoS_params = {
            "A": self.__A,
            "B":self.__B,
            "sigma": self.__sigma,
            "alpha": self.__static_alpha,
            "S0": self.__S0,
            "n0": self.__n0,
            "e0": self.__e0,
            "EF0": self.__EF0,
        }

        return EoS_params

    def print_EoS_params(self):
        ''' Only to check the values. '''

        left_width = 15
        rigt_width = 15
        EoS_params = self.get_EoS_params()

        print(" ","EOS PARAMETERS".center(left_width+rigt_width, '-'))
        print(" ","="*(left_width+rigt_width))

        for item, value in  EoS_params.items():
            print(" ",item.ljust(left_width, '.')+" | "+'{:8.6e}'.format(value).ljust(rigt_width))
        print("\n")

    #-----------------------------------------------------------------------------------

    def kP_of_kN(self, kN):

        kP_num = (
            (kN**2*self.__c**2 + self.__mN**2.0*self.__c**4 - self.__mE**2.0*self.__c**4.0)**2.0
            -(
                (2.0*self.__mP**2.0*self.__c**4.0)*
                (kN**2*self.__c**2 + self.__mN**2.0*self.__c**4 + self.__mE**2.0*self.__c**4.0)
            )
            +(self.__mP**4.0*self.__c**8.0)
        )**(1.0/2.0)

        kP_denom = (
            2.0*self.__c*(kN**2*self.__c**2 + self.__mN**2.0*self.__c**4)**(1.0/2.0)
        )

        kP = kP_num/kP_denom

        return kP

    #-----------------------------------------------------------------------------------

    def get_nP(self, kP):
        nP = kP**3.0/(3.0*self.__pi**2.0*self.__hbar**3.0)
        return nP

    #-----------------------------------------------------------------------------------

    def get_nN(self, kN):
        nN = kN**3.0/(3.0*self.__pi**2.0*self.__hbar**3.0)
        return nN
    #-----------------------------------------------------------------------------------
    
    def Ntot(self, kN, ntot_i = 0.0):
        
        kP = self.kP_of_kN(kN)
        nP = self.get_nP(kP)
        nN = self.get_nN(kN)
        ntot = nP + nN 

        return ntot - ntot_i

    #-----------------------------------------------------------------------------------

    def get_kN(self, ntot_curr):
        rootObj = Root()
        kN_curr = rootObj.Bisection(self.Ntot, ntot_curr, self.k_lower, self.k_upper)
        return kN_curr

    #-----------------------------------------------------------------------------------

    def get_alpha(self, ntot_curr):
        kN_curr = self.get_kN(ntot_curr)
        kP_curr = self.kP_of_kN(kN_curr)

        nP = self.get_nP(kP_curr)
        nN = self.get_nN(kN_curr)
        alpha = (nN-nP)/(nN+nP)

        return alpha

    #-----------------------------------------------------------------------------------

    def get_u(self,barP_i):
        rootObj = Root()
        u = rootObj.Bisection(self.barP, barP_i, self.u_lower, self.u_upper)
        return u

    #-----------------------------------------------------------------------------------

    def get_density(self, barP_i):
        u = self.get_u(barP_i)
        rho = u*self.__n0*self.__mN
        return rho
        
    #-----------------------------------------------------------------------------------

    def F(self, u):
        return u

    #-----------------------------------------------------------------------------------

    def dF_du(self, u):
        return 1.0

    #-----------------------------------------------------------------------------------

    def S(self, u):

        S = (
            (2.0**(2.0/3.0)-1)*self.__EF0*(u**(2.0/3.0)-self.F(u))
            + self.__S0*self.F(u)
        )

        return S

    #-----------------------------------------------------------------------------------
    def dS2_du(self, u): # = u^2 * dS/du

        dS2_du = (
            (2.0**(2.0/3.0)-1)*self.__EF0
            *((2.0/3.0)*u**(5.0/3.0)-u**2*self.dF_du(u))
            + u**2*self.__S0*self.dF_du(u)
        )

        return dS2_du

    #-----------------------------------------------------------------------------------

    def barP(self,u, barP_i=0.0):

        if self.__is_dynamic_alpha == True:
            n_curr = self.__n0*u
            alpha = self.get_alpha(n_curr)
        else:
            alpha = self.__static_alpha

        P = self.__n0*(
            (2.0/3.0)*self.__EF0*u**(5.0/3.0)
            + (self.__A/2.0)*u**2
            + (self.__B*self.__sigma/(1+self.__sigma))*u**(1.0+self.__sigma)
            + (alpha**2)*self.dS2_du(u)
        )

        barP = (P/self.__e0)

        return (barP-barP_i)


    #-----------------------------------------------------------------------------------

    def barE(self, barP_i):

        u = self.get_u(barP_i)
        if self.__is_dynamic_alpha == True:
            n_curr = self.__n0*u
            alpha = self.get_alpha(n_curr)
        else:
            alpha = self.__static_alpha
        #print("alpha: {}".format(alpha))

        E = self.__n0*(
            (self.__eN0*u)
            + (self.__EF0*u**(5.0/3.0))
            + ((self.__A/2.0)*u**2)
            + ((self.__B/(1+self.__sigma))*u**(1.0+self.__sigma))
            + alpha**2.0*u*self.S(u)
        )

        barE = (E/self.__e0)

        return barE


#-----------------------------------------------------------------------------------

def main():
    pass
if __name__ == "__main__":
    main()
