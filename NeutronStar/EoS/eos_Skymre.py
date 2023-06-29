#-----------------------------------------------------------------------------------

from StarParams.params import Parameters as Prm 
from Root.root import RootFinder as Root

#-----------------------------------------------------------------------------------

class EoSSkymre:
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
        self.__eN0  = self.__mN*self.__c**2 # ([erg])Rest mass energy of neutron.

        '''EoS specific constants.'''
        self.__t0 = 1024.1 * conv_factors["MEV_TO_ERG"] * (10**-13)**3
        self.__t3 = 14600.8 * conv_factors["MEV_TO_ERG"] * (10**-13)**6
        self.__n0 = 0.16/(10**(-13))**3 # [#/(fm^3)] --> [#/(cm^3)]
        self.__e0 = self.__mN**4*self.__c**5/(3.*self.__pi**2*self.__hbar**3) # [erg/cm^-3] # This quantity is arbitrary
        self.__EF0 = (
            (3.0/5.0)*(1.0/(2.0*self.__mN))
            *(3.0*self.__pi**2*self.__hbar**3*self.__n0/2.0)**(2.0/3.0)
            )
        self.u_lower = 0.0 # u = (n/n0)
        self.u_upper = 200.0

    #-----------------------------------------------------------------------------------
    def get_EoS_params(self):
        EoS_params = {
            "t0": self.__t0,
            "t3": self.__t3,
            "n0": self.__n0,
            "e0": self.__e0,
        }
        return EoS_params
    #-----------------------------------------------------------------------------------
    
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

    def get_u(self, barP_i):
        rootObj = Root()
        u = rootObj.Bisection(self.barP, barP_i, self.u_lower, self.u_upper)
        return u

    #-----------------------------------------------------------------------------------
    
    def get_density(self, barP_i):
        u = self.get_u(barP_i)
        rho = u*self.__n0*self.__mN
        return rho
        
    #-----------------------------------------------------------------------------------

    def barP(self, u,  barP_i=0.0):

        P = self.__n0*(
            ((2.0/3.0)*2.0**(2.0/3.0)*self.__EF0)*u**(5.0/3.0)
            -((1.0/4.0)*self.__t0*self.__n0)*u**2.0
            +((1.0/12.0)*self.__t3*self.__n0**2.0)*u**3.0
        )
        
        barP = P/self.__e0 

        return (barP - barP_i)

    #-----------------------------------------------------------------------------------

    def barE(self, barP_i):
        u = self.get_u(barP_i)
        
        E = self.__n0*(
            (self.__eN0)*u
            +(2.0**(2.0/3.0)*self.__EF0)*u**(5.0/3.0)
            -((1.0/4.0)*self.__t0*self.__n0)*u**2.0
            +((1.0/24.0)*self.__t3*self.__n0**2.0)*u**3.0
        )

        barE = (E/self.__e0)

        return barE

#-----------------------------------------------------------------------------------
