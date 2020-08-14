import numpy as np 

#-----------------------------------------------------------------------------------

class EoS_Asym:

    def __init__(self):

        self.__mN = 939.57 # MeV
        self.__A  = -118.2 # MeV
        self.__B  = 65.39 # MeV
        self.__sigma = 2.112 # dimension less
        self.__alpha = 1.0   # Dimension less
        self.__S0  = 30.0
        self.__EF0 =  22.1 # MeV
        self.__n0  = 0.16#122.346e4 # MeV^4

    #-----------------------------------------------------------------------------------

    def F(self, u):
        return u
    
    #-----------------------------------------------------------------------------------
    
    def dF_du(self, u):
        return 1.0

    #-----------------------------------------------------------------------------------

    def S(self, u):

        S = (
            (2.0**(2.0/3.0)-1)*(5.0/3.0)*self.__EF0*(u**(2.0/3.0)-self.F(u))
            + self.__S0*self.F(u)
            )

        return S
    
    #-----------------------------------------------------------------------------------

    def dS_du(self, u):

        EF0 = self.__EF0
        S0 = self.__S0
        
        dS_du = (
            (2.0**(2.0/3.0))*(5.0/3.0)*EF0*((2.0/3.0)*u**(-1.0/3.0)-self.dF_du(u))
            + S0*self.dF_du(u)
        )

        return dS_du

    #-----------------------------------------------------------------------------------
    

    def E(self, u):

        mN = self.__mN
        n0  = self.__n0
        A   = self.__A
        B   = self.__B 
        EF0 = self.__EF0
        sigma = self.__sigma 
        alpha = self.__alpha 

        E = n0*(
            u*mN
            + (EF0*u**(5.0/3.0)) 
            + ((A/2.0)*u**2) 
            + ((B/(1+sigma))*u**(1.0+sigma)) 
            + alpha**2 * self.S(u)
            )
        
        return E 

    #-----------------------------------------------------------------------------------

    def P(self, u):

        n0  = self.__n0
        EF0 = self.__EF0
        A   = self.__A
        B   = self.__B 
        sigma = self.__sigma 
        alpha = self.__alpha 

        P = n0*(
            (2.0/3.0)*EF0*u**(5.0/3.0)
            + (A/2.0)*u**2
            + (B*sigma/(1+sigma))*u**(1.0+sigma)
            + (alpha**2)*(u**2)*self.dS_du(u)
        )

        return P


