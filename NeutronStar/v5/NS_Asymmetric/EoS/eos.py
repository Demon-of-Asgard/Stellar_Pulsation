try:
    from StarParams.params import Params as Prm 
except:
    print("Error wile importing 'StarParams.params' @ eos.py")
    exit(1)
#-----------------------------------------------------------------------------------


class Eos_Asym():

    def __init__(self):

        _prm = Prm()
        conv_factors = _prm.get_conv_factors()

        self.__mN = 1.6749e-24  # [gm]
        self.__A = 118.2 * conv_factors["MEV_TO_ERG"] # A[MeV] --> A[erg]
        self.__B = 65.39 * conv_factors["MEV_TO_ERG"] # B[MeV] --> B[erg]
        self.__sigma = 2.112 # dimension less
        self.__alpha = 1.0   # Dimension less

    #-----------------------------------------------------------------------------------

    def get_normalized_parameters(self):

        _prm = Prm()
        star_params = _prm.get_star_params()

        c   = star_params["c"]
        e0  = star_params["e0"] # [erg/cm^3]
        n0  = star_params["n0"] # [#/cm^3]
        EF0 = star_params["EF0"]

        eN0 = self.__mN*c**2

        barEN0 = eN0*(n0/e0)
        barA   = self.__A *(n0/e0)
        barB   = self.__B *(n0/e0)
        barEF0 = EF0*(n0/e0)
        sigma  = self.__sigma 
        alpha  = self.__alpha

        dict_normalized_coeffs = {
            "barEN0":barEN0,
            "barA":barA,
            "barB":barB,
            "barEF0":barEF0,
            "sigma":sigma,
            "alpha":alpha,
        }

        return dict_normalized_coeffs

    #-----------------------------------------------------------------------------------

    def barE(self,u):

        eos_params = self.get_normalized_parameters()

        barEN0 = eos_params["barEN0"]
        barA   = eos_params["barA"]
        barB   = eos_params["barB"]
        sigma  = eos_params["sigma"]
        barEF0 = eos_params["barEF0"]
        alpha  = eos_params["alpha"]
        
        '''
        barA = -118.2
        barB = 65.39
        sigma = 2.112
        barEF0 = 22.1
        barEN0 = 938.0
        '''
    

        barE = (barEN0*u) + (barEF0*u**(5.0/3.0)) + ((barA/2.0)*u**2) + ((barB/(1+sigma))*u**(1.0+sigma)) 

        return barE
    #-----------------------------------------------------------------------------------

    def barP(self, u):

        eos_params = self.get_normalized_parameters()

        barA   = eos_params["barA"]
        barB   = eos_params["barB"]
        sigma  = eos_params["sigma"]
        barEF0 = eos_params["barEF0"]
        alpha  = eos_params["alpha"]

        '''
        barA = -118.2
        barB = 65.39
        sigma = 2.112
        barEF0 = 22.1
        '''


        barP = ((2.0/3.0)*barEF0)*u**(5.0/3.0) + (barA/2.0)*u**2 + (barB*sigma/(1.0+sigma))*u**(1.0+sigma)

        return barP 

    #-----------------------------------------------------------------------------------

    def dbarP_du(self, u):
        pass 
    
    #-----------------------------------------------------------------------------------

    def eos(self, u):
        _prm = Prm()
        star_params = _prm.get_star_params()

        barE = 0
        return barE

#-----------------------------------------------------------------------------------

def main():
    pass
if __name__ == "__main__":
    main()