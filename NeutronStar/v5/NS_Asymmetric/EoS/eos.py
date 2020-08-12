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
        _prm.print_conv_factors()

        self.__A = 118.2 * conv_factors["MEV_TO_ERG"] # A[MeV] -- > A[erg]
        self.__B = 65.39 * conv_factors["MEV_TO_ERG"] # B[MeV] --> B[erg]
        self.__sigma = 2.112 # dimension less
        self.__alpha = 1.0 # Dimension less

    #-----------------------------------------------------------------------------------

    def get_normalized_coefficiants(self):

        _prm = Prm()
        star_params = _prm.get_star_params()

        e0 = star_params["e0"] # [erg/cm^3]
        n0 = star_params["n0"] # [#/cm^3]
        EF0 = star_params["EF0"]

        barA = self.__A *(n0/e0)
        barB = self.__B *(n0/e0)
        barEF0 = EF0*(n0/e0)
        sigma = self.__sigma 
        alpha = self.__alpha

        dict_normalized_coeffs = {
            "barA":barA,
            "barB":barB,
            "barEF0":barEF0,
            "sigma":sigma,
            "alpha":alpha,
        }

        return dict_normalized_coeffs
    #-----------------------------------------------------------------------------------



    def barE(self,u):
        pass


    def barP(self, u):
        pass

    def dbarP_du(self, u):
        pass 
    

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