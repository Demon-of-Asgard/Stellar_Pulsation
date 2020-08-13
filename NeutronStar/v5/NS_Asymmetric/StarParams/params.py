try:
    import numpy as np
except ImportError:
    print("ImportError-Numpy @ params")
    exit(1)

class Params:
    def __init__(self):
        
        '''Primary Constants'''
        self.__pi = np.pi
        self.__c  = 2.9979e10  # [cm/s]
        self.__Ms = 1.989e33  # [gm]
        self.__h  = 6.62607004e-27  # [erg s]
        self.__hbar = self.__h/(2.0*self.__pi)  # [erg s]
        self.__MeV_to_Erg =  1.602176487E-6 # Dimesion less
        self.__mN = 1.6749e-24  # [gm]
        self.__me = 9.10e-28  # [gm]
        self.__R0 = 1.476  # [km]
        self.__n0 = 1.60e38 # [cm^-3] := n0 = 0.16/(fm^3)

        '''Conversion factors'''
        self.__INV_MEV_TO_CM = (1.0/self.__MeV_to_Erg)*self.__hbar*self.__c
        self.__INV_MEV_TO_SEC = (1.0/self.__MeV_to_Erg)*self.__hbar
        self.__MEV_TO_ERG = self.__MeV_to_Erg
        self.__MEV_TO_GM = self.__MeV_to_Erg/(self.__c**2)

        #Check later 
        '''Derived constants'''
        self.__e0N = self.__mN*self.__c**2 # Rest mass energy of neutron.
        self.__e0 = self.__mN**4*self.__c**5/(3.*self.__pi**2*self.__hbar**3) # This quantity is arbitrary
        self.__EF0 = (3.0/5.0)*(1.0/(2.0*self.__mN))*(3.0*self.__pi**2*self.__hbar**3*self.__n0/2.0)**(2.0/3.0) # <KE> of neutron for k_F(n0)

        '''Initial values'''
        self.__barM0 = 0.0
        self.__barP0 = 0.01
        self.__dr = 1.0e-4
        self.__r0 = self.__dr

    #-------------------------------------------------------------

    def get_conv_factors(self):

        dict_conv_factors = {
            "INV_MEV_TO_CM":self.__INV_MEV_TO_CM,
            "INV_MEV_TO_SEC":self.__INV_MEV_TO_SEC,
            "MEV_TO_ERG":self.__MEV_TO_ERG,
            "MEV_TO_GM":self.__MEV_TO_GM,

            "CM_TO_INV_MEV":1.0/self.__INV_MEV_TO_CM,
            "SEC_TO_INV_MEV": 1.0/self.__INV_MEV_TO_SEC,
            "ERG_TO_MEV":1.0/self.__MEV_TO_ERG,
            "GM_TO_MEV":1.0/self.__MEV_TO_GM,
        }
        return dict_conv_factors

    #-------------------------------------------------------------
    
    def print_conv_factors(self):

        left_width = 15
        rigt_width = 15

        conversion_factors = self.get_conv_factors()

        print(" ","CONVERSION FACTORS".center(left_width+rigt_width, '-'))
        print(" ","="*(left_width+rigt_width))

        for item, value in  conversion_factors.items():
            print(" ",item.ljust(left_width, '.')+" | "+'{:8.6e}'.format(value).ljust(rigt_width))
            #print(" "+"_"*(left_width+rigt_width+4))
        print("\n")

    #-----------------------------------------------------------

    def get_star_params(self):

        star_params = {
            "c":self.__c,
            "pi":self.__pi,
            "R0":self.__R0,
            "Ms":self.__Ms,
            "n0":self.__n0,
            "e0N":self.__e0N,
            "e0" : self.__e0,
            "EF0":self.__EF0,
        }

        return star_params

    #-------------------------------------------------------------
    
    def print_star_params(self):
        ''' Only to check the values. '''

        left_width = 15
        rigt_width = 15
        star_params = self.get_star_params()

        print(" ","STAR PARAMS".center(left_width+rigt_width, '-'))
        print(" ","="*(left_width+rigt_width))

        for item, value in  star_params.items():
            print(" ",item.ljust(left_width, '.')+" | "+'{:8.6e}'.format(value).ljust(rigt_width))
        print("\n")

    #-------------------------------------------------------------
    
    def get_star_init_values(self):

        r0 = self.__r0
        dr = self.__dr
        barP0 = self.__barP0
        barM0 = self.__barM0

        init_values = {

            "r0":r0,
            "dr":dr,
            "barP0":barP0,
            "barM0":barM0,
        }
        return init_values


    #-------------------------------------------------------------

    def print_star_init_values(self):
        ''' Only to check the values. '''

        left_width = 15
        rigt_width = 15
        star_init_values = self.get_star_init_values()

        print(" ","STAR INITIAL VALUES".center(left_width+rigt_width, '-'))
        print(" ","="*(left_width+rigt_width))

        for item, value in  star_init_values.items():
            print(" ",item.ljust(left_width, '.')+" | "+'{:8.6e}'.format(value).ljust(rigt_width))
        print("\n")

#-------------------------------------------------------------
def main():
    obj = Params()
    print(obj.get_star_params())

#-------------------------------------------------------------

if __name__ == "__main__":
    main()

#-------------------------------------------------------------
