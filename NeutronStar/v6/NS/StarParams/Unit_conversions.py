''' Unit conversions. '''

import sys

class UnitConverter:

    def __init__(self):

        self.pi = 3.1415926535 # Dimension less
        self.h = 6.62607004e-27  # [erg s]
        self.c = 3.0e10 # cm/s
        self.hbar = self.h/(2.0*self.pi)  # [erg s]
        self.MeV_by_Erg =  1.602176487E-6 # Dimesion less

    def get_conv_factors(self):

        INV_MEV_TO_CM = self.hbar*self.c/self.MeV_by_Erg
        INV_MEV_TO_SEC = self.hbar/self.MeV_by_Erg
        MEV_TO_ERG = self.MeV_by_Erg

        dict_conv_factors = {
            "INV_MEV_TO_CM":INV_MEV_TO_CM,
            "INV_MEV_TO_SEC":INV_MEV_TO_SEC,
            "MEV_TO_ERG":MEV_TO_ERG,
            "CM_TO_INV_MEV":1.0/INV_MEV_TO_CM,
            "SEC_TO_INV_MEV": 1.0/INV_MEV_TO_SEC,
            "ERG_TO_MEV":1.0/MEV_TO_ERG,
        }
        return dict_conv_factors

    def convert(self, value,from_unit, to_unit):

        if from_unit == "MeV^-1" and to_unit == "cm":
            new_value = value*self.INV_MEV_TO_CM

        if from_unit == "MeV^-1" and to_unit == "s":
            new_value = value*self.INV_MEV_TO_SEC

        if from_unit == "MeV" and to_unit == "erg":
            new_value = value*self.MEV_TO_ERG

        return new_value
       
def main(value, from_unit, to_unit):
    
    boundary  = "="*50

    print("",boundary, "\n")

    conv_obj = UnitConverter()
    dict_conv_factors = conv_obj.get_conv_factors()

    for index in dict_conv_factors:
        print(" | {} : {:8.6e}".format(index, dict_conv_factors[index]))
    
    new_value = conv_obj.convert(value, from_unit, to_unit)

    print("\n | {:6.4e} [{}] is {:6.4e} [{}]".format(value, from_unit, new_value, to_unit))

    print("",boundary, "\n")



if __name__ == "__main__":

    try:
        value = float(sys.argv[1])
    except:
        print("Missing the value. ")
        exit(1)
    
    try:
        from_unit = sys.argv[2]
    except:
        print("Missing the 'from-unit'.")
        exit(2)
    
    try:
        to_unit = sys.argv[3]
    except:
        print("Missing the 'to-unit'.")
        exit(3)

    main(value, from_unit, to_unit)