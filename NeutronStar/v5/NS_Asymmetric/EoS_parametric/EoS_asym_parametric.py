from StarParams.params import Parameters as Prm 

class EoSAsymParametric():

    def __init__(self):
        self.__A = 0.8642

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
        
    def get_EoS_params(self):
        EoS_params = {
            "A":self.__A,
            "e0":5.487718e+36,
        }
        return EoS_params
     #-----------------------------------------------------------------------------------

    def barE(self, barP):
        barE = self.__A*barP**(1.0/2.0)
        return barE
        
#-----------------------------------------------------------------------------------
