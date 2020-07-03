try:
    import numpy as np 
except ImportError:
    print("Numpy-ImportError. Please install numpy to the current venv. Thank you!")

from EoS.eos import eqn_state # Reminder: The class eqn_state inherits the class params

class NS(eqn_state):

    def __init__(self):
        super().__init__()
    
    def Pfunc(self):
        print("In Pfunc()")
        self.eos()

    def Mfunc(self):
        pass
    
    def buildNS(self):
        params = self.get_params()
        for index in params:
            print(" |{} : {}".format(index,params[index]))
        
        self.Pfunc()
            
