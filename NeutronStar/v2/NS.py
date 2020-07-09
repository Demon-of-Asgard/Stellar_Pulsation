#!/home/demon/anaconda3/bin/python3

#---------------------------------------------------------------------
try:
    import sys
except ImportError:
    print("sys Import error.")
try:
    import numpy as np 
except ImportError:
    print("Import error. \n Install numpy to the current venv.")
    
from BuildNS.buildNS import NS

#---------------------------------------------------------------------

def main(startype):
    Star = NS()
    Star.construct()
    # if build_stat == 0:
    #     return(0)
    # else:
    #     return(1)
    return(0)

#---------------------------------------------------------------------

if __name__ == "__main__":

    startype = " "

    try:
        startype = sys.argv[1]
    except:
        startype = "NS"
    
    retval = main(startype)

    if retval == 0:
        print("Success")
    else:
        print("Failed to execute the programme.")

#---------------------------------------------------------------------
