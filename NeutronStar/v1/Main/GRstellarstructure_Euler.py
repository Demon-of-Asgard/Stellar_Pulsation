#!/home/demon/anaconda3/bin/python3

''' This program construct the mass and pressure distribution within a neutron star using TOV equations.
    Refs:
        1)  Neutronstar for undergraduate.
            Authors: Richard R. Silbar, Sanjay Reddy.
            Arxiv: 0309041
            Url: https://arxiv.org/abs/nucl-th/0309041
        
        2)  Compact stars for undergraduates.
            Authors: I. Segert, M. Hampel, C. Greiner, J. S. Beilich.
            Arxiv: 0506417v1
            Url: https://arxiv.org/abs/astro-ph/0506417
'''

#-------------------------------------------------------------------------------

import sys
import numpy as np
import matplotlib.pyplot as plt
from starparams.star_params import params as prm 
from BuildStar.build_star import build_star
from Plot.plot import plot

#-------------------------------------------------------------------------------

def getparams(type):
    ''' Get the predefined parameters depending on the type (relativistic or non-relativistic) of the
        star. '''

    sp = prm() #Instantiating the params class object
    if type == "rl":
        params = sp.get_rlparams()
    if type == "nrl":
        params = sp.get_nrlparams()
    return params

#-------------------------------------------------------------------------------

def main(type, outfname):
    params = getparams(type)
    outputfile = build_star(params, outfname)
    plot(outputfile)

#-------------------------------------------------------------------------------

if __name__ == "__main__":
    '''Handles the sys.argv[] and make call to the function main().'''

    try:
        type = sys.argv[1]
    except:
        type = "rl"
        
    if type == "nrl":
        print("\n <::Cooking a non-relativistic star::>")
    else:
        print("\n <::Cooking a relativistic star::>")

    try:
        outfname = sys.argv[2]
        print("\n Set output file name to '{}'".format(outfname))
    except:
        outfname = "./Output/output.dat"
        print("\n Set output file name to '{}'".format(outfname))
    
    

    main(type, outfname)

#-------------------------------------------------------------------------------
