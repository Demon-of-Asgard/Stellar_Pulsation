#!/home/demon/anaconda3/bin/python3

import random as rd
import numpy as np 

hitlist = [10, 100, 1000, 10000, 100000, int(1e6), int(1e7)]

for hits in hitlist:
    hitincircle = 0
    for i in range(hits):
        x = rd.uniform(0.0, 1.0)
        y = rd.uniform(0.0, 1.0)

        if x**2+y**2 <= 1.0:
            hitincircle += 1

    estimatedPI = 4.0*hitincircle/hits
    err = 100.*np.sqrt((estimatedPI - np.pi)**2)/np.pi
    print("Hits: {} Estimated pi: {} Percentage Error: {}".format(hits, estimatedPI, err))