'''
Created on 13 dic 2018

@author: alfonso
'''
from compute_triplets.compute import ComputeTriplets
import matplotlib.pyplot as plt
from numpy import where, isnan

r12 = 10.
r13 = 20.
binSize = 5.
nbar = 1.e-2

def go(r12, r13, binSize, nbar):
    
    cT = ComputeTriplets()
    print "Doing stuff..."
    points = cT.multiple_shells(r12, r13, binSize, nbar)
    triplets = cT.getTriplets(points)
    #index = where(isnan(triplets["r23"]))
    #print triplets["r12"][index]
    #print triplets["r13"][index]
    #print triplets["mu"][index]
    #binnedTriplets = cT.getBinnedTriplets(triplets, r12=10, r13=20, binSize=5)
    #cT.plotBinnedTriplets(binnedTriplets)
    #plt.show()
    print "Done!"

if __name__ == '__main__':
    go(r12, r13, binSize, nbar)
    pass