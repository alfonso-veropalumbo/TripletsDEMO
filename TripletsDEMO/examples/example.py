'''
Created on 13 dic 2018

@author: alfonso
'''
from compute_triplets.compute import ComputeTriplets
import matplotlib.pyplot as plt

r12 = 10.
r13 = 20.
binSize = 5.
nbar = 1.e-2

    
def go(r12, r13, binSize, nbar):

    cT = ComputeTriplets()
    points = cT.multiple_shells(r12, r13, binSize, nbar)
    print "Doing stuff..."
    triplets = cT.getTriplets(points)
    binnedTriplets = cT.getBinnedTriplets(triplets, r12=10, r13=20, binSize=5)
    print "Done!"
    cT.plotBinnedTriplets(binnedTriplets)
    plt.show()

if __name__ == '__main__':
    #profile.run('go(r12, r13, binSize, nbar)')
    go(r12, r13, binSize, nbar)
    pass