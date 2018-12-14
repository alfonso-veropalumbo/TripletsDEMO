'''
Created on 13 dic 2018

@author: alfonso
'''

from generate_data.generate import GenerateData
from numpy import zeros, cos, sin, arccos, sqrt, power, pi, where, array,\
    histogram
from scipy.special import legendre
import matplotlib.pyplot as plt

class ComputeTriplets(GenerateData):
    '''
    classdocs
    '''


    def __init__(self):
        '''
        Constructor
        '''
        
    def Pl(self, mu, poly):
        return poly(mu)
    
    def limits(self, r12, r13, binSize):
        '''
        
        :param r12:
        :param r13:
        :param binSize:
        '''
        r12_min = r12
        r12_max = r12+binSize
        r13_min = r13
        r13_max = r13+binSize
    
        r23_min = abs(r12_max-r13_min)
        if ( r12_min == r13_min):
            r23_min = 0 

        r23_max = r12_max+r13_max
        nbins = int((r23_max-r23_min)/binSize)
        
        return r12_min, r12_max, r13_min, r13_max, r23_min, r23_max, nbins
    
    def getTriplets(self, points):
        '''
        
        :param points: 
        '''
        
        nShells = len(points.keys())

        bin1_out = []
        bin2_out = []
        bin3sq = []
        costheta = []
        
        #Isosceles
        for s1 in range(nShells):
            N12 = len(points[s1])
            bin1sq = power(points[s1][:,0], 2)
            for i in range(N12):
                for j in range(N12):
                    bin1_out.append(points[s1][i,0])
                    bin2_out.append(points[s1][j,0])
                    if (i==j):
                        mu=1.
                    else:
                        mu = cos(points[s1][i,1])*cos(points[s1][j,1])+\
                             sin(points[s1][i,1])*sin(points[s1][j,1])*cos(points[s1][i,2]-points[s1][j,2])
                    costheta.append(mu)
                    bin3sq.append(bin1sq[i]+bin1sq[j]-2.*mu*points[s1][i,0]*points[s1][j,0])
 
        #Non Isosceles
        for s1 in range(nShells):
            N12 = len(points[s1])
            bin1sq = power(points[s1][:,0], 2)
            for s2 in range(s1+1, nShells):
                N13 = len(points[s2])
                bin2sq = power(points[s2][:,0], 2)
                for i in range(N12):
                    for j in range(N13):
                        bin1_out.append(points[s1][i,0])
                        bin2_out.append(points[s2][j,0])
                        mu = cos(points[s1][i,1])*cos(points[s2][j,1])+\
                             sin(points[s1][i,1])*sin(points[s2][j,1])*cos(points[s1][i,2]-points[s2][j,2])
                        costheta.append(mu)
                        bin3sq.append(bin1sq[i]+bin2sq[j]-2.*mu*points[s1][i,0]*points[s2][j,0])
        
        theta = arccos(costheta)
        bin3 = sqrt(bin3sq)
        
        triangles = {}
        triangles["r12"] = array(bin1_out)
        triangles["r13"] = array(bin2_out)
        triangles["r23"] = array(bin3)
        triangles["mu"] = array(costheta)
        triangles["theta"] = array(theta)
        return triangles
    
    def getBinnedTriplets(self, triangles, r12=10., r13=50.,  binSize=5, nb = 10, norders=11):
        '''
        
        :param triangles: dictionary with triangles
        :param r12: the first shell side (lower edge)
        :param r13: the second shell side (lower edge)
        :param binSize: the bin size
        :param nb: number of bins for mu/theta binned triangles
        :param norders: number of legendre polynomials to be use
        '''

        r12_min, r12_max, r13_min, r13_max, r23_min, r23_max, nbins = self.limits(r12, r13, binSize)
    
        ww = where((triangles["r12"]>r12_min) & (triangles["r12"]<r12_max) &\
                   (triangles["r13"]>r13_min) & (triangles["r13"]<r13_max))

        triplets_mu, edges = histogram(triangles["mu"][ww], bins=nb, range=[-1, 1]) 
        bins_mu = array([0.5*(edges[i+1]+edges[i]) for i in range(nb)])
  
        triplets_theta, edges = histogram(triangles["theta"][ww], bins=nb, range=[0, pi]) 
        bins_theta = array([0.5*(edges[i+1]+edges[i]) for i in range(nb)])
  
        triplets_r23, edges = histogram(triangles["r23"][ww], bins=nbins, range=[r23_min, r23_max]) 
        bins_r23 = array([0.5*(edges[i+1]+edges[i]) for i in range(nbins)])
    
        tl = zeros(norders)
        for l in range(norders):
            tl[l] += sum(0.5*(2*l+1)*legendre(l)(triangles["mu"][ww]))
                         
        binnedTriplets = {}
        
        binnedTriplets["mu"] = array([bins_mu, triplets_mu])
        binnedTriplets["theta"] = array([bins_theta, triplets_theta])
        binnedTriplets["r23"] = array([bins_r23, triplets_r23])
        binnedTriplets["l"] = array([range(norders), tl])
        
        return binnedTriplets
    
    def plotBinnedTriplets(self, binnedTriplets):
        '''
        
        :param binnedTriplets:
        '''
        
        figure = plt.figure(figsize=(20, 5))

        ax1 = figure.add_subplot(141)
        ax1.set_xlabel(r"$\cos(\theta))$")
        ax1.set_ylabel(r"$\mathcal{T}(\cos(\theta))$")
    

        ax2 = figure.add_subplot(142)
        ax2.set_xlabel(r"$\theta$")
        ax2.set_ylabel(r"$\mathcal{T}(\theta)$")  

        ax3 = figure.add_subplot(143)
        ax3.set_xlabel(r"$r_{23}$")
        ax3.set_ylabel(r"$\mathcal{T}(r_{23})$")

        ax4 = figure.add_subplot(144)
        ax4.set_xlabel(r"$l$")
        ax4.set_ylabel(r"$\mathcal{T}_l$")
        
        if isinstance(binnedTriplets, (list,)):
            nData = len(binnedTriplets)
            for i in range(nData):
                ax1.plot(binnedTriplets[i]["mu"][0], binnedTriplets[i]["mu"][1], 'o')
                ax2.plot(binnedTriplets[i]["theta"][0], binnedTriplets[i]["theta"][1], 'o')
                ax3.plot(binnedTriplets[i]["r23"][0], binnedTriplets[i]["r23"][1], 'o')
                ax4.plot(binnedTriplets[i]["l"][0], binnedTriplets[i]["l"][1], '-')
        elif isinstance(binnedTriplets, (dict,)):
            ax1.plot(binnedTriplets["mu"][0], binnedTriplets["mu"][1], 'o')
            ax2.plot(binnedTriplets["theta"][0], binnedTriplets["theta"][1], 'o')
            ax3.plot(binnedTriplets["r23"][0], binnedTriplets["r23"][1], 'o')
            ax4.plot(binnedTriplets["l"][0], binnedTriplets["l"][1], '-')
      
        figure.tight_layout() 
        
        return figure
