'''
Created on 12 dic 2018

@author: alfonso
'''

import numpy as np

class GenerateData(object):
    '''
    This class generate data in two or
    multiple spherical shells.
    '''


    def __init__(self):
        '''
        Constructor
        '''
    def shell_volume(self, rMin, rMax):
        '''
        Computes the volume of a spherical shell
        :param rMin: minimum scale [Mpc/h]
        :param rMax: maximum radius [Mpc/h]
        '''
        return 4.*np.pi/3*(np.power(rMax,3)-np.power(rMin, 3))

    def two_shells(self, r12_min, r12_max, r13_min, r13_max, nbar):   
        '''
        Generate points in two spherical shells centered at the origin.   
        Return a dictionary whose keys and values are respectively
        the bin values and an array with spherical coordinates of 
        the points in each bin
        
        :param r12_min: first shell minimum radius [Mpc/h]
        :param r12_max: first shell maximum radius [Mpc/h]
        :param r13_min: second shell minimum radius [Mpc/h]
        :param r13_max: second shell maximum radius [Mpc/h]
        :param nbar: density[(Mpc/h)^-3]
        '''
        V12 = self.shell_volume(r12_min, r12_max)
        N12 = int(V12*nbar)

        V13 = self.shell_volume(r13_min, r13_max)
        N13 = int(V13*nbar)

        #Objects in the first shell
        bin1 = np.random.random(N12)*(r12_max-r12_min)+r12_min
        theta1 = np.random.random(N12)*np.pi
        phi1 = np.random.random(N12)*2*np.pi
        
        #Objects in the second shell
        bin2 = np.random.random(N13)*(r13_max-r13_min)+r13_min
        theta2 = np.random.random(N13)*np.pi
        phi2 = np.random.random(N13)*2*np.pi 
        
        points = {}
        points[0] = np.array([bin1, theta1, phi1]).T
        points[1] = np.array([bin2, theta2, phi2]).T
        
        return points

    def multiple_shells(self, rMin, rMax, binSize, nbar):   
        '''
        Generate points in spherical shells centered at the origin   
        from a minimum to a maximumx radius with an user-defined 
        bin size.
        Return a dictionary whose keys and values are respectively
        the bin values and an array with spherical coordinates of 
        the points in each bin

        :param rMin: minimum scale [Mpc/h]
        :param rMax: maximum scale [Mpc/h]
        :param binSize: bin size [Mpc/h]
        :param nbar: density[(Mpc/h)^-3]
        '''
        
        r1 = rMin
        r2 = rMin+binSize
        bb = 0
        points = {}
        while (r1<rMax+binSize):
            VV = self.shell_volume(r1, r2)
            NN = int(VV*nbar)

            rr = np.random.random(NN)*(r2-r1)+r1
            theta = np.random.random(NN)*np.pi
            phi = np.random.random(NN)*2*np.pi
            
            points[bb] = np.array([rr, theta, phi]).T
            bb += 1
            r1 = r2
            r2 += binSize
    
        return points