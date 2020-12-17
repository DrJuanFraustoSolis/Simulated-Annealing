# universe.py
#
# Copyright 2007       Frank Eisenmenger, U.H.E. Hansmann,
#                      Jan H. Meinke, Sandipan Mohanty
#
import smmp
from math import *

class Universe:
    """Describes the environment of the simulation. Including the temperature,
    the solvent, and the force field. This class keeps track of all the 
    molecules.
    """
    
    def __init__(self, T = 300, ff = 'ecepp3', st = 1, libdir = 'SMMP/'):
        self.setTemperature(T)
        if ff == 'ecepp3':
            smmp.epar_l.flex = 0
            smmp.epar_l.sh2 = 0
            smmp.epar_l.epsd = 0
            smmp.epar_l.ientyp=0
        elif ff == 'ecepp2':
            smmp.epar_l.flex = 0
            smmp.epar_l.sh2 = 1
            smmp.epar_l.epsd = 0
            smmp.epar_l.ientyp=0
        elif ff == 'flex':
            smmp.epar_l.flex = 1
            smmp.epar_l.sh2 = 0
            smmp.epar_l.epsd = 0
            smmp.epar_l.ientyp=0
        smmp.epar_l.tesgrd = 0
        smmp.isolty.itysol = st
        smmp.init_energy(libdir)
        self.__objects = []
        self.__boxSize = 10
        smmp.updchois.upchswitch = 0
        smmp.updchois.rndord = 0
        smmp.updchois.bgsprob = 0
        
    def add(self, m):
        self.__objects.append(m)
        
    def energy(self):
        """Return the total energy of this universe."""
        return smmp.energy()
    
    def rgyr(self):
        """Radius of gyration of the entire system."""
        return smmp.rgyr(0)[0]
        
    def endToEnd(self):
        return smmp.rgyr(0)[1]
        
    def hbond(self):
        """Returns the total number of hydrogen bonds in the system"""
        mhb = 0
        for i in self.__objects:
            mhb = mhb + i.hbond()
        mhb = mhb + smmp.interhbond()
        return mhb
    
    def helix(self):
        return smmp.helix()
    
    def sheet(self):
        """Checks if we have any parallel or antiparallel beta-sheets in the
        system.
        """
        minStrandContent = 2
        minHydrogenBonds = 3
        res = []
        for i in range(0, len(self.__objects)):
            if self.__objects[i].strand() >= minStrandContent:
                for j in range(i, len(self.__objects)):
                    if smmp.h_bond.mmhb[self.__objects[i].id()][self.__objects[j].id()] >= minHydrogenBonds:
                        if self.__objects[j].strand() >= minStrandContent:
                            print "Found a candidate for a sheet between %s and %s." % (i, j)
                            dir1 = self.__objects[i].directionVector()
                            dir2 = self.__objects[j].directionVector()
                            sum = 0
                            for k in range(0, 3):
                                sum += dir1[k] * dir2[k]
                            if sum > 0.7:
                                res.append([i, j, 1])
                            elif sum < -0.7:
                                res.append([i, j, -1])
        return res 
    
    def temperature(self):
        """Returns the current temperature of this Universe."""
        return 1.0 / ( smmp.bet.beta * 1.98773e-3 )
    ## @brief Sets the temperature value of this Universe.
    #  Be aware that SMMP is not thread save. Changing the temperature may
    #  lead to strange behavior in the current sweep.
    #  @pre T > 0
    #
    #  @param T temperature in Kelvin.
    def setTemperature(self, T):
        if T> 0:
            smmp.bet.beta = 1.0 / ( T * 1.98773e-3 )
        
    def boxSize(self):
        return self.__boxSize

    def objects(self):
        """Returns the list of interacting things in the universe."""
        res = self.__objects
        return res
        
    def save(self):
        """Saves the state of this Universe."""
        pass
