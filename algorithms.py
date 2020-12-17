# algorithms.py
#
# Copyright 2007       Frank Eisenmenger, U.H.E. Hansmann,
#                      Jan H. Meinke, Sandipan Mohanty
#
import copy
import random
import smmp
from math import *
from universe import *
from protein import *

class CanonicalMonteCarlo:
    """An implementation of a canonical Monte Carlo run using Metropolis weights.
    After a Monte-Carlo sweep over all internal variables, this implementation 
    calculates the interactions between all Protein s in the Universe
    """
    
    def __init__(self, myUniverse, nequi, sweeps, weight =  lambda x : x * smmp.bet.beta):
        self.__defSweeps = sweeps
        self.__defEqui = nequi
        self.__univ = myUniverse
        self.__acc = 0
        self.__performedSweeps = 0
        self.__minEnergy = myUniverse.energy()
        self.__l = 0.1
        self.__minCounter = 0
        self.__weight = weight
        
    ## Fraction of accepted Monte Carlo moves of the total number of attempted moves.
    #  @pre At least one move must have been performed before this function can be
    #       called
    #  @return fraction of accepted moves as a number between 0 and 1.
    def acceptanceRate(self):
        return float(self.__acc) / (self.__performedSweeps * smmp.mol_par.nvr)
    ## The lowest energy ever seen at the end of a Monte Carlo move.
    #  @return lowest energy seen at the end of a sweep, sweep
    def minimumEnergy(self):
        return self.__minEnergy, self.__minCounter
    ## Number of sweeps performed so far
    def performedSweeps(self):
        return self.__performedSweeps
        
    def run(self):
        """A complete Monte Carlo run including initial equilibration."""
        self.sweeps(self.__defSweeps)
        
    def equilibrate(self):
        """Do nequi sweeps to eqilibrate the system."""
        self.sweeps(self.__defEqui)
        
    def sweeps(self, n):
        """Perform n sweeps."""
        for i in range(0, n):
            self.sweep()
        
    def sweep(self):
        """Performs a single Monte Carlo sweep."""
        eol = self.__univ.energy()
        eol, self.__acc = smmp.metropolis(eol, self.__acc, self.__weight)
        #eol = self.metropolis()
        self.__performedSweeps += 1
        if eol < self.__minEnergy:
            self.__minEnergy = eol
            self.__minCounter = self.__performedSweeps
            smmp.outpdb(0, "best.pdb")

    def metropolis(self):
        """Implementation of the Metropolis algorithm."""
        eol = self.__univ.energy()
        # Update the internal degrees of freedom first
        for i in range(0, smmp.mol_par.nvr):
            jv = smmp.var_i.idvr[i] - 1
            vrol = smmp.var_r.vlvr[jv]
            dv=smmp.var_r.axvr[jv] * (random.random() - 0.5)
            smmp.var_r.vlvr[jv] = (vrol + dv + pi) % (2 * pi) - pi
            etrial = self.__univ.energy()
            dW = self.__weight(etrial) - self.__weight(eol)
            if dW < 0:
                eol = etrial
                self.__acc += 1
            else:
                rd =random.random()
                dW = min(dW, 100)
                if  exp(-dW) > rd:
                    eol = etrial
                    self.__acc +=1
                else:
                    smmp.var_r.vlvr[jv] = vrol
        # Update the global degrees of freedom
        return eol
