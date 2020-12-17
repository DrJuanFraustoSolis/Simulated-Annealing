#!/usr/bin/env python
# main.py
#
# Copyright 2007       Frank Eisenmenger, U.H.E. Hansmann,
#                      Jan H. Meinke, Sandipan Mohanty

## \mainpage pySMMP: Python bindings to SMMP v. 1.1.
#    SMMP is a FORTRAN library for simulating proteins. It provides 
#    minimization routines and various Monte Carlo methods. SMMP uses an all-atom
#    representation of the protein. It stores the configuration within the 
#    standard geometry model. In the standard geometry model the internal degrees
#    of freedom, e.g., the dihedral angles, the bond lengths, etc. are stored and
#    the cartesian coordinates are calculated on the fly. For more details on 
#    the model SMMP uses and the functions it provides see the manual (manual.ps).
#   
#    pySMMP builds on SMMP. It provides python bindings to the functions and 
#    global variables, but it uses classes to represent the Proteins and 
#    algorithms. 
#    
#    \section sec_tutorial A Tutorial: A canonical Monte Carlo simulation of protein A
#
#    The basic setup is the same for most simulations. This tutorial takes you
#    step by step to a running simulation. First we need to import the library
#    \code 
#         import smmp
#     \endcode
#    This makes the smmp routines available to the Python interpreter. Next we
#    create the universe by initializing a universe.Universe object by first 
#    importing the universe package and then setting the object:
#    \code
#        import universe
#        myUniverse = universe.Universe(T=300)  
#    \endcode
#    This initializes the smmp library as well. We also set the initial 
#    temperature of the universe to 300K. The Universe object initializes the 
#    SMMP library and keeps track of global parameters, e.g, the temperature T.
#    After you created a Universe, you can initialize a Protein. Proteins can
#    be read from a sequence and a variable file or from a PDB file.
#    \code
#        import protein
#        
#        protA = protein.Protein("1bdd.seq", "1bdd.var")
#        myUniverse.add(protA)
#    \endcode
#    Adding the protein to myUniverse makes the universe aware of it.
#    Now you are ready to run a simulation. You can also querry the protein for 
#    various properties, e.g., 
#    \code
#        print protA.rgyr()
#    \endcode
#    prints the radius of gyration of the protein.
#    A canonical Monte Carlo algorithm is available in the algorithms package.
#    \code
#        import algorithms
#        myMC = algorithms.CanonicalMonteCarlo(myUniverse, 1000, 10000)
#    \endcode
#    This sets up a simulation with 1000 steps for equilibration and a default 
#    length of 10000 sweeps.
#    The quickest way to run the simulation is by calling the run method of the 
#    algorithm object myMC to do the entire Monte Carlo run.
#    \code
#        myMC.run()
#    \endcode

import sys, time
from optparse import OptionParser
from protein import *
from universe import *
from algorithms import CanonicalMonteCarlo

if __name__ == "__main__":
# Set up command line options    
    usage = "usage: %prog [options] seq-filename"
    parser = OptionParser(usage = usage, version="%prog 20050805")
    parser.add_option("-v", "--variables", dest="varFileName",
                  help="read amino-acid conformation from FILE", metavar="FILE")
    parser.add_option("--equi", type="int", dest="nequi", default=100,
                  help="number of sweeps for equilibration", metavar="n")
    parser.add_option("--sweeps", type="int", dest="sweeps", default=10000,
                  help="number of sweeps for annealing", metavar="n")
    parser.add_option("--mes", type="int", dest="mes", default=10,
                  help="number of sweeps between measurements", metavar="n")
    parser.add_option("--Tmax", type="float", dest="tmax", default=1000,
                  help="maximum temperature for annealing", metavar="T")
    parser.add_option("--Tmin",  type="float", dest="tmin", default=100, 
                  help="maximum temperature for annealing",metavar="T")

    (options, args) = parser.parse_args()
    if len(args) < 1:
        parser.print_help()
        sys.exit()
    seqFileName = args[0]
    
# Start simulation setup.    
    myUniv = Universe(T=300, st=0)   
    p1 = Protein(seqFileName, options.varFileName)
    p1.setOrigin((0, 0, 0), (0, 0, 0))
    myUniv.add(p1)
    T = 300.0
    smmp.outpdb(1, "start.pdb")
    print "E of starting configuration.", myUniv.energy()
    ts = time.time()
    smmp.minim(1, options.sweeps, 1.0E-12)

    smmp.outpdb(1, "equi.pdb")
    
    t1 = time.clock()
    print myUniv.energy(), myUniv.rgyr(), myUniv.helix(), p1.hbond(), (t1-ts)
    p1.saveVariables("final.var")
    smmp.outpdb(0, "final.pdb")
    te = time.time()
    print "The calculation took %ss." % (te - ts)
 
    
