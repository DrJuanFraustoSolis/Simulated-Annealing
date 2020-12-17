#!/usr/bin/env python
# protein.py
#
# Copyright 2007       Frank Eisenmenger, U.H.E. Hansmann,
#                      Jan H. Meinke, Sandipan Mohanty
#
import sys
import smmp
from math import *

crd = 57.2957795130823

def nursvr(ivr):
    res = 0
    for i in range(smmp.mol_i.ntlml - 1, -1, -1):
        ifirs = smmp.mol_i.irsml1[i] - 1
        if ivr + 1 > smmp.res_i.ivrrs1[ifirs] and smmp.mol_i.nvrml[i] > 0:
            for j in range(smmp.mol_i.irsml2[i] - 1, ifirs - 1, -1):
                if ivr + 1 > smmp.res_i.ivrrs1[j]  and smmp.res_i.nvrrs[j] > 0:
                    res = j + 1
                    return res
    return res

## Provides read-only access to SMMP's atoms. Usually, you won't initialize
#  an Atom directly since you need to know the internal array index of the
#  Atom in SMMP. You can get a list of all the Atoms in a Protein by using 
#  protein.Protein.getAtoms().
#  Various properties are accessible. See the functions' documentations for
#  details.
class Atom:
    """Atoms are the basic chemical units."""
    #TODO: insert namespace
    
    ## Initialize the atom with its ID. The ID is the array index of the atom
    #  in SMMP. 
    #       @param id index of atom in SMMP
    def __init__(self, id):
        self.__id = id
    
    def position(self):
        """Cartesian coordinates of this atom."""
        res = (smmp.atm_r.xat[self.__id], smmp.atm_r.yat[self.__id], smmp.atm_r.zat[self.__id])
        return res
    
    def charge(self):
        """Partial charge of this atom."""
        return smmp.atm_r.cgat[self.__id]
    
    def radius(self):
        """Van-der-Waals Radius of this atom."""
        return smmp.solvent.rvdw[self.__id]
        
    def name(self):
        return ''.join([smmp.atm_c.nmat[i][self.__id] for i in range(0,4)]).strip()
        
    def __str__(self):
        return "%s at (%s, %s, %s)" % (self.name(), self.position()[0],  self.position()[1], self.position()[2])

## Amino acids bound within a protein are often called residues. Strictly
#  speaking the residue is the part that's connected to the Calpha and not part
#  of the backbone, but we'll include the backbone atoms as well. Just like Atoms
# Residues reference the corresponding data in SMMP and are not independent. 
class Residue:
    """An amino acid within a Protein."""
    
    ## Initialize the residue with the corresponding index in SMMP.
    #  @param id index of residue in SMMP
    def __init__(self, id):
        self.__id = id

    def atoms(self):
        """Returns the list of atoms belonging to this residue."""
        first = smmp.res_i.iatrs1[self.__id]
        last = smmp.res_i.iatrs2[self.__id]
        res = [Atom(i) for i in range(first - 1, last - 1)]
        return res
    
        
## Proteins are constructed from a sequence file and a configuration file or 
#  read from a PDB file.
class Protein:
    """A protein is a polypeptide made of the 20 naturally occuring amino 
    acids.
    """
    
    def __init__(self, seq, var=' '):
        """The structure of the protein is read in from a file. 
        Pass the name of a textfile with a sequence string and--optionally--the name of a
        variable file that describes the conformation. Or give the name of pdb
        file---it must end with 'pdb'---that contains both information."""
        
        grpn = 'nh2'
        grpc = 'cooh'
        prevID = smmp.mol_i.ntlml - 1
        if seq[-3:].lower() == 'pdb':
            if smmp.mol_i.ntlml:
                print "Sorry, but I can only read the first protein from a pdb file." 
                raise ValueError
            smmp.init_molecule(0, grpn, grpc, seq)
        else:
            smmp.init_molecule(1, grpn, grpc, seq, var)
        self.__id = smmp.mol_i.ntlml - 1
        if self.__id - prevID > 1:
            print "Sorry, I don't know how to deal with multiple proteins in one sequence file, yet."
            raise NotImplementedError
    
    def setOrigin(self, x, a):
        """Sets the position of the N terminal of the protein in the global 
        coordinate system and the orientation of the protein. Expects two 
        tuples. The first contains the coordinates of the first atom (x, y,
        and z) and the second contains three angles.
        """
        #TODO: Implement type checking
        smmp.mol_r.gbpr[0][self.__id] = x[0]
        smmp.mol_r.gbpr[1][self.__id] = x[1]
        smmp.mol_r.gbpr[2][self.__id] = x[2]
        smmp.mol_r.gbpr[3][self.__id] = a[0]
        smmp.mol_r.gbpr[4][self.__id] = a[1]
        smmp.mol_r.gbpr[5][self.__id] = a[2]

    def origin(self):
        r = [smmp.mol_r.gbpr[i][self.__id] for i in range(0, 3)]
        a = [smmp.mol_r.gbpr[i][self.__id] for i in range(3, 6)]
        return r, a
    
    def getDihedral(self, angle, residue):
        try:
            res = int(residue)
            idx = res
        except ValueError:
            res = residue.split()
            res[0] = res[0].lower()
            res[1] = int(res[1])
            idx = self.findResidue(res[0], res[1])
        varIdx = ''
        for i in range(smmp.res_i.ivrrs1[idx], smmp.res_i.ivrrs1[idx] + smmp.res_i.nvrrs[idx]):
            varName = ''.join(list([str(smmp.var_c.nmvr[j][i]) for j in range(0, 3)]))
            if varName == angle:
                varIdx = i
                break
        if varIdx:
            return smmp.var_r.vlvr[varIdx] * crd
        raise KeyError, '%s not found for %s.' % (angle, residue)
    
    ## Residues are given by their 3-letter abbreviation and their index. 
    #  For example, findResidue('gly', 5) returns the index of the fifth occurence
    #  of glycine in the protein. If the residue is not found an exception is raised.
    #
    #  @param name 3-letter abbreviation of the residue
    #  @param idx index of occurence of the amino acid given in name.
    
    def findResidue(self, name, idx = 1):
        """Find the idxth residue of type name and return the index."""
        res = ''.join(self.residues()).lower()
        name = name.lower()
        j = 0
        try:
            for i in range(0, idx):
                k = res.index(name, j)
                j = j + k
        except:
            raise KeyError, "%s%s not found." % (name, idx)
        if j >=len(res):
            raise KeyError, "%s%s not found." % (name, idx)
        return j / 3
        
    def atoms(self):
        """Returns the list of atoms belonging to this protein."""
        first = smmp.res_i.iatrs1[smmp.mol_i.irsml1[self.__id] - 1]
        last = smmp.res_i.iatrs2[smmp.mol_i.irsml2[self.__id] - 1]
        res = [Atom(i) for i in range(first - 1, last - 1)]
        return res
    
    def residues(self):
        """Returns a list of the amino acids in this protein."""
        seq = [str(c[0]) for c in smmp.res_c.seq]
        first = smmp.mol_i.irsml1[self.__id]- 1
        last = smmp.mol_i.irsml2[self.__id] - 1
        return ''.join(seq).split()[first:last]
        
    def id(self):
        return self.__id
        
    def strand(self):
        """Returns the number of residues in a beta-strand-like configuration."""
        maxDev = 30.0
        phiSheet = -150.0
        psiSheet = 150.0
        strandCount = 0
        ifivr = smmp.mol_i.ivrml1[self.__id]
        i = 0
        iv = smmp.var_i.idvr[i]
        while iv < ifivr:
            i += 1
            iv = smmp.var_i.idvr[i]
        for j in range(i, ifivr + smmp.mol_i.nvrml[self.__id]):
            iv = smmp.var_i.idvr[j] - 1
            name =''. join([str(smmp.var_c.nmvr[iv * 3 + k][0]) for k in range(0,3)])
            if name == 'phi':
                phi = smmp.var_r.vlvr[iv] * crd
                psi = smmp.var_r.vlvr[smmp.var_i.idvr[j+1] - 1] * crd
                if abs(phi - phiSheet) < maxDev and abs(psi - psiSheet) < maxDev:
                    strandCount += 1
        return strandCount
    
    def directionVector(self):
        """Calculates the normalized end-to-end vector for this molecule. 
        The first and last c_alpha carbon is used as reference."""
        atoms = self.atoms()
        beginning = []
        for a in atoms:
            if a.name().strip() == 'ca':
                if beginning:
                    end = a.position()
                else:
                    beginning = a.position()
        res = [end[i] - beginning[i] for i in range(0,3)] 
        sum = sqrt(reduce(lambda x, y: x + y ** 2, res))
        res = map(lambda x: x / sum, res)
        return res
        
    def contacts(self):
        """Returns a triplet containing the total number of contacts, the number
        of native contacts, and the hemming distance. The latter two quantities 
        are only meaningful if a reference structure has been defined.
        """
        return smmp.contacts()

    def hbond(self):
        return smmp.hbond(self.__id + 1, 0)
        
    def energy(self):
        return smmp.energy()
    
    def rgyr(self):
        return smmp.rgyr(self.__id + 1) 
    
    def rmsdfun(self):
        return smmp.rmsdfun()
        
    def savePDB(self, file=""):
       """Write the current configuration to file using PDB formatting."""
       smmp.outpdb(self.__id, file)
       
    def saveVariables(self, file=""):
        """Output the current conformation in the variable format. The output 
        can be read in again and used to restart/continue a simulation.
        
        @param file if empty write to standard out
        """
        # Same functionality as smmp.outvar(self.__id + 1, 0)
        if not file:
            out = sys.stdout
        else:
            out = open(file, 'w')
        for i in range(smmp.mol_i.ivrml1[self.__id] - 1 , smmp.mol_i.ivrml1[self.__id] + smmp.mol_i.nvrml[self.__id] - 1):
            if not smmp.var_l.fxvr[i]:
                name = ''.join([str(smmp.var_c.nmvr[i * 3 + j][0]) for j in range(0,3)])
                print >>out, nursvr(i), ':', name, ':', smmp.var_r.vlvr[i] * crd
            else:
                it=smmp.var_i.ityvr[i]
                if it == 3:
                    print >>out, nursvr(i), ':', name, ':', smmp.var_r.vlvr[i] * crd, ' &'
        
    def zimmer(self):
        first = smmp.mol_i.irsml1[self.__id][0] - 1
        last = smmp.mol_i.irsml2[self.__id][0] - 1
        nrs = smmp.mol_i.irsml2[smmp.mol_i.ntlml - 1][0] - 1
        smmp.zimmer(last)
        return smmp.zimme.zimm[first:last]

    def __len__(self):
        return smmp.mol_i.irsml2[self.__id]  - smmp.mol_i.irsml1[self.__id]


if __name__ == "__main__":
    libdir = 'SMMP/'

    smmp.epar_l.flex = 0
    smmp.epar_l.sh2 = 0
    smmp.epar_l.epsd = 0
    smmp.epar_l.tesgrd = 0
    smmp.isolty.itysol = 0
    smmp.init_energy(libdir)
    
    p1 = Protein('EXAMPLES/1bdd.pdb')
    for i in p1.atoms():
        print i
    print p1.getDihedral('phi', 11)
