# -*- coding: utf-8 -*-
"""
Created on Thu Feb 22 09:30:04 2018

@author: rachael
"""
from __future__ import absolute_import, division, print_function
import argparse,numpy as np
import pdb
from martini22_ff import martini22
from warnings import warn

class Bead:
    """
    A class that holds all necessary information about a particular
    atom or bead
    
    ----------
    Attributes
    ----------
    resNo: int
        number of corresponding residue
    resName: string
        name of corresponding residue
    name: string
        bead/atom name
    number: int
        bead/atom number
    pos: numpy vector of floats
        bead/atom Cartesian coords
    vel: numpy vector of floats
        bead/atom velocity
    btype: string
        bead/atom type
    charge: float
        charge on the bead, by default 0
    structure: string
        secondary structure, by default random coil
    """
    def __init__(self,resNo,resname,name,number,pos,vel,btype):
        self.resNo = resNo
        self.resname = resname
        self.name = name
        self.number = number
        self.pos = pos
        self.vel = vel
        self.btype = btype
        self.charge = 0.0
        self.structure = 'C'
        
    def __str__(self):
        """
        Create a string representation as this bead should be written out
        in an itp file
        """
        s = '{:>5}{:>6}{:>6}{:>6}{:>6}{:>6}{:>8.4}; {}\n'.format(self.number,
                                                                 self.btype,
                                                                 self.resNo,
                                                                 self.resname,
                                                                 self.name,
                                                                 self.number,
                                                                 self.charge,
                                                                 self.structure)
        return s
        
class Bond:
    """
    A class that keeps track of the 2-body interactions for a topology
    (bonds and constraints)
    
    ----------
    Attributes
    ----------
    ainds: [bead,bead]
        the constituent atoms
    params: [int,num1,num2]
        the bond parameters (int is the bond type, num1 is either the bond 
        length or the table lookup number, num2 is either the bond strength
        or 1)
    notes: free form string
        any comment to go after the end, defaults to ''
    """    
    def __init__(self,ainds,params,notes = ''):
        """
        Create a bond given two atoms and some parameters
        """
        self.ainds = ainds
        self.params = params
        self.notes = notes
    
    def __str__(self):
        """
        Create a string representation as this bond should be written out
        in an itp file
        """
        s = ''
        for a in self.ainds:
            s += '\t'
            s += str(a.number)
        for p in self.params:
            s += '\t'
            s += str(p)
        s += '; '+ str(self.notes) + '\n'

        return s
        
    def contains(self,ind):
        """
        Check whether any of the atoms in the bond has a particular index
        
        ----------
        Parameters
        ----------
        ind: int
            index of atom to query
            
        -------
        Returns
        -------
        c: bool
            True if atom with index ind found, false else
        """
        for a in self.ainds:
            if a.number == ind:
                return True
        return False
    
class Angle(Bond):
    """
    A class that keeps track of the 3-body interactions (angles) for a
    topology file
    
    -------
    Parents
    -------
    Bond (same as bond except parameters and ainds are longer lists)
    
    """
    def __init__(self,ainds,params,notes=''):
        Bond.__init__(self,ainds,params,notes)
    

class Dihedral(Angle):
    """
    A class that tracks the 4-body interactions (dihedrals) for a topology file
    
    -------
    Parents
    -------
    Angle -- identical except it has long ainds and params 
    """  
    def __init__(self,ainds,params,notes=''):
        Bond.__init__(self,ainds,params,notes)
        
class Itp:
    """
    A class that holds different parts of the .itp file format ready for 
    writing

    ----------
    Attributes
    ----------
    moltype:  [string,int]
        name, exclusions
    chemName: string
        short version of chemistry name (ie DFAG)
    atomlist: list of beads containing structural information about each atom
    bondlist: list of bonds containing structural information about each bond
    conlist: list of constraints containing structural information about each 
        one
    anglist: list of angles containing structural information about each one
    dihlist: list of dihedrals containing structural information about each one
    
    """
    def __init__(self,chemName,moltype,atomlist,bondlist,conlist,anglist,
                 dihlist):
        """
        Initialize attributes
        """
        self.chemName = chemName
        self.moltype = moltype
        self.atomlist = atomlist
        self.bondlist = bondlist
        self.conlist = conlist
        self.anglist = anglist
        self.dihlist = dihlist
        
    def write(self,fname):
        """
        write out an itp file
        
        ----------
        Parameters
        ----------
        fname: string
            name of file to be written to
        """
        fid = open(fname,'w')
        fid.write('; MARTINI (martini22) Coarse Grained topology file for \
                  "Protein"\n')
        fid.write('; written by createMartiniModel for ' + self.chemName+\
                  ' chemistry\n')
        fid.write('\n')
        fid.write('[ moleculetype ]\n')
        fid.write('; Name         Exclusions\n')
        fid.write('{0}\t\t{1}\n\n'.format(self.moltype[0],self.moltype[1]))
        fid.write('[ atoms ]\n')
        for atom in self.atomlist:
            fid.write(str(atom))
        fid.write('\n')
        fid.write('[ bonds ]\n')
        for bond in self.bondlist:
            fid.write(str(bond))
        fid.write('\n')
        fid.write('[ constraints ]\n')
        for constraint in self.conlist:
            fid.write(str(constraint))
        fid.write('\n')
        fid.write('[ angles ]\n')
        for angle in self.anglist:
            fid.write(str(angle))
        fid.write('\n')
        fid.write('[ dihedrals ]\n')
        for dih in self.dihlist:
            fid.write(str(dih))
        fid.write('\n')
        fid.close()
    

class Gro:
    """
    A class that holds the different parts of the .gro file format ready for
    writing
    """
    def __init__(self,title,atomno,atoms,box):
        """
        Initializes a Gro object.
        
        ----------
        Parameters
        ----------
        title: string
            Information about the file
        atomno: int
            number of atoms in the system
        box: list of floats, length 3
            the box vectors of the system (assume cubic box)
        atoms: list of Bead objects
            containing necessary information
        """
        self.title = title
        self.atomno = atomno
        self.box  = box
        self.atoms = atoms
    
    def write(self,filename):
        fid = open(filename,'w')
        fid.write(self.title+'\n')
        fid.write(str(self.atomno)+'\n')
        bformat = "%5d%-5s%5s%5d%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f"
        for bead in self.atoms:
            beadline = bformat % (bead.resNo,bead.resname,bead.name,
                                  bead.number,bead.pos[0],bead.pos[1],
                                  bead.pos[2],bead.vel[0],bead.vel[1],
                                  bead.vel[2])
            fid.write(beadline+'\n')
        fid.write('{0} {1} {2}\n'.format(self.box[0],self.box[1],self.box[2]))
        fid.close()

class Blist:
    """
    A list containing a list of Bonds with some added useful features
    
    ----------
    Attributes
    ----------
    
    """
    def __init__(self):
        self.entries = []
    
    def __getitem__(self,key):
        return self.entries[key]
        
    def __setitem__(self,key,item):
        self.entries[key] = item
    
    def __len__(self):
        return len(self.entries)
        
    def __str__(self):
        return str([str(entry) for entry in self.entries])
        
    def append(self,entry):
        """
        Append a bond to the list
        
        ----------
        Parameters
        ----------
        entry: a Bond, Angle, or Dihedral object
        
        """
        self.entries.append(entry)
    
    def remove(self,entry):
        """
        Remove a bond from the list
        
        ----------
        Parameters
        ----------
        entry: a Bond, Angle, or Dihedral object
        
        """
        self.entries.remove(entry)
    
    def removeByIndex(self,ID):
        """
        Removes all entries in the list containing an atom with number ind
        
        ----------
        Parameters
        ----------
        ID: int
        """
        
        newEntries = []
        for entry in self.entries:
            if not entry.contains(ID):
                newEntries.append(entry)
        self.entries = newEntries        
        
        
class Topology:
    """
    contains all requisite information for a single-molecule system
    
    ----------
    Attributes
    ----------
    title: string
            Information about the file
    atomno: int
            number of atoms in the system
    box: list of floats, length 3
            the box vectors of the system (assume cubic box)
    moltype:  [string,int]
        name, exclusions
    chemName: string
        short version of chemistry name (ie DFAG)
    atomlist: list of beads containing structural information about each atom
    bondlist: list of bonds containing structural information about each bond
    conlist: list of constraints containing structural information about each 
        one
    anglist: list of angles containing structural information about each one
    dihlist: list of dihedrals containing structural information about each one
    """
    def __init__(self):
        self.title = ''
        self.atomlist = []
        self.atomno = 0
        self.box = [0.,0.,0.]
        self.moltype = ['Protein',1]
        self.chemName = 'DXXX'
        self.bondlist = Blist()
        self.conlist = Blist()
        self.anglist = Blist()
        self.dihlist = Blist()
    

class DXXXTopology(Topology):
    """
    Base topology, which inherits from DAAA as the simplest possible homodimer
    Can replace residues in the side chains
    
    ----------
    Attributes
    ----------
     title: string
            Information about the file
    atomno: int
            number of atoms in the system
    box: list of floats, length 3
            the box vectors of the system (assume cubic box)
    moltype:  [string,int]
        name, exclusions
    chemName: string
        short version of chemistry name (ie DFAG)
    atomlist: list of beads containing structural information about each atom
    bondlist: list of bonds containing structural information about each bond
    conlist: list of constraints containing structural information about each 
        one
    anglist: list of angles containing structural information about each one
    dihlist: list of dihedrals containing structural information about each one
    Itp: an Itp object
        mostly for writing
    Gro: a Gro object
        mostly for writing
    """
    def __bangle__(self,spline,nbeads):
        """
        Helper function that creates the parameters for a bond, constraint,
        or angle
        
        ----------
        Parameters
        ----------
        spline: a list of strings
        nbeads: the number of beads to go into the creation (2 for bond, 
        constraint, 3 for angle)
        
        -------
        Returns
        -------
        beads: a list of Beads
            the participating beads in the bonded interaction
        params: a list of the parameters
        notes: any comments
        """
        beads = []
        for bi in range(0,nbeads):
            beads.append(self.atomlist[int(spline[bi])-1])
        params = []
        if ';' in spline:
            cind = spline.index(';')
        else:
            cind = len(spline)
        for ci in range(nbeads,cind):
            params.append(spline[ci])
        notes = ''
        for nind in range(cind+1,len(spline)):        
            notes+= spline[nind]+' '
      
        return (beads,params,notes)
    
    def clearVars(self):
        """
        reset all variables to empty, primarily used for testing purposes
        """
        Topology.__init__(self)
        
    def __init__(self,itpname,groname):
        """
        Initialize topology from a base itp file (preferably DFAG) and grofile
        
        ----------
        Parameters
        ----------
        fname: string
            location of topology file to initialize from
        """
        Topology.__init__(self)
        fid = open(itpname)
        topology = fid.readlines()
        fid.close()
        fid = open(groname)
        grofile = fid.readlines()
        grodata = grofile[2:(len(grofile)-1)]
        box = grofile[len(grofile)-1].split()
        self.box = [float(box[0]),float(box[1]),float(box[2])]
        fid.close()
        flag = 'notes'
        for line in topology:
            if line.strip() == '[ moleculetype ]':
                flag = 'moltype'
                continue
            if line.strip() == '[ atoms ]':
                flag = 'atoms'
                continue
            if line.strip() == '[ bonds ]':
                flag = 'bonds'
                continue
            if line.strip() == '[ constraints ]':
                flag = 'cons'
                continue
            if line.strip() == '[ angles ]':
                flag = 'angles'
                continue
            if line.strip() == '[ dihedrals ]':
                flag = 'dihs'
                continue
            if flag == 'notes':
                self.title += line
            else:
                if line[0] == ';' or line == '\n':
                    continue
                else:
                    spline = line.split()
                    if flag == 'moltype':
                        self.moltype[0] = spline[0]
                        self.moltype[1] = int(spline[1])
                    elif flag == 'atoms':
                        resNo = int(spline[2])
                        resname = spline[3]
                        beadname = spline[4]
                        beadno = int(spline[0])
                        beadtype = spline[1]
                        gline = grodata[beadno-1]
                        spgline = gline.split()
                        pos = np.array([float(spgline[3]),float(spgline[4]),
                                        float(spgline[5])])
                        vel = np.array([0.,0.,0.])
                        A = Bead(resNo,resname,beadname,beadno,pos,vel,
                                 beadtype)
                        self.atomlist.append(A)
                    elif flag == 'bonds':
                        (beads,params,notes) = self.__bangle__(spline,2)
                        B = Bond(beads,params,notes)
                        self.bondlist.append(B)
                    elif flag == 'cons':
                        (beads,params,notes) = self.__bangle__(spline,2)
                        C = Bond(beads,params,notes)
                        self.conlist.append(C)
                    elif flag == 'angles':
                        (beads,params,notes) = self.__bangle__(spline,3)
                        A = Angle(beads,params,notes)
                        self.anglist.append(A)
                    elif flag == 'dihs':
                        (beads,params,notes) = self.__bangle__(spline,4)
                        D = Dihedral(beads,params,notes)
                        self.dihlist.append(D)
        self.Itp = Itp(self.chemName,self.moltype,self.atomlist,self.bondlist,
                       self.conlist,
                       self.anglist,self.dihlist)
        title = 'This file was created by createMartiniModel for the DXXX-OPV3-XXXD system with side residues PHE, ALA, and GLY'
        self.Gro = Gro(title,len(self.atomlist),self.atomlist,self.box)
    
    def write(self,fname):
        """
        Write out an itp file, a top file and a gro file corresponding to 
        the system
        
        ----------
        Parameters
        ----------
        fname: string
            the base name to use for both fname.itp and fname.gro
        """
        self.Gro.write(fname+'.gro')
        self.Itp.write(fname+'.itp') 
        top = open(fname+'.top','w')
        top.write('#include "martini.itp"\n')
        top.write('#include "{}"\n'.format(fname+'.itp'))
        top.write('[ system ]\n\n')
        top.write('; name\n')
        top.write('{} system\n\n'.format(self.chemName))
        top.write('[ molecules ]\n\n')
        top.write('; name \t number\n\n')
        top.write('{} \t {}\n'.format(self.moltype[0],self.moltype[1]))
    
    def resSwap(self,name,structure,resID):
        """
        Swap out one residue for another.  This is the meat of this whole
        structure.
        
        ----------
        Parameters
        ----------
        name: string
            name of the residue to be swapped in
        structure: string
            the 2ndary structure of the residue to be swapped in
        resID: int
            the residue number for the residue to be swapped out
        
        -----
        Notes
        -----
        (1) first index is start index
        (2) remove all atoms, bonds, constraints, angles and dihedrals 
            containing the atoms in the IDs list.  Make sure to renumber
            the atoms that come after these atoms
        (3) add all atoms, bonds, constraints, angles, and dihedrals
            necessary by the definitions found in Martini 2.2 FF.
            This includes non-explicit but generic BB angles & bonds
            (bonds should connect backbone bead to res+1 and res-1;
            angles should connect backbone bead to res+1,res+2,res-1,res-2)
            Once again make sure to renumber the atoms that come after these 
            atoms.
        """
        
        ind0 = self.removeRes(resID)
        
        self.addRes(name,structure,ind0)
    
        
    def removeRes(self,resID):
        """
        Remove the residue comprised of the atoms with the given resID
        
        ----------
        Parameters
        ----------
        IDs: list of ints
            list of indices
        
        -------
        Returns
        -------
        ind0: int
            index of the first atom ID where the residue started
        -----    
        Notes
        -----
        * find and remove all bonds, constraints, angles, dihedrals containing
        atoms with the given indices
        * find and remove all beads with the given indices
        * renumber all beads with indices larger than the given indices
        """
        IDs = []
        for atom in self.atomlist:
            if atom.resNo == resID:
                IDs.append(atom.number)
                
        for ID in IDs:
            self.bondlist.removeByIndex(ID)
            self.conlist.removeByIndex(ID)
            self.anglist.removeByIndex(ID)
            self.dihlist.removeByIndex(ID)
        
        reind = 0
        
        natomlist = []
        for atom in self.atomlist:
            if atom.number not in IDs:
                atom.number += reind
                if atom.resNo > resID:
                    atom.resNo -= 1
                natomlist.append(atom)
            else:
                reind -= 1
        self.atomlist = natomlist
        return min(IDs)
        

def main():
    """
    A function that creates and writes out itp, top, and gro files for
    use with Gromacs, for DXXX peptides. 
    Parameters are given via the command-line.
    
    ----------
    Parameters
    ----------
    residues: list of three strings
        the three amino acids of the DXXX side chain, in order
    symmetry: bool
        if symmetry is true, we are of the form DXYZ-OPV3-ZYXD
        else, we are of the form DXYZ-OPV3-XYZD
        
    -------
    Returns
    -------
    None, but prints out two things
    itp: an itp Gromacs file
    gro: a Gromacs gro file 
    """
    parser = argparse.ArgumentParser(description='get DXXX peptide structure')
    parser.add_argument('residues',metavar='R',nargs=3)
    parser.add_argument('symmetry',metavar='S',type=bool)
    parser.add_argument('pdbname',metavar='P')
    parser.add_argument('itpname',metavar='I')
    args = parser.parse_args()
    residues = args.residues
    symmetry = args.symmetry
    pdbname = args.pdbname
    itpname = args.itpname
    (Itp,Gro) = constructTopology(residues,symmetry)
    Gro.write(pdbname)
    Itp.write(itpname)

if __name__ == "__main__":
    main()