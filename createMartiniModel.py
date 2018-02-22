# -*- coding: utf-8 -*-
"""
Created on Thu Feb 22 09:30:04 2018

@author: rachael
"""
import argparse

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
    """
    def __init__(self,resNo,resname,name,number,pos,vel,btype):
        self.resNo = resNo
        self.resname = resname
        self.name = name
        self.number = number
        self.pos = pos
        self.vel = vel
        self.btype = btype
    

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

class Itp:
    """
    A class that holds the different parts of an itp file for construction
    and writing
    
    ----------
    Attributes
    ----------
    
    """        


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
    Pdb.write(pdbname)
    Itp.write(itpname)

if __name__ == "__main__":
    main()