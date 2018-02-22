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
            s += str(a)
        for p in self.params:
            s += '\t'
            s += str(p)
        s += '; '+ str(self.notes) + '\n'

        return s
    
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
        Bond.__init__(ainds,params,notes)
    

class Dihedral(Angle):
    """
    A class that tracks the 4-body interactions (dihedrals) for a topology file
    
    -------
    Parents
    -------
    Angle -- identical except it has 
    """  
    def __init__(self,ainds,params,notes=''):
        Bond.__init__(ainds,params,notes)
        
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
        fid.write('{0}\t\t{1}\n'.format(self.moltype[0],self.moltype[1]))
        fid.write('[ atoms ]\n')
        for atom in self.atomlist:
            fid.write(str(atom))
        fid.write('[ bonds ]\n')
        for bond in self.bondlist:
            fid.write(str(bond))
        
        fid.write('[ constraints ]\n')
        for constraint in self.conlist:
            fid.write(str(constraint))
            
        fid.write('[ angles ]\n')
        for angle in self.anglist:
            fid.write(str(angle))
            
        fid.write('[ dihedrals ]\n')
        for dih in self.dihlist:
            fid.write(str(dih))
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