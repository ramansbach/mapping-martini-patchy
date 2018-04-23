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
        
    ------
    Raises
    ------
    Warning, if the bond is created between two atoms whose residues are
    not numbered consecutively.  Eventually, this should possibly be ported
    to an error.
    """    
    def __init__(self,ainds,params,notes = '',fromCreateNewBonds=False):
        """
        Create a bond given two atoms and some parameters
        """
        if len(ainds) == 2:
            if abs(ainds[0].resNo - ainds[1].resNo) != 1 and (not fromCreateNewBonds):
                #pdb.set_trace()
                warn('Bond created between disjoint residues. Incorrectly ordered topology file?')
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
    
    def containsRes(self,resID):
        """
        Check whether any of the atoms in the bond belongs to a particular
        residue
        
        ----------
        Parameters
        ----------
        resID: int
            the resID of the residue we are checking for
        """
        for a in self.ainds:
            if a.resNo == resID:
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
    def __init__(self,ainds,params,notes='',fromCreateNewBonds=False):
        Bond.__init__(self,ainds,params,notes,fromCreateNewBonds)
    

class Dihedral(Angle):
    """
    A class that tracks the 4-body interactions (dihedrals) for a topology file
    
    -------
    Parents
    -------
    Angle -- identical except it has long ainds and params 
    """  
    def __init__(self,ainds,params,notes='',fromCreateNewBonds=False):
        Bond.__init__(self,ainds,params,notes,fromCreateNewBonds)
        
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
        
    def __add__(self,other):
        """
        Add two Blists together by concatenating their entries
        """
        newBlist = Blist()
        newBlist.entries = self.entries+other.entries
        return newBlist
        
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
     
        ind = 0
        for entry in self.entries:
            if not entry.contains(ID):
                newEntries.append(entry)
            
            ind += 1
        self.entries = newEntries   
    
    def insertByResIndex(self,newentries,resID):
        """
        Inserts entries between entries containing atoms with resID - 1 and 
        those with atoms with resID + 1
        
        ----------
        Parameters
        ----------
        newentries: Blist
            the new entries to be inserted
            
        resID: int
            the residue ID location to insert the entries
        """
        startInd = -1

        currind = 0
        #pdb.set_trace()
        for entry in self.entries:
            if entry.containsRes(resID - 1):
                startInd = currind
            currind += 1
        #pdb.set_trace()
        nentries = []
        for i in range(0,startInd+1):
            nentries.append(self.entries[i])
        #pdb.set_trace()
        for entry in newentries:
            nentries.append(entry)
        #pdb.set_trace()
        for i in range(startInd+1,len(self.entries)):
            nentries.append(self.entries[i])
        #pdb.set_trace()
        self.entries = nentries
        
        
        
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
    
    def findBBinRes(self,resID):
        """
        A function that locates the BB residue inside the atomlist with the
        given residue ID.  Assumes there is only one
        
        ----------
        Parameters
        ----------
        resID: int
            the residue ID to search for the backbone bead of
        
        -------
        Returns
        -------
        beadID: int
            the index of the bead in the atomlist that is the backbone bead
            for the given residue, None if there is no backbone bead
        """
        beadID = 0
        for atom in self.atomlist:
            if atom.resNo == resID and atom.name == 'BB':
                return beadID
            beadID += 1
        return None
    
    def createNewBonds(self,params,shiftID,bbAtom,resID):
        """
        Helper function that creates bonds, constraints, angles or dihedrals
        
        ----------
        Parameters
        ----------
        params: list
            These parameters are determined by the particular type of thing
            in general
        shiftID: int
            essentially what kind of bonded interaction: 1 = bond or constraint
            , 2 = angle, 3 = dihedral (it's how far out to take your window) of
            BB beads that are included
        bbAtom: Bead
            central bead
        resID: eventual ID of residue containing central bead
        -------
        Returns
        -------
        newBonds: Blist
            list of  
        """
        
        newBonds = Blist()
        finalResID = self.atomlist[len(self.atomlist)-1].resNo
        span = range(max(1,resID-shiftID),min(finalResID,
                         resID+shiftID)+1)
        #pdb.set_trace()
        for j in range(0,len(span)-shiftID):
            lrinds = span[j:j+shiftID+1]
            ainds = []
            notes = ''
            for lri in lrinds:
                if lri == resID:
                    currbead = bbAtom
                else:
                    currBeadID = self.findBBinRes(lri)
                    currbead = self.atomlist[currBeadID]
                notes += currbead.resname + '-'
                ainds.append(currbead)
                
            if shiftID == 1:
                newBonds.append(Bond(ainds,params,notes,True))
            elif shiftID == 2:
                newBonds.append(Angle(ainds,params,notes,True))
            elif shiftID == 3:
                newBonds.append(Dihedral(ainds,params,notes,True))
            else:
                warn('No such bonded interaction, returning empty Blist.')
        return newBonds
    
    def createRes(self,name,structure,bbID):
        """
        Creates a residue to be inserted based on the Martini 2.2 force field
        
        ----------
        Parameters
        ----------
        name: string
            the residue name as it appears in the Martini 2.2 ff lookup table
        structure: string
            the residue secondary structure as it appears in the Martini 2.2
            ff lookup table
        bbID: index of current BB bead where the residue is to be replaced
        -------
        Returns
        -------
        resAtoms: list of atoms
            all atoms contained in the residue
        resBonds: Blist of bonds
        resCons: Blist of constraints
        resAngs: Blist of angles
        resDihs: Blist of dihedrals
        
        -----
        Notes
        -----
        1) find type and position of BB bead
        2) find type and positions of SC beads
        3) bonds/constraints
        4) angles
        5) dihedrals
        """
        resAtoms = []
        resBonds = Blist()
        resCons = Blist()
        resAngs = Blist()
        resDihs = Blist()
        ff = martini22()
        ss = ['F','E','H','1','2','3','T','S','C']
        ssind = ss.index(structure)
        #add backbone bead
        if ff.bbtyp.has_key(name):
            bbbeadtype = ff.bbtyp[name][ssind]
        else:
            bbbeadtype = ff.bbdef[ssind]
        bbAtom = self.atomlist[bbID]
        bbbead = Bead(0,name,'BB',0,bbAtom.pos,bbAtom.vel,bbbeadtype)
        resAtoms.append(bbbead)
        
        bbl = ff.bbldef[ssind]
        bbk = ff.bbkb[ssind]
        #add backbone bonds connecting left and right
        #pdb.set_trace()
        if bbk is not None:
            newBBBonds = self.createNewBonds([1,bbl,bbk],1,bbbead,bbAtom.resNo)
            resBonds = resBonds + newBBBonds
        else:
            newBBBonds = self.createNewBonds([1,bbl],1,bbbead,bbAtom.resNo)
            resCons = resCons + newBBBonds
                
        #add backbone angles connecting left and right
        if ff.bbatyp.has_key(name):
            bba = ff.bbatyp[name][ssind]
        else:
            bba = ff.bbadef[ssind]
        bbka = ff.bbka[ssind]
        newAngs = self.createNewBonds([2,bba,bbka],2,bbbead,bbAtom.resNo)
        resAngs = resAngs + newAngs
        #if they exist, add backbone dihedrals connecting left and right
        if ssind < len(ff.bbddef):
            bbd = ff.bbbdef[ssind]
            bbkd = ff.bbkd[ssind]
            newDihs = self.createNewBonds([1,bbd,bbkd],3,bbbead,bbAtom.resNo)
            resDihs = resDihs + newDihs
           
        
        #add side chain beads         
        if len(ff.sidechains[name]) > 0:
            
            scpos = self.getSCPos(ff.sidechains[name][1],bbAtom.pos,
                                      len(ff.sidechains[name][0]))
            #pdb.set_trace()
            for sci in range(len(ff.sidechains[name][0])):
                scbead = Bead(0,name,'SC'+str(sci+1),sci+1,scpos[sci,:],
                              np.array([0.,0.,0.]),ff.sidechains[name][0][sci])
                resAtoms.append(scbead)
        
            #add bonded interactions containing SC beads
            bondConnect = ff.connectivity[name]
            bondAttributes = ff.sidechains[name]
            for bsetind in range(len(bondConnect)):
                for bind in range(len(bondConnect[bsetind])):
                    currBInds = bondConnect[bsetind][bind]
                    currBParams = bondAttributes[bsetind+1][bind]
                    currBeads = []
                    currNotes = ''
                    #pdb.set_trace()
                    for cbind in range(len(currBInds)):
                        currBeads.append(resAtoms[currBInds[cbind]])
                        currNotes += resAtoms[currBInds[cbind]].resname +'-'
                    if len(currBInds) == 2:
                        
                        if currBParams[1] is not None: #bond
                            currParams = [1,currBParams[0],currBParams[1]]
                            resBonds.append(Bond(currBeads,currParams,
                                                 currNotes))
                        else: #constraint
                            currParams = [1,currBParams[0]]
                            resCons.append(Bond(currBeads,currParams,
                                                currNotes))
                    elif len(currBInds) == 3:#angle
                        currParams = [2,currBParams[0],currBParams[1]]
                        resAngs.append(Angle(currBeads,currParams,currNotes))
                    elif len(currBInds) == 4:#dihedral
                        currParams = [2,currBParams[0],currBParams[1]]
                        resDihs.append(Dihedral(currBeads,currParams,
                                                currNotes))
                    else:#???
                        warn("Unknown bonded interaction. Not adding SC bonds.")
                        
        return (resAtoms,resBonds,resCons,resAngs,resDihs)
    
    def getSCPos(self,scs,bbPos,nBeads):
        """
        get the starting positions for the SC beads
        
        ----------
        Parameters
        ----------
        scs: Side chain bond parameters from force field
        bbPos: position of backbone bead
        nBeads: int
            number of beads in the SC
        
        -------
        Returns
        -------
        scpos: numpy array, B x 3
            the positions of the B SC beads
            
        -----
        Notes
        -----
        Right now this is pretty hacky and assumes that anything really bad
        will be sorted out during equilibration
        
        We just define everything as taking place in yz plane
        """
        scpos = np.zeros((nBeads,3))
        if nBeads < 3:
            prevPos = bbPos
            for i in range(len(scs)):
                pos = prevPos + scs[i][0]*np.array([0.,0.,1.])
                scpos[i,:] = pos
            #pdb.set_trace()
        elif nBeads == 3:
            scpos[0,:] = bbPos + scs[0][0]*np.array([0.,0.,1.])
            scpos[1,:] = scpos[0,:] + scs[1][0]*np.array([0.,np.cos(np.pi/3),
                                                          np.sin(np.pi/3)])
            scpos[2,:] = scpos[0,:] + scs[2][0]*np.array([0.,-np.cos(np.pi/3),
                                                          np.sin(np.pi/3)])
            #pdb.set_trace()
        else:
            scpos[0,:] = bbPos + scs[0][0]*np.array([0.,0.,1.])
            scpos[1,:] = scpos[0,:] + scs[1][0]*np.array([0.,np.cos(np.pi/3),
                                                          np.sin(np.pi/3)])
            scpos[2,:] = scpos[0,:] + scs[2][0]*np.array([0.,-np.cos(np.pi/3),
                                                          np.sin(np.pi/3)])
            scpos[3,:] = scpos[2,:] + scs[1][0]*np.array([0.,np.cos(np.pi/3),
                                                          np.sin(np.pi/3)])
            #pdb.set_trace()
        return scpos

class ResidueTopology(Topology):
    """
    Topology for a single amino acid residue, rather than an OPV3 core.
    Initialized as a single alanine for compatibility
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
        Topology.__init__(self)
        self.chemName = 'X'
        resNo = 1
        resname = 'ALA'
        name = 'BB'
        number = 1
        pos = np.array([0.,0.,0.])
        vel = np.array([0.,0.,0.])
        btype = 'P4'
        Ala = Bead(resNo,resname,name,number,pos,vel,btype)
        self.atomlist.append(Ala)
        self.atomno = 1
        self.title = 'Single amino acid residue created from Martini 2.2 FF'
    
    def changeAA(self,resname):
        """
        Changes the alanine into whatever other amino acid is given.
        
        ----------
        Parameters
        ----------
        resname: string
            name of the residue to swap in
        """
        (resAtoms,resBonds,resCons,resAngs,resDihs) = \
                                                 self.createRes(resname,'C',0)
        self.atomlist = resAtoms
        self.bondlist = resBonds
        self.conlist = resCons
        self.anglist = resAngs
        self.dihlist = resDihs
        for i in range(len(self.atomlist)):
            self.atomlist[i].number += 1
            self.atomlist[i].resNo = 1
        
            
        
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
        title = 'This file was created by createMartiniModel for a single res'
        gro = Gro(title,len(self.atomlist),self.atomlist,self.box)
        gro.write(fname+'.gro')
        itp = Itp(self.chemName,self.moltype,self.atomlist,self.bondlist,
                       self.conlist,
                       self.anglist,self.dihlist)
        itp.write(fname+'.itp') 
        top = open(fname+'.top','w')
        top.write('#include "martini.itp"\n')
        top.write('#include "{}"\n'.format(fname+'.itp'))
        top.write('[ system ]\n\n')
        top.write('; name\n')
        top.write('{} system\n\n'.format(self.chemName))
        top.write('[ molecules ]\n\n')
        top.write('; name \t number\n\n')
        top.write('{} \t {}\n'.format(self.moltype[0],self.moltype[1]))

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
                if line[0] == ';' or line == '\n' or line[0] == '#':
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
                        try:
                            (beads,params,notes) = self.__bangle__(spline,2)
                        except:
                            pdb.set_trace()
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
        title = 'This file was created by createMartiniModel for the DXXX-OPV3-XXXD system with side residues PHE, ALA, and GLY'
        gro = Gro(title,len(self.atomlist),self.atomlist,self.box)
        gro.write(fname+'.gro')
        itp = Itp(self.chemName,self.moltype,self.atomlist,self.bondlist,
                       self.conlist,
                       self.anglist,self.dihlist)
        itp.write(fname+'.itp') 
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
        IDs = []
        for atom in self.atomlist:
            if atom.resNo == resID:
                IDs.append(atom.number-1)
        #pdb.set_trace()
        ind0 = min(IDs)
        newRes = self.createRes(name,structure,ind0)
        #pdb.set_trace()
        self.removeRes(resID)
        #pdb.set_trace()
        
        self.addRes(newRes,resID,ind0)
    
        
    def removeRes(self,resID):
        """
        Remove the residue comprised of the atoms with the given resID
        
        ----------
        Parameters
        ----------
        IDs: list of ints
            list of indices
             
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
    
    def addRes(self,newRes,resID,ind):
        """
        Add a residue to the system at the given index
        
        ----------
        Parameters
        ----------
        newRes: structure containing (list,Blist,Blist,Blist,Blist)
            list of atoms
            Blist of bonds
            Blist of constraints
            Blist of angles
            Blist of dihedrals
        resID: int
            index of resID
        ind0: where to insert atoms
            
        -----
        Notes
        -----
            1) add atoms to atomlist, renumbering everything after them
            2) insert bonds, constraints, angles, and dihedrals
                * note to check: am I gonna fuck myself if I stick these
                  at the ends of their respective lists? Because I think
                  that might not play well with removeRes, but I'll have to
                  check
        """
        resAtoms = newRes[0]
        resBonds = newRes[1]
        resCons = newRes[2]
        resAngs = newRes[3]
        resDihs = newRes[4]
        newatoms = []        
        for i in range(0,ind):
            newatoms.append(self.atomlist[i])
        for i in range(len(resAtoms)):
            resAtoms[i].number = i+ind+1
            resAtoms[i].resNo = resID
            newatoms.append(resAtoms[i])
        for i in range(ind,len(self.atomlist)):
            atom = self.atomlist[i]
            atom.number = i+len(resAtoms)+1
            atom.resNo += 1
            newatoms.append(atom)
        #pdb.set_trace()
        self.atomlist = newatoms
        self.bondlist.insertByResIndex(resBonds,resID)
        self.conlist.insertByResIndex(resCons,resID)
        self.anglist.insertByResIndex(resAngs,resID)
        self.dihlist.insertByResIndex(resDihs,resID)
        

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
