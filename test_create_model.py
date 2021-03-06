# -*- coding: utf-8 -*-
"""
Created on Thu Feb 22 10:08:49 2018

@author: rachael

Unit tests for createMartiniModel
"""
import numpy as np,numpy.testing as npt
from createMartiniModel import *
from martini22_ff import martini22
import pdb
def test_Bead():
    """
    test bead construction
    """
    B = Bead(1,'PHE','BB',1,np.array([0.,0.,0.]),np.array([0.,0.,0.]),'SC5')
    assert B.resNo ==1
    assert B.resname == 'PHE'
    assert B.name == 'BB'
    assert B.number == 1
    assert B.btype == 'SC5'
    npt.assert_array_equal(B.pos,np.array([0.,0.,0.]))
    npt.assert_array_equal(B.vel,np.array([0.,0.,0.]))

def test_write_gro():
    """
    test gro file construction/writing
    """
    B = Bead(1,'PHE','BB',1,np.array([0.,0.,0.]),np.array([0.,0.,0.]),'SC5')
    G = Gro('testing gro construction',1,[B],[30.3,30.3,30.3])
    G.write('test.gro')
    

def test_DXXX_bangle():
    """
    test helper function for DXXXTopology
    """
    btest = ' 1     2      8         3     1 ; GLY(C)-COA(C)\n'
    ctest = ' 2    3      1   0.20940  1 ; COP 28\n'
    atest = '   1    2    3      8         5     1 ; COP-COP-VIN\n'
    Top = DXXXTopology('homoAA.itp','homoA.gro')
    Top.clearVars()
    Top.atomlist.append(Bead(1,'TES','BB',1,np.array([0.,0.,0.]),np.array([0.,0.,0.]),'P4'))
    Top.atomlist.append(Bead(1,'TES','SC1',2,np.array([0.,0.,0.]),np.array([0.,0.,0.]),'SP4'))
    Top.atomlist.append(Bead(1,'TES','BB',3,np.array([0.,0.,0.]),np.array([0.,0.,0.]),'P4'))
    (beads,params,notes) = Top.__bangle__(btest.split(),2)
    (cons,cparams,cnotes) = Top.__bangle__(ctest.split(),2)
    (angs,aparams,anotes) = Top.__bangle__(atest.split(),3)
    assert beads[0].number == 1
    assert beads[0].resNo == 1
    assert beads[0].resname == 'TES'
    assert beads[0].name == 'BB'
    assert beads[0].btype == 'P4'
    
    assert beads[1].number == 2
    assert beads[1].resNo == 1
    assert beads[1].resname == 'TES'
    assert beads[1].name == 'SC1'
    assert beads[1].btype == 'SP4'
    
    assert beads[1] == cons[0]
    assert angs[0] == beads[0]
    assert angs[1] == cons[0]
    assert angs[2] == cons[1]
    
    assert params[0] == '8'
    assert params[1] == '3'
    assert params[2] == '1'
    
    assert cparams[0] == '1'
    assert cparams[1] == '0.20940'
    assert cparams[2] == '1'
    
    assert aparams[0] == '8'
    assert aparams[1] == '5'
    assert aparams[2] == '1'
    
    assert notes.strip() == 'GLY(C)-COA(C)'
    assert cnotes.strip() == 'COP 28'
    assert anotes.strip() == 'COP-COP-VIN'
  
    
def test_Bonds():
    """
    create some bonds of different types and test them
    """
    Top = DXXXTopology('homoAA.itp','homoA.gro')
    blist = Top.bondlist
    clist = Top.conlist
    alist = Top.anglist
    dlist = Top.dihlist
    
    bcounts = 0
    for bond in blist:
        if bond.contains(10):
            bcounts+=1
    assert bcounts == 2
    ccounts = 0
    for con in clist:
        if con.contains(9):
            ccounts +=1
    assert ccounts==2
    acounts = 0
    for ang in alist:
        if ang.contains(14):
            acounts +=1
    assert acounts == 3
    dcounts = 0
    for dih in dlist:
        if dih.contains(9):
            dcounts += 1
    assert dcounts == 4*9

def test_Blist():
    """
    test the Blist object for containing bonds
    """
    B1 = Bead(1,'TES','T',1,np.array([0.,0.,0.]),np.array([0.,0.,0.]),'T')
    B2 = Bead(1,'TES','T',2,np.array([0.,0.,0.]),np.array([0.,0.,0.]),'T')
    B3 = Bead(1,'TES','T',3,np.array([0.,0.,0.]),np.array([0.,0.,0.]),'T')
    Blst = Blist()
    B12 = Bond([B1,B2],[])
    B23 = Bond([B2,B3],[])
    B13 = Bond([B1,B3],[])
    Blst.append(B12)
    Blst.append(B23)
    Blst.append(B13)
    assert Blst[0] == B12
    Blst.remove(B12)
    assert Blst[0] == B23
    Blst.append(B12)
    Blst.removeByIndex(3)
    assert Blst[0] == B12
    
def test_DXXXTopology():
    """
    Test base initialization of DXXX Topology from file
    """
    Top = DXXXTopology('homoAA.itp','homoA.gro')
    btest = Top.atomlist[4]
    assert btest.resNo == 3
    assert btest.resname == 'COP'
    assert btest.number == 5
    assert btest.name == 'BB'
    assert btest.btype == 'SC5'
    assert btest.charge == 0.0
    assert btest.structure == 'C'
    assert Top.conlist[0].ainds[1] == btest
    assert Top.conlist[0].ainds[0].number == 3
    assert float(Top.conlist[0].params[1]) ==  0.2094
    dtest = Top.dihlist[len(Top.dihlist)-1]
    assert dtest.ainds[0] == Top.atomlist[5]
    assert dtest.ainds[1] == Top.atomlist[6]
    assert dtest.ainds[2] == Top.atomlist[8]
    assert dtest.ainds[3] == Top.atomlist[9]
    assert dtest.params[0] == '9'
    assert dtest.params[1] == '90'
    assert dtest.params[2] == '-0.006706'
    assert dtest.params[3] == '4'
    assert len(dtest.params)==4
    
def test_DXXXTopology_write():
    """
    Make sure there are no serious bugs in the write out portion
    """
    Top = DXXXTopology('DFAG.itp','DFAG.gro')
    Top.write('DFAG_auto')
    
def test_removeRes_back():
    """
    Test removing residues from a system
    """
    Top1 = DXXXTopology('homoAA.itp','homoA.gro')
    
    for i in range(2,10):
        resID = 2
        Top1.removeRes(resID)
    assert len(Top1.atomlist) == 1
    assert len(Top1.bondlist) == 0
    assert len(Top1.anglist) == 0
    assert len(Top1.conlist) == 0
    assert len(Top1.dihlist) == 0
    assert Top1.atomlist[0].number == 1
    assert Top1.atomlist[0].resNo == 1
    assert Top1.atomlist[0].resname == 'ALA'
    
def test_removeRes_front():
    """
    Test removing residues from a system
    """
    Top1 = DXXXTopology('homoAA.itp','homoA.gro')
    
    for i in range(1,9):
        resID = 1
        Top1.removeRes(resID)
    assert len(Top1.atomlist) == 1
    assert len(Top1.bondlist) == 0
    assert len(Top1.anglist) == 0
    assert len(Top1.conlist) == 0
    assert len(Top1.dihlist) == 0
    assert Top1.atomlist[0].number == 1
    assert Top1.atomlist[0].resNo == 1
    assert Top1.atomlist[0].resname == 'ALA'

    
    
def test_removeRes_midsingle():
    Top1 = DXXXTopology('homoAA.itp','homoA.gro')
    bnum = len(Top1.bondlist)
    cnum = len(Top1.conlist)
    anum = len(Top1.anglist)
    dnum = len(Top1.dihlist)
    atnum = len(Top1.atomlist)
    Top1.removeRes(4)
    assert len(Top1.atomlist) == atnum - 1
    assert len(Top1.bondlist) == bnum - 2
    assert len(Top1.conlist) == cnum
    assert len(Top1.anglist) == anum - 3
    assert len(Top1.dihlist) == dnum - 4*9
    assert Top1.atomlist[4].number == 5   
    assert Top1.atomlist[4].resname == 'COP'
    assert Top1.atomlist[5].number == 6
    assert Top1.atomlist[5].resname == 'BEN'
    assert Top1.bondlist[1].ainds[0].number == 2
    assert Top1.bondlist[1].ainds[1].number == 3
    assert Top1.bondlist[2].ainds[0].number == 8
    assert Top1.bondlist[2].ainds[1].number == 9
    
def test_removeRes_midsingle2():
    Top1 = DXXXTopology('DFAG.itp','DFAG.gro')
    bnum = len(Top1.bondlist)
    cnum = len(Top1.conlist)
    anum = len(Top1.anglist)
    dnum = len(Top1.dihlist)
    atnum = len(Top1.atomlist)
    Top1.removeRes(8)
    assert len(Top1.atomlist) == atnum - 3
    assert len(Top1.bondlist) == bnum - 2
    assert len(Top1.conlist) == cnum - 3
    assert len(Top1.anglist) == anum - 4
    assert len(Top1.dihlist) == dnum - 5*9
    assert Top1.atomlist[12].number == 13
    assert Top1.atomlist[12].resname == 'VIN'
    assert Top1.atomlist[13].number == 14
    assert Top1.atomlist[13].resname == 'VIN'
    assert Top1.bondlist[5].ainds[0].number == 12
    assert Top1.bondlist[5].ainds[1].number == 13
    assert Top1.bondlist[5].ainds[0].resname == 'COP'
    assert Top1.bondlist[5].ainds[1].resname == 'VIN'
    assert Top1.bondlist[6].ainds[0].number == 14
    assert Top1.bondlist[6].ainds[1].number == 17
    assert Top1.bondlist[6].ainds[0].resname == 'VIN'
    assert Top1.bondlist[6].ainds[1].resname == 'COP'
    
def test_removeRes_mid1():
    """
    Test removing residues from a system
    """
    Top1 = DXXXTopology('DFAG.itp','DFAG.gro')
    
    for i in range(11):
        resID = 1
        Top1.removeRes(resID)
        '''        
        print("remaining: ")
        for atom in Top1.atomlist:
            print(atom)
        print("-----------\n")
        '''
        #pdb.set_trace()
    #pdb.set_trace()
    for i in range(2):
        resID = 3
        Top1.removeRes(resID)
        #if i == 0:
            #pdb.set_trace()
    assert len(Top1.atomlist) == 2
    assert len(Top1.bondlist) == 1
    assert len(Top1.anglist) == 0
    assert len(Top1.conlist) == 0
    assert len(Top1.dihlist) == 0
    
    assert Top1.atomlist[0].number == 1
    assert Top1.atomlist[0].resNo == 1
    assert Top1.atomlist[0].resname == 'GLY'
    assert Top1.atomlist[1].number == 2
    assert Top1.atomlist[1].resNo == 2
    assert Top1.atomlist[1].resname == 'ALA'
    assert Top1.atomlist[0] == Top1.bondlist[0].ainds[0]
    assert Top1.atomlist[1] == Top1.bondlist[0].ainds[1]
    
def test_removeRes_mid2():
    """
    Test removing residues from a system
    """
    Top1 = DXXXTopology('DFAG.itp','DFAG.gro')
    
    for i in range(12):
        resID = 1
        Top1.removeRes(resID)
        '''
        print("remaining: ")
        for atom in Top1.atomlist:
            print(atom)
        print("-----------\n")
        '''
        #pdb.set_trace()
    #pdb.set_trace()
    for i in range(1):
        resID = 3
        Top1.removeRes(resID)
        #if i == 0:
            #pdb.set_trace()
    assert len(Top1.atomlist) == 5
    assert len(Top1.bondlist) == 2
    assert len(Top1.anglist) == 3
    assert len(Top1.conlist) == 3
    assert len(Top1.dihlist) == 9
    
    assert Top1.atomlist[0].number == 1
    assert Top1.atomlist[0].resNo == 1
    assert Top1.atomlist[0].resname == 'ALA'
    for i in range(1,4):
        assert Top1.atomlist[i].number == i+1
        assert Top1.atomlist[i].resNo == 2
        assert Top1.atomlist[i].resname == 'PHE'
    assert Top1.atomlist[0] == Top1.bondlist[0].ainds[0]
    assert Top1.atomlist[1] == Top1.bondlist[0].ainds[1]
    
def test_getSCPosARG():
    """
    test getting the SC position for a few key residues
    """
    ff = martini22()
    scsARG = ff.sidechains['ARG'][1]
    Top = Topology()
    scposARG = Top.getSCPos(scsARG,np.array([0.,0.,0.]),2)
    npt.assert_array_equal(scposARG,np.array([[0.,0.,0.33],[0.,0.,0.34]]))
    
def test_getSCPosPHE():
    """
    test getting SC positions for phenylalanine
    """
    ff = martini22()
    scsPHE= ff.sidechains['PHE'][1]
    Top = Topology()
    scposPHE = Top.getSCPos(scsPHE,np.array([0.,0.,0.]),3)
    npt.assert_array_almost_equal(scposPHE,np.array([[0.,0.,0.31],
                                                     [0.,0.135,0.27*np.sqrt(3)/2+0.31],
                                                     [0.,-0.135,0.27*np.sqrt(3)/2+0.31]]))

def test_getSCPosTRP():
    """
    test getting SC positions for TRP
    """                              
    ff = martini22()
    scs= ff.sidechains['TRP'][1]
    Top = Topology()
    scpos = Top.getSCPos(scs,np.array([0.,0.,0.]),4)
    npt.assert_array_almost_equal(scpos,np.array([[0.,0.,0.3],
                                                     [0.,0.135,0.27*np.sqrt(3)/2+0.3],
                                                     [0.,-0.135,0.27*np.sqrt(3)/2+0.3],
                                                     [0.,0.,0.27*np.sqrt(3)+0.3]]))       

def check_file_equivalency(fname1,fname2):
    """helper function for making sure files are the same"""
    f1 = open(fname1)
    ft1 = open(fname2)
    f1r = f1.readlines()
    f1.close()
    ft1r = ft1.readlines()
    ft1.close()
    for lind in range(len(f1r)):
        assert(f1r[lind]==ft1r[lind])
                                                         
def test_change_PHE_to_ALA():
    """
    Test swapping residues
    """
    Top1 = DXXXTopology('DFAG.itp','DFAG.gro')
    Top1.resSwap('ALA','C',2)
    Top1.write('DAAG_half_test')
    #check that this is the same as one that we have visually confirmed looks ok
    for suffix in ['.gro','.itp','.top']:
        check_file_equivalency('DAAG_half_test'+suffix,
                               'DAAG_half_test_checked'+suffix)
    
def test_change_DFAG_to_DAAG():
    Top1 = DXXXTopology('DFAG.itp','DFAG.gro')
    Top1.resSwap('ALA','C',2)
    Top1.resSwap('ALA','C',14)
    Top1.write('DAAG_test')
    for suffix in ['.gro','.itp','.top']:
        check_file_equivalency('DAAG_test'+suffix,
                               'DAAG_test_checked'+suffix)
    
def test_change_ALA_to_PHE():
    Top1 = DXXXTopology('DFAG.itp','DFAG.gro')
    Top1.resSwap('PHE','C',3)
    Top1.write('DFFG_half_test')
    for suffix in ['.gro','.itp','.top']:
        check_file_equivalency('DFFG_half_test'+suffix,
                               'DFFG_half_test_checked'+suffix)

def test_change_DFAG_to_DFFG():
    Top1 = DXXXTopology('DFAG.itp','DFAG.gro')
    Top1.resSwap('PHE','C',3)
    Top1.resSwap('PHE','C',13)
    Top1.write('DFFG_test')
    for suffix in ['.gro','.itp','.top']:
        check_file_equivalency('DFFG_test'+suffix,
                               'DFFG_test_checked'+suffix)
                               
def test_change_PHE_to_TRP():
    Top1 = DXXXTopology('DFAG.itp','DFAG.gro')
    Top1.resSwap('TRP','C',2)
    Top1.write('DTAG_half_test')
    
def test_change_DFAG_to_DTAG():
    Top1 = DXXXTopology('DFAG.itp','DFAG.gro')
    Top1.resSwap('TRP','C',2)
    Top1.resSwap('TRP','C',14)
    Top1.write('DTAG_test')
    
def test_change_ALA_to_VAL():
    Top1 = DXXXTopology('DFAG.itp','DFAG.gro')
    Top1.resSwap('VAL','C',3)
    Top1.write('DFVG_half_test')    
    
def test_change_DFAG_to_DFVG():
    Top1 = DXXXTopology('DFAG.itp','DFAG.gro')
    Top1.resSwap('VAL','C',3)
    Top1.resSwap('VAL','C',13)
    Top1.write('DFVG_test')
    
def test_change_DFAG_to_DTVG():
    Top1 = DXXXTopology('DFAG.itp','DFAG.gro')
    Top1.resSwap('TRP','C',2)
    Top1.resSwap('TRP','C',14)
    Top1.write('DTVG_midway_test')
    Top1.resSwap('VAL','C',3)
    Top1.write('DTVG_2midway_test')
    Top1.resSwap('VAL','C',13)
    Top1.write('DTVG_test')
