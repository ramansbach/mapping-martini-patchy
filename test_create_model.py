# -*- coding: utf-8 -*-
"""
Created on Thu Feb 22 10:08:49 2018

@author: rachael

Unit tests for createMartiniModel
"""
import numpy as np,numpy.testing as npt
from createMartiniModel import *

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