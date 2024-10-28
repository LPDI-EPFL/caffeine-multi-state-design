#!/usr/bin/python
from dynamicAS.DynamicAS import *
from bp.BPsolver import *
from config.config import * 
from ematrix.EMatrix import * 
from util.createMatrix import *
import time
import sys
#if len(sys.argv) != 2:
#  print "Usage: "
#  print sys.argv[0]+" {energyMatrixFile}"
#  sys.exit()

matrix_directory_bcl_bh3_groupB = "/Users/pablo/cur/msd_data/02b-energyMatrices_groupB/bclxl_bh3/3fdl_bclxl_bh3_min_0001"
matrix_directory_bcl_bh3_groupA = "/Users/pablo/cur/msd_data/02a-energyMatrices_groupA/bclxl_bh3/3fdl_bclxl_bh3_min_0001"

bcl_bh3_emat = EMatrix(matrix_directory_bcl_bh3_groupA)
# Allowed AAs: 
print "Allowed AAs: "
res_id_tups = bcl_bh3_emat.getSortedResIds()
for res_ix in range(len(res_id_tups)):
  print res_id_tups[res_ix]
  print bcl_bh3_emat.singleBody[res_id_tups[res_ix]]['_order']
(numResidues,rotsPerRes, eIntraMat, ePairMat, map_rots_to_seq) = bcl_bh3_emat.getReducedRotamerOnlyMatrices()
#dynas = DynamicAS(numResidues,rotsPerRes, eIntraMat, ePairMat)
#node = dynas.doAS()
#gmec = node.computeConfSoFar()
#map_resix_to_protix = bcl_bh3_emat.getSortedResIds()
#for res_ix in range(len(gmec)):
#  print `map_resix_to_protix[res_ix]`+" "+map_rots_to_seq[res_ix][gmec[res_ix]]

## Compute the BP energy
bpsolver = BPsolver(numResidues, rotsPerRes, eIntraMat, ePairMat)
bpsolver.optimize()
energy, lmec = bpsolver.computeEnergy()
map_resix_to_protix = bcl_bh3_emat.getSortedResIds()
print "BP energy = "+`energy`
for res_ix in range(len(lmec)):
  print `map_resix_to_protix[res_ix]`+" "+map_rots_to_seq[res_ix][lmec[res_ix]]

