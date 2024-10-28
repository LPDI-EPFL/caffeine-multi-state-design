#!/usr/bin/python
# This class computes a lower bound on the energy of a protein design problem, based on 
# Edge MPLP (NMPLP)  by globerson and Jaakkola
# Initialize messages to any value
# Iterate until convergence
from math import *
import random
import ipdb
from util.createMatrix import *
from util.mymath import *
from config import *
import time
class MPLPorig:
  def __init__(self,numRes, rotsPerRes, intraEmat, pairEmat):
    self.intraEmat = intraEmat
    self.pairEmat = pairEmat
    self.numRes = numRes
    self.rotsPerRes = rotsPerRes
    self.interactionGraph = create2DResMatrix(numRes, True)
    self.interactionGraph[0][4] = False
    self.interactionGraph[4][0] = False
    self.interactionGraph[0][0] = False
    self.interactionGraph[1][1] = False
    self.interactionGraph[2][2] = False
    self.interactionGraph[3][3] = False
    self.interactionGraph[4][4] = False
    self.unifiedEmat = self.mergeIntraAndPairMats()

  def getNeighbors(self, resI):
    neighbors = []
    for resJ in range(self.numRes):
      if self.interactionGraph[resI][resJ]:
        neighbors.append(resJ)
    print "Neighbors "+`resI`+" = "+`neighbors`
    return neighbors



  # Computes the low-energy bound using the EMPLP algorithm
  # availableRots is a "numRes*rots" matrix that tells us which 
  # rotamers at each position are available.  Only those rotamers
  # will be considered for the calculation.
  # Also checks if the solution is integer
  def optimizeEMPLP(self, availableRots,iterations):
    # Order of updates.  For now, we just use the 0...numRes order
    nodeOrder = range(self.numRes)
#    random.shuffle(nodeOrder)
    # l is the message matrix for EMPLP; we set all messages to zero
    self.l = create3DMsgMat(self.numRes,self.rotsPerRes,0.0)

    # The belief on each rotamer
    self.belief = create2DRotMatrix(self.numRes, self.rotsPerRes, 0.0)
    
    # EMPLP algorithm.  
    # The algorithm should loop until convergence.  For now we do it I times. 
    # Complexity of EMPLP: O(I*numRes*numRes*rotsPerRes*rotsPerRes*numRes) 
    for iteration in range(iterations):
      for resI_index in range(self.numRes):
        ipdb.set_trace()
        resI = nodeOrder[resI_index]
        for resJ in self.getNeighbors(resI): 
            # We first update self.l[resJ][resI][rotIR] and immediately after we update self.l[resI][resJ][rotJS] 
            for rotIR in availableRots[resI]:
#            for rotIR in range(self.rotsPerRes[resI]):
              self.belief[resI][rotIR] -= self.l[resJ][resI][rotIR]
              incMesgsToRotIR_fromAllExcept_resJ = self.belief[resI][rotIR] # - self.l[resJ][resI][rotIR] We already substracted the message from resJ.
              msgsFromRotsAtJ_to_rotIR = []
              for rotJS in availableRots[resJ]:
#                incMesgsToRotJS_fromAllExcept_resI = self.belief[resJ][rotJS] - self.l[resI][resJ][rotJS]
                msgsFromRotsAtJ_to_rotIR.append(self.belief[resJ][rotJS] - self.l[resI][resJ][rotJS] + self.unifiedEmat[resI][rotIR][resJ][rotJS])

              self.l[resJ][resI][rotIR] = -0.5*incMesgsToRotIR_fromAllExcept_resJ + 0.5*min(msgsFromRotsAtJ_to_rotIR)
              self.belief[resI][rotIR] += self.l[resJ][resI][rotIR]

            # Now we update self.l[resI][resJ][rotJS] 
            for rotJS in availableRots[resJ]:
#            for rotJS in range(self.rotsPerRes[resJ]):
              self.belief[resJ][rotJS] -= self.l[resI][resJ][rotJS]
              incMesgsToRotJS_fromAllExcept_resI = self.belief[resJ][rotJS]
              msgsFromRotsAtI_to_rotJS = []
              for rotIR in availableRots[resI]:
#                incMesgsToRotIR_fromAllExcept_resJ = self.belief[resI][rotIR] - self.l[resJ][resI][rotIR]
                msgsFromRotsAtI_to_rotJS.append(self.belief[resI][rotIR] - self.l[resJ][resI][rotIR] + self.unifiedEmat[resI][rotIR][resJ][rotJS])

              self.l[resI][resJ][rotJS] = (-0.5)*incMesgsToRotJS_fromAllExcept_resI + 0.5*min(msgsFromRotsAtI_to_rotJS)
              self.belief[resJ][rotJS]+= self.l[resI][resJ][rotJS]
    Ebound = 0.0
    # To check whether the solution is integer, we must check whether there is only one minimum to each belief.
    sol_is_integer = True
    for resI in nodeOrder:
      # This can be done more efficiently but I'm lazy.
      minForResI = min(self.belief[resI])
      print "min for res"+`resI`+" = "+`minForResI`
      countEqual = 0
      for rotIRbelief in self.belief[resI]:
        if eq(rotIRbelief, minForResI):
          countEqual += 1
      if countEqual >1:
        sol_is_integer = False
#        pdb.set_trace()
      Ebound += minForResI
   
    return Ebound, sol_is_integer

  # From the intra and pair matrices we create a single one by dividing the intra energy
  #   of a vertex by (n-1) and adding it to the edge.  
  def mergeIntraAndPairMats(self):
    unifiedEmat = create4DRotMatrix(self.numRes, self.rotsPerRes, 0.0)
    for resI in range(self.numRes):
      for rotIR in range(self.rotsPerRes[resI]): 
        for resJ in self.getNeighbors(resI):
          for rotJS in range(self.rotsPerRes[resJ]): 
            unifiedEmat[resI][rotIR][resJ][rotJS] = self.pairEmat[resI][rotIR][resJ][rotJS]
            unifiedEmat[resI][rotIR][resJ][rotJS] += self.intraEmat[resI][rotIR]/float(len(self.getNeighbors(resI)))
            unifiedEmat[resI][rotIR][resJ][rotJS] += self.intraEmat[resJ][rotJS]/float(len(self.getNeighbors(resJ)))
            unifiedEmat[resJ][rotJS][resI][rotIR] = unifiedEmat[resI][rotIR][resJ][rotJS]
    return unifiedEmat
  
  def computeEboundForAllRots(self):
    for resI in range(self.numRes): 
      for rotIR in range(self.rotsPerRes[resI]):
        print "Res"+`resI`+"Rot"+`rotIR`+"="+`self.computeEboundForRot(resI, rotIR)`

  def computeEboundForRot(self, resI, rotIR):
    Ebound = 0.0
    for resJ in range(self.numRes):
      if(resJ != resI):
        Ebound += min(self.belief[resJ])
      else:
        Ebound += self.belief[resI][rotIR]
    return Ebound
