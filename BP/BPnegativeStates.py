from math import *
from ematrix.EMatrix import *
import util.mymath 
import ipdb
  
# Instances of this class contains the messages between all pairs of rotamers
class Messages:
  def __init__(self, ematrix, availableAAs, statename=""):
    self.maxProduct = True
    self.msgMat = {}
    self.newMsgMat = {}
    self.ematrix = ematrix
    self.availableAAs = availableAAs
    self.statename = statename
      
    for senderRes in ematrix.getSortedResIds(): 
      for senderAA in availableAAs[senderRes]:
        self.msgMat[(senderRes, senderAA)] = {}
        self.newMsgMat[(senderRes, senderAA)] = {}
        for recvRes in self.ematrix.getNeighbors(senderRes):
          for recvAA in availableAAs[recvRes]:
            self.newMsgMat[(senderRes,senderAA)][(recvRes,recvAA)] = {}
            self.msgMat[(senderRes,senderAA)][(recvRes,recvAA)] = {}
            for recvRot in ematrix.getRots(recvRes, recvAA):
              self.newMsgMat[(senderRes,senderAA)][(recvRes,recvAA)][recvRot] = 1.0
              self.msgMat[(senderRes,senderAA)][(recvRes,recvAA)][recvRot] = 1.0

  # resI and resJ are (res_id_tup,AA) tuples. rotR is an index (potentially a tuple).
  def getMsg_resI_to_resJrotR(self, resI, resJ, rotR):
    return self.msgMat[resI][resJ][rotR]

  # resJ is (res_id_tup,AA) tuples. rotR is an index (potentially a tuple)
  def get_all_msgs_to_resJrotR(self,resJ, rotR):
    msgs_to_resJrotR = {}
    res_id_tupJ = resJ[0]
    aaJ = resJ[1]
    for res_id_tupI in self.ematrix.getNeighbors(res_id_tupJ):
      for aaI in self.availableAAs[res_id_tupI]:
        msgs_to_resJrotR[(res_id_tupI,aaI)] = self.getMsg_resI_to_resJrotR((res_id_tupI,aaI),resJ,rotR)
    return msgs_to_resJrotR

  # resJ is (res_id_tup,AA) tuples. rotR is an index (potentially a tuple).
  def setMsg_resI_to_resJrotR(self, resI, resJ, rotR, msg):
    self.newMsgMat[resI][resJ][rotR] = msg

  # sendMsgsToAllNeighbors writes in the message matrix the message to every rotamer from "senderRes" (a tuple defined by both residue id and amino acid); that 
  # is, senderRes sends the messages from residue "i" to rotamer "r" at residue "j", for all "r" and "j"
  def sendMsgsToAllNeighbors(self, senderRes_id_tup, senderRes_aa):
    # recvRes: residue "j" to whom the message is directed to 
    for recvRes_id_tup in self.ematrix.getNeighbors(senderRes_id_tup):
      for recvAA in self.availableAAs[recvRes_id_tup]:
        # recvRot: the specific rotamer at which the message is sent to
        for recvRot in self.ematrix.getRots(recvRes_id_tup, recvAA):
          msg = 0.0
          # Our message will be the sum-product or max-product of all the messages from the rotamers at senderRes_id_tup,senderRes_aa
          for senderRot in self.ematrix.getRots(senderRes_id_tup, senderRes_aa):
            E = self.ematrix.getIntraE(senderRes_id_tup, senderRes_aa, senderRot) + \
                    self.ematrix.getPairE(senderRes_id_tup, senderRes_aa, senderRot, recvRes_id_tup, recvAA, recvRot)
            constant = exp(-E)
#            if constant < util.mymath.SMALL_FLOAT: 
              ## Energy is too high; setting to a small float
#              constant = util.mymath.SMALL_FLOAT

            prodIncMesgs = constant
            for resK_id_tup in self.ematrix.getNeighbors(senderRes_id_tup):
              for aaK in self.availableAAs[resK_id_tup]:
                if resK_id_tup != recvRes_id_tup or aaK != recvAA:
                  msg = self.getMsg_resI_to_resJrotR((resK_id_tup,aaK), (senderRes_id_tup, senderRes_aa), senderRot)
                  if msg > 0.0:
                    prodIncMesgs *= msg
            if(self.maxProduct and msg < prodIncMesgs):
              msg = prodIncMesgs
            # sum product BP
            else:
              if not self.maxProduct:
                msg += prodIncMesgs
          self.setMsg_resI_to_resJrotR((senderRes_id_tup,senderRes_aa), (recvRes_id_tup,recvAA), recvRot, msg)

  # Once all messages have been updated, we can commit the messages on a new matrix.  
  def commitMessages(self):
    # Compute the partition function matrix Z
    Z = {}
    for senderRes_id_tup in self.ematrix.getSortedResIds():
      for senderAA in self.availableAAs[senderRes_id_tup]: 
        Z[(senderRes_id_tup,senderAA)] = {}
        for recvRes_id_tup in self.ematrix.getNeighbors(senderRes_id_tup):
          for recvAA in self.availableAAs[recvRes_id_tup]:
            Z[(senderRes_id_tup,senderAA)][(recvRes_id_tup,recvAA)] = 0.0
            for recvRot in self.ematrix.getRots(recvRes_id_tup,recvAA):
              # For every pair of sender+receiver residues, Z[senderRes][recvRes] contains the summation of all the messages from 
              #  senderRes to all rotamers at recvRes
              Z[(senderRes_id_tup,senderAA)][(recvRes_id_tup,recvAA)] += \
                self.newMsgMat[(senderRes_id_tup,senderAA)][(recvRes_id_tup,recvAA)][recvRot]
    
    for senderRes_id_tup in self.ematrix.getSortedResIds(): 
      for senderAA in self.availableAAs[senderRes_id_tup]:
        for recvRes_id_tup in self.ematrix.getNeighbors(senderRes_id_tup):
          for recvAA in self.availableAAs[recvRes_id_tup]:
            for recvRot in self.ematrix.getRots(recvRes_id_tup, recvAA):
              # NOTE: If the partition function of a rotamer is zero, then we assume that all messages are zero. 
              myZ = Z[(senderRes_id_tup,senderAA)][(recvRes_id_tup,recvAA)] 
              if myZ > 0.0:
                self.msgMat[(senderRes_id_tup,senderAA)][(recvRes_id_tup,recvAA)][recvRot] = \
                    self.newMsgMat[(senderRes_id_tup,senderAA)][(recvRes_id_tup,recvAA)][recvRot] / \
                    Z[(senderRes_id_tup,senderAA)][(recvRes_id_tup,recvAA)]
              else:
                self.msgMat[(senderRes_id_tup,senderAA)][(recvRes_id_tup,recvAA)][recvRot] = 0.0

class BProtamer:
  myState = 1.0
  def __init__(self, myRes, myRot, myIntEnergy):
    self.myResAndType = myRes
    self.myRot = myRot
    self.myIntEnergy = myIntEnergy
  def updateState(self, myIncomingMsgs):
    self.myState = exp(-self.myIntEnergy)
    for res_id_and_type in myIncomingMsgs.keys():
      if res_id_and_type != self.myResAndType:
        myNeighborBeliefAboutMe = myIncomingMsgs[res_id_and_type]
        self.myState *= myNeighborBeliefAboutMe

# Probability of a pair of rotamers 
class BProtamerpair:
  myState = 1.0
  myPair={}
  def __init__(self, res_id1, rot1, intraE1, res_id2, rot2, intraE2, pairE):
    self.res_id1 = res_id1
    self.rot1 = rot1
    self.intraE1 = intraE1
    self.res_id2 = res_id2
    self.rot2 = rot2
    self.intraE2 = intraE2
    self.pairE = pairE
  def updateState(self, incomingMsgI, incomingMsgJ):
    self.myState = exp(-self.intraE1 - self.intraE2 - self.pairE)
    for res_idk in incomingMsgI.keys():
      if res_idk != self.res_id1 and res_idk != self.resid2:
        self.myState *= incomingMsgI[res_idk]*incomingMsgJ[res_idk]
  
# This class solves a variation of the protein design problem using Belief Propagation.
# The variation is that every (res_id,aa) becomes a new residue. For example, a problem where
# there were 4 residue positions, with 4 amino acids at each position and four rotamers for each 
# amino acid, now becomes a problem with 16 residues and 4 rotamers per position. The pairwise interactions
# are maintained.
# It implements max-product and sum-product.
#     sum-product: add all messages coming from all rotamers at i.  
#     max-product: send only the max-message
# ematrix is an object of type Ematrix.
class BPnegativeStatesSolver:
  def __init__(self, ematrix, assignedAAs, statename="", homodimer=False):
    # Initialize class that manages messages
    self.iterations = 50
    self.availableAAs = {}
    ## Compute the amino acids that will be allowed during the search.
    for res_id_tup in ematrix.getSortedResIds():
      self.availableAAs[res_id_tup] ={}
      # If there is only one amino acid allowed for this position in the ematrix
      if len(ematrix.getAllowedAAs(res_id_tup)) == 1:
        # Then assignedAAs is ignored. 
        # NOTE:This is a limitation of the current implementation that forces negative states 
        # to have either fixed AA identity or the same identities as the positive states.
        self.availableAAs[res_id_tup] = ematrix.getAllowedAAs(res_id_tup)
      elif res_id_tup in assignedAAs:  # AA has been assigned.
        self.availableAAs[res_id_tup] = assignedAAs[res_id_tup]
      else: # Otherwise assign all others.
        self.availableAAs[res_id_tup] = ematrix.getAllowedAAs(res_id_tup)

    self.messages = Messages(ematrix, self.availableAAs, statename)
    self.ematrix = ematrix
    # In Belief propagation protein design, each residue position is a variable
    # and it can be in several states, each of which represents a rotamer at that position.
    # Initialize variable states matrix, a (numRes x rotsPerRes) matrix that contains, 
    #  for each variable, the probability for each of its states.  Initially, all 
    #  states have a probability of 1. Note that here the residue index is (res_id_tup, aa).
    self.variableStates = {}
    for res_id_tup in self.ematrix.getSortedResIds():
      for aa in self.availableAAs[res_id_tup]:
        self.variableStates[(res_id_tup,aa)] = {}
        for rot in self.ematrix.getRots(res_id_tup, aa):
          self.variableStates[(res_id_tup,aa)][rot] = BProtamer((res_id_tup, aa), rot, self.ematrix.getIntraE(res_id_tup,aa,rot))

    self.pairVariableStates = {}
    for res_id_tup1 in self.ematrix.getSortedResIds():
      for aa1 in self.availableAAs[res_id_tup1]:
        self.pairVariableStates[(res_id_tup1,aa1)] = {}
        for rot1 in self.ematrix.getRots(res_id_tup1,aa1):
          self.pairVariableStates[(res_id_tup1,aa1)][rot1] = {}
          for res_id_tup2 in self.ematrix.getNeighbors(res_id_tup1):
            for aa2 in self.availableAAs[res_id_tup2]:
              if homodimer:
                resi2 = res_id_tup2[1]
                resi1 = res_id_tup1[1]
                if resi2 == resi1 and aa1 != aa2:
                    continue
              self.pairVariableStates[(res_id_tup1,aa1)][rot1][(res_id_tup2,aa2)] = {}
              for rot2 in self.ematrix.getRots(res_id_tup2,aa2):
                self.pairVariableStates[(res_id_tup1,aa1)][rot1][(res_id_tup2,aa2)][rot2] = \
                BProtamerpair((res_id_tup1,aa1), rot1, (res_id_tup2,aa2), rot2, \
                  self.ematrix.getIntraE(res_id_tup1,aa1,rot1), \
                  self.ematrix.getIntraE(res_id_tup2,aa2,rot2), \
                  self.ematrix.getPairE(res_id_tup1,aa1,rot1,res_id_tup2,aa2,rot2))
        
    
  def optimize(self):
    # The algorithm should loop until convergence.  For now we do it self.iterations times. 
    for i in range(self.iterations):
      for curRes_id_tup in self.ematrix.getSortedResIds():
        for curRes_aa in self.availableAAs[curRes_id_tup]:
          # We update the messages to everybody else for this residue.  
          self.messages.sendMsgsToAllNeighbors(curRes_id_tup, curRes_aa)
          
      # In the current version we update asynchronously.  
      self.messages.commitMessages()
  
      for curRes_id_tup in self.ematrix.getSortedResIds():
        for curRes_aa in self.availableAAs[curRes_id_tup]:
          for curRot in self.ematrix.getRots(curRes_id_tup, curRes_aa):
            # We first update the max marginal of the rotamer 
            self.variableStates[(curRes_id_tup, curRes_aa)][curRot].updateState(self.messages.get_all_msgs_to_resJrotR((curRes_id_tup, curRes_aa),curRot))

  def computeSubtrahend(self):

    pairwiseEnergy = 0.0
    internalEnergy = 0.0
    maxState = {}
    for res_id_tup in self.ematrix.getSortedResIds():
      maxState[res_id_tup] = {}
      for aa in self.availableAAs[res_id_tup]:
        maxValue = float("-inf")
        for rot in self.ematrix.getRots(res_id_tup, aa):
          if self.variableStates[(res_id_tup,aa)][rot].myState > maxValue:
            maxState[res_id_tup][aa] = rot
            maxValue = self.variableStates[(res_id_tup,aa)][rot].myState

    totalEnergy = 0.0
    subtrahendIntraMat = {}
    subtrahendPairMat = {}
    for res_id_tup1 in self.ematrix.getSortedResIds():
      subtrahendIntraMat[res_id_tup1] = {}
      for aa1 in self.availableAAs[res_id_tup1]:
        rot1 = maxState[res_id_tup1][aa1]
        subtrahendIntraMat[res_id_tup1][aa1] = \
          self.ematrix.getIntraE(res_id_tup1,aa1,rot1)
        totalEnergy += subtrahendIntraMat[res_id_tup1][aa1] 

    for res_id_tup1 in self.ematrix.getSortedResIds():
      subtrahendPairMat[res_id_tup1] = {}
      for res_id_tup2 in self.ematrix.getNeighbors(res_id_tup1):
        subtrahendPairMat[res_id_tup1][res_id_tup2] = {}
        for aa1 in self.availableAAs[res_id_tup1]:
          subtrahendPairMat[res_id_tup1][res_id_tup2][aa1] = {}
          for aa2 in self.availableAAs[res_id_tup2]:
            rot1 = maxState[res_id_tup1][aa1]
            rot2 = maxState[res_id_tup2][aa2]
            subtrahendPairMat[res_id_tup1][res_id_tup2][aa1][aa2] = \
              self.ematrix.getPairE(res_id_tup1,aa1,rot1,res_id_tup2,aa2,rot2)
            totalEnergy += subtrahendPairMat[res_id_tup1][res_id_tup2][aa1][aa2]/2.0

    print "Total energy of representative rotamers in negative state is: "+`totalEnergy`

    self.subtrahendIntraMat = subtrahendIntraMat
    self.subtrahendPairMat = subtrahendPairMat
    return subtrahendIntraMat, subtrahendPairMat

#    print "pairwiseEnergyInteger = "+`pairwiseEnergy`
#    print "internalEnergyInteger = "+`internalEnergy`   
#    print "totalEnergyInteger= "+`pairwiseEnergy+internalEnergy`   

  # This function returns the res_id_tup of the residue that should be expanded next.
  # This is a heuristic for the dynamic expansion of residues. 
  def computeNextResidueToExpandBasedOnVariableStatesRMS(self):
    # Compute a partition function for each state
    Z = {}
    for res_id_tup in self.ematrix.getSortedResIds():
      for aa in self.availableAAs[res_id_tup]:
        Z[(res_id_tup,aa)] = 0.0
        for rot in self.ematrix.getRots(res_id_tup,aa):
          Z[(res_id_tup,aa)] += self.variableStates[(res_id_tup,aa)][rot].myState
    
    res_id_with_min_rms = ""
    min_state_rms = float("inf")
    for res_id_tup in self.ematrix.getSortedResIds():
      for aa in self.availableAAs[res_id_tup]:
        stateValues = [x.myState/Z[(res_id_tup,aa)] for x in self.variableStates[(res_id_tup,aa)]] 
        rms = mymath.rootMeanSquaresForAll(stateValues)
        if rms < min_state_rms: 
          res_id_with_min_rms = res_id_tup
          min_state_rms = rms

    return res_id_with_min_rms

  ## For debugging purposes
#  def getSubtrahendEnergyOfSequence(self, sequence):
#    energy = 0.0
#    for res_id_tup1 in self.ematrix.getSortedResIds():
#      energy += self.subtrahendIntraMat[res_id_tup1][sequence[res_id_tup1]]
#      for res_id_tup2 in self.ematrix.getNeighbors(res_id_tup1):
#        energy += self.subtrahendPairMat[res_id_tup1][res_id_tup2][sequence[res_id_tup1]][sequence[res_id_tup2]]/2.0
#    return energy

  # Also for debugging
  def getLMECandEnergy(self, sequence):
    maxState = {}
    totalEnergy = 0.0
    for res_id_tup in self.ematrix.getSortedResIds():
      maxState[res_id_tup] = {}
      aa = sequence[res_id_tup]
      maxValue = float("-inf")
      for rot in self.ematrix.getRots(res_id_tup, aa):
        if self.variableStates[(res_id_tup,aa)][rot].myState > maxValue:
          maxState[res_id_tup][aa] = rot
          maxValue = self.variableStates[(res_id_tup,aa)][rot].myState

    totalEnergy = 0.0
    for res_id_tup1 in self.ematrix.getSortedResIds():
      aa1 = sequence[res_id_tup1]
      rot1 = maxState[res_id_tup1][aa1]
      print "Internal energy of "+`(res_id_tup1,aa1,rot1)`+"="+`self.ematrix.getIntraE(res_id_tup1,aa1,rot1)`
      totalEnergy += self.ematrix.getIntraE(res_id_tup1,aa1,rot1)

    for res_id_tup1 in self.ematrix.getSortedResIds():
      for res_id_tup2 in self.ematrix.getNeighbors(res_id_tup1):
        aa1 = sequence[res_id_tup1]
        aa2 = sequence[res_id_tup2]
        rot1 = maxState[res_id_tup1][aa1]
        rot2 = maxState[res_id_tup2][aa2]
        print "Pair energy of "+`((res_id_tup1,aa1,rot1),(res_id_tup2,aa2,rot2))`+"="+\
            `self.ematrix.getPairE(res_id_tup1,aa1,rot1,res_id_tup2,aa2,rot2)`
        totalEnergy += self.ematrix.getPairE(res_id_tup1,aa1,rot1,res_id_tup2,aa2,rot2)/2.0

    print "Total energy of representative rotamers in negative state is: "+`totalEnergy`


#    print "pairwiseEnergyInteger = "+`pairwiseEnergy`
#    print "internalEnergyInteger = "+`internalEnergy`   
#    print "totalEnergyInteger= "+`pairwiseEnergy+internalEnergy`   

  # This function returns the res_id_tup of the residue that should be expanded next.
  # This is a heuristic for the dynamic expansion of residues. 
  def computeNextResidueToExpandBasedOnVariableStatesRMS(self):
    # Compute a partition function for each state
    Z = {}
    for res_id_tup in self.ematrix.getSortedResIds():
      for aa in self.availableAAs[res_id_tup]:
        Z[(res_id_tup,aa)] = 0.0
        for rot in self.ematrix.getRots(res_id_tup,aa):
          Z[(res_id_tup,aa)] += self.variableStates[(res_id_tup,aa)][rot].myState
    
    res_id_with_min_rms = ""
    min_state_rms = float("inf")
    for res_id_tup in self.ematrix.getSortedResIds():
      for aa in self.availableAAs[res_id_tup]:
        stateValues = [x.myState/Z[(res_id_tup,aa)] for x in self.variableStates[(res_id_tup,aa)]] 
        rms = mymath.rootMeanSquaresForAll(stateValues)
        if rms < min_state_rms: 
          res_id_with_min_rms = res_id_tup
          min_state_rms = rms

    return res_id_with_min_rms

  ## For debugging purposes
  def getSubtrahendEnergyOfSequence(self, sequence):
    energy = 0.0
    for res_id_tup1 in self.ematrix.getSortedResIds():
      energy += self.subtrahendIntraMat[res_id_tup1][sequence[res_id_tup1]]
      for res_id_tup2 in self.ematrix.getNeighbors(res_id_tup1):
        energy += self.subtrahendPairMat[res_id_tup1][res_id_tup2][sequence[res_id_tup1]][sequence[res_id_tup2]]/2.0
    return energy
  ## For debugging purposes
#  def getLMECandEnergy(self, sequence):
#    energy = 0.0
#    for res_id_tup1 in self.ematrix.getSortedResIds():
#      energy += self.subtrahendIntraMat[res_id_tup1][sequence[res_id_tup1]]
#      for res_id_tup2 in self.ematrix.getNeighbors(res_id_tup1):
#        energy += self.subtrahendPairMat[res_id_tup1][res_id_tup2][sequence[res_id_tup1]][sequence[res_id_tup2]]/2.0
#    return energy
#
 #   energy = 0.0
#    for res_id_tup1 in self.ematrix.getSortedResIds():
#      energy += self.subtrahendIntraMat[res_id_tup1][sequence[res_id_tup1]]
#      for res_id_tup2 in self.ematrix.getNeighbors(res_id_tup1):
#        energy += self.subtrahendPairMat[res_id_tup1][res_id_tup2][sequence[res_id_tup1]][sequence[res_id_tup2]]/2.0
#    return energy

