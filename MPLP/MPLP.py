# This class computes a lower bound on the energy of a protein design problem, based on 
# Edge MPLP (NMPLP)  by globerson and Jaakkola
# Initialize messages to any value
# Iterate until convergence
from math import *
import random
from util.createMatrix import *
import ipdb
from util.mymath import *
import config.config
import time
class MPLP:
  debug = False
  doPruning = False
  # subtrahend is a dict where the first dimension keys are the negative state from which the terms come from. 
  # The second field is either 'intra' or 'pair'. The 3-5 fields are residue and aa or residue1, residue2, aa1, aa2
  def __init__(self, ematrix, subtrahend={}):
    self.ematrix = ematrix
    self.subtrahend = subtrahend

  # Computes the low-energy bound using the EMPLP algorithm
  # availableRots is a "numRes*rots" matrix that tells us which 
  # rotamers at each position are available.  Only those rotamers
  # will be considered for the calculation.
  # Also checks if the solution is integer
  def optimizeEMPLP(self, assignedAAs={}, assignedRots={}):
    availableRots = {}
    totalRots = 0
    ## Compute availableAAs
    availableAAs = {}
    ## Compute the amino acids that will be allowed during the search.
    for res_id_tup in self.ematrix.getSortedResIds():
      availableAAs[res_id_tup] = {}
      # If there is only one amino acid allowed for this position in the ematrix
      if len(self.ematrix.getAllowedAAs(res_id_tup)) == 1:
        # Then assignedAAs is ignored. 
        # NOTE:This is a limitation of the current implementation that forces negative states 
        # to have either fixed AA identity or the same identities as the positive states.
        availableAAs[res_id_tup] = self.ematrix.getAllowedAAs(res_id_tup)
      elif res_id_tup in assignedAAs:  # AA has been assigned.
        availableAAs[res_id_tup] = assignedAAs[res_id_tup]
      else: # Otherwise assign all others.
        availableAAs[res_id_tup] = self.ematrix.getAllowedAAs(res_id_tup)


    for res_id_tupI in self.ematrix.getSortedResIds():
      availableRots[res_id_tupI] = {}
      if res_id_tupI in assignedRots: 
        assert len(availableAAs[res_id_tupI]) == 1
        for aaI in availableAAs[res_id_tupI]:
          availableRots[res_id_tupI][aaI] = [assignedRots[res_id_tupI][aaI]]
          totalRots += 1
      else:
        for aaI in availableAAs[res_id_tupI]:
          availableRots[res_id_tupI][aaI] = self.ematrix.getRots(res_id_tupI, aaI)
          totalRots+= len(availableRots[res_id_tupI][aaI])
    # Prune space with DEE to accelerate MPLP.
    isRotPruned = self.DEE(availableAAs, availableRots, self.doPruning)
    availableRots = {}
    totalRots =0
    for res_id_tupI in self.ematrix.getSortedResIds():
      availableRots[res_id_tupI] = {}
      for aa in isRotPruned[res_id_tupI]:
        availableRots[res_id_tupI][aa] = []
        for rot in isRotPruned[res_id_tupI][aa]:
          if not isRotPruned[res_id_tupI][aa][rot]:
            availableRots[res_id_tupI][aa].append(rot)
            totalRots += 1

    for res_id_tupI in self.ematrix.getSortedResIds():
      availableRots[res_id_tupI] = {}
      if res_id_tupI in assignedRots: 
        for aaI in availableAAs[res_id_tupI]:
          availableRots[res_id_tupI][aaI] = [assignedRots[res_id_tupI][aaI]]
      else:
        for aaI in availableAAs[res_id_tupI]:
          availableRots[res_id_tupI][aaI] = self.ematrix.getRots(res_id_tupI, aaI)
    for res_id_tupI in self.ematrix.getSortedResIds():
      availableRots[res_id_tupI] = {}
      for aaI in availableAAs[res_id_tupI]:
        if res_id_tupI in assignedRots and aaI in assignedRots[res_id_tupI]:
          availableRots[res_id_tupI][aaI] = [assignedRots[res_id_tupI][aaI]]
        else:
          availableRots[res_id_tupI][aaI] = self.ematrix.getRots(res_id_tupI, aaI)
    # l is the message matrix for EMPLP; we set all messages to zero
    self.l = {}
    for res_id_tupI in self.ematrix.getSortedResIds():
      self.l[res_id_tupI] = {}
      for res_id_tupJ in self.ematrix.getNeighbors(res_id_tupI):
        self.l[res_id_tupI][res_id_tupJ] = {}
        for aaJ in availableAAs[res_id_tupJ]:
          for rotJS in availableRots[res_id_tupJ][aaJ]:
            self.l[res_id_tupI][res_id_tupJ][(aaJ,rotJS)] = 0.0

    # The belief on each rotamer
    self.belief = {}
    for res_id_tupJ in self.ematrix.getSortedResIds():
      self.belief[res_id_tupJ] = {}
      for aaJ in availableAAs[res_id_tupJ]:
        for rotJS in availableRots[res_id_tupJ][aaJ]:
          self.belief[res_id_tupJ][(aaJ,rotJS)] = 0.0
    
    
    # EMPLP algorithm.  
    # The algorithm should loop until convergence.  For now we do it I times. 
    # Complexity of EMPLP: O(I*numRes*numRes*rotsPerRes*rotsPerRes*numRes) 
    #for iteration in range(config.config.MPLP_ITERATIONS):
    for iteration in range(1):
      for res_id_tupI in self.ematrix.getSortedResIds():
        for res_id_tupJ in self.ematrix.getNeighbors(res_id_tupI):
          for aaI in availableAAs[res_id_tupI]:
            # We first update self.l[resJ][resI][rotIR] and immediately after we update self.l[resI][resJ][rotJS] 
            for rotIR in availableRots[res_id_tupI][aaI]:
              self.belief[res_id_tupI][(aaI,rotIR)] -= self.l[res_id_tupJ][res_id_tupI][(aaI,rotIR)]
              incMesgsToRotIR_fromAllExcept_resJ = self.belief[res_id_tupI][(aaI,rotIR)] # - self.l[resJ][resI][rotIR] We already substracted the message from resJ.
              msgsFromRotsAtJ_to_rotIR = []
              for aaJ in availableAAs[res_id_tupJ]:
                for rotJS in availableRots[res_id_tupJ][aaJ]:
                  msgsFromRotsAtJ_to_rotIR.append(self.belief[res_id_tupJ][(aaJ,rotJS)] - self.l[res_id_tupI][res_id_tupJ][(aaJ,rotJS)] \
                        + self.unifiedEmat(res_id_tupI, aaI, rotIR, res_id_tupJ, aaJ, rotJS, availableAAs))

              self.l[res_id_tupJ][res_id_tupI][(aaI,rotIR)] = -0.5*incMesgsToRotIR_fromAllExcept_resJ + 0.5*min(msgsFromRotsAtJ_to_rotIR)
              self.belief[res_id_tupI][(aaI,rotIR)] += self.l[res_id_tupJ][res_id_tupI][(aaI,rotIR)]

          for aaJ in availableAAs[res_id_tupJ]:
            # Now we update self.l[resI][resJ][rotJS] 
            for rotJS in availableRots[res_id_tupJ][aaJ]:
              self.belief[res_id_tupJ][(aaJ,rotJS)] -= self.l[res_id_tupI][res_id_tupJ][(aaJ,rotJS)]
              incMesgsToRotJS_fromAllExcept_resI = self.belief[res_id_tupJ][(aaJ,rotJS)]
              msgsFromRotsAtI_to_rotJS = []
              for aaI in availableAAs[res_id_tupI]:
                for rotIR in availableRots[res_id_tupI][aaI]:
                  msgsFromRotsAtI_to_rotJS.append((self.belief[res_id_tupI][(aaI,rotIR)] - self.l[res_id_tupJ][res_id_tupI][(aaI,rotIR)] +\
                    self.unifiedEmat(res_id_tupI,aaI,rotIR,res_id_tupJ,aaJ,rotJS, availableAAs)))

              self.l[res_id_tupI][res_id_tupJ][(aaJ,rotJS)] = (-0.5)*incMesgsToRotJS_fromAllExcept_resI + 0.5*min(msgsFromRotsAtI_to_rotJS)
              self.belief[res_id_tupJ][(aaJ,rotJS)] += self.l[res_id_tupI][res_id_tupJ][(aaJ,rotJS)]
    Ebound = 0.0
    EcontributionPerResidue = {}
    gmec = {}
    # Compute ebound and gmec
    for res_id_tupI in self.ematrix.getSortedResIds():
      minForResI = float("inf")
      for aaI in availableAAs[res_id_tupI]:
        for rotIR in availableRots[res_id_tupI][aaI]:
          if self.belief[res_id_tupI][(aaI,rotIR)] < minForResI:
            minForResI = self.belief[res_id_tupI][(aaI,rotIR)]
            gmec[res_id_tupI] = (aaI,rotIR)
      EcontributionPerResidue[res_id_tupI] = minForResI
      Ebound += minForResI
    
    if self.subtrahend != {}:
      for statename in self.subtrahend.keys():
        subtrahend_residual = self.getSubtrahendResidual(statename)
        Ebound -= subtrahend_residual
    return Ebound, EcontributionPerResidue

  # Returns a single term from the ematrix containing part of the intra terms and the pairwise term
  def unifiedEmat(self, res_id_tupI, aaI, rotIR, res_id_tupJ, aaJ, rotJS, availableAAs):
    unifE = self.ematrix.getPairE(res_id_tupI, aaI, rotIR, res_id_tupJ, aaJ, rotJS)
    neighborsI = float(len(self.ematrix.getNeighbors(res_id_tupI)))
    neighborsJ = float(len(self.ematrix.getNeighbors(res_id_tupJ)))
    unifE += self.ematrix.getIntraE(res_id_tupI, aaI, rotIR)/neighborsI
    unifE += self.ematrix.getIntraE(res_id_tupJ, aaJ, rotJS)/neighborsJ
    # The subtrahends have been defined for this problem.
    if self.subtrahend != {}:
      for statename in self.subtrahend.keys():
        if res_id_tupI in self.subtrahend[statename]['intra'] and aaI in self.subtrahend[statename]['intra'][res_id_tupI]:
          unifE -= self.subtrahend[statename]['intra'][res_id_tupI][aaI]/neighborsI
          unifE -= self.getSubtrahendResidual(statename, res_id_tupI, aaI)/neighborsI
        if res_id_tupJ in self.subtrahend[statename]['intra'] and aaJ in self.subtrahend[statename]['intra'][res_id_tupJ]:
          unifE -= self.subtrahend[statename]['intra'][res_id_tupJ][aaJ]/neighborsJ
          unifE -= self.getSubtrahendResidual(statename, res_id_tupJ, aaJ)/neighborsJ
        if res_id_tupI in self.subtrahend[statename]['pair'] and res_id_tupJ in self.subtrahend[statename]['pair'][res_id_tupI] \
            and aaI in self.subtrahend[statename]['intra'][res_id_tupI] and  aaJ in self.subtrahend[statename]['intra'][res_id_tupJ]:
          unifE -= self.subtrahend[statename]['pair'][res_id_tupI][res_id_tupJ][aaI][aaJ]

    return unifE

  # Compute a new residue to expand based on the belief of the MPLP optimization.
  # Returns the residue with the lowest RMS on its belief. 
  def getHeuristicNextResidueToExpandBasedOnBeliefs(self, assignedRots={}):
    min_rms = float("inf")
    min_rms_res_id_tup = ""
    for res_id_tup in self.ematrix.getSortedResIds():
      if res_id_tup not in assignedRots:
        Z=0
        # First compute a "partition function" to normalize values for all residues between 0 and 1.
        for rot in self.belief[res_id_tup].keys():
          Z += self.belief[res_id_tup][rot]
        if Z == 0: 
            continue
        rms_for_this_residue = 0
        for rot in self.belief[res_id_tup].keys():
          rms_for_this_residue += (self.belief[res_id_tup][rot]/Z) ** 2
        rms_for_this_residue = rms_for_this_residue ** 0.5
        if rms_for_this_residue < min_rms:
          min_rms = rms_for_this_residue
          min_rms_res_id_tup = res_id_tup
    # If no residue was selected, then return any residue. 
    if min_rms_res_id_tup == "":
      for res_id_tup in self.ematrix.getSortedResIds():
        if res_id_tup not in assignedRots:
          min_rms_res_id_tup = res_id_tup
    return min_rms_res_id_tup
     
#  def getBestSequenceAccordingToBeliefs:
#    for res_id_tup in self.ematrix.getSortedResIds():


  # Some residues in the subtrahend matrices may not be included in the ematrix. 
  # For positive/negative design, this energy has to be accounted in some way. 
  #   This function computes that value. 
  # If res_id_tup is a specific residue in self.ematrix, then this function returns the sum of 
  #   pairwise energies between res_id_tup and its neighbors that are not in self.ematrix.
  # If res_id_tup is not defined, then this function returns the subtrahend energies and pairwise
  # energies between all residues that are not in self.ematrix.
  def getSubtrahendResidual(self, statename, res_id_tupI=None, aaI=None):
    subtrahend_residual_energy = 0.0
    if res_id_tupI != None:
      for res_id_tupJ in self.subtrahend[statename]['pair'][res_id_tupI].keys():
        if res_id_tupJ not in self.ematrix.getSortedResIds(): ## If J is not defined in the MPLP matrix.
          for aaJ in self.subtrahend[statename]['pair'][res_id_tupI][res_id_tupJ][aaI].keys():
            subtrahend_residual_energy += self.subtrahend[statename]['pair'][res_id_tupI][res_id_tupJ][aaI][aaJ]
    else:
      for res_id_tupI in self.subtrahend[statename]['intra'].keys():
        if res_id_tupI not in self.ematrix.getSortedResIds(): ## Only residues that are not present in the ematrix state
          for aaI in self.subtrahend[statename]['intra'][res_id_tupI].keys():
            subtrahend_residual_energy += self.subtrahend[statename]['intra'][res_id_tupI][aaI]
            for res_id_tupJ in self.subtrahend[statename]['pair'][res_id_tupI].keys():
              if res_id_tupJ not in self.ematrix.getSortedResIds(): ## And their interaction to other residues not in the ematrix state
                assert(len(self.subtrahend[statename]['pair'][res_id_tupI][res_id_tupJ][aaI].keys()) == 1)
                for aaJ in self.subtrahend[statename]['pair'][res_id_tupI][res_id_tupJ][aaI].keys():
                  subtrahend_residual_energy += self.subtrahend[statename]['pair'][res_id_tupI][res_id_tupJ][aaI][aaJ]/2.0

    return subtrahend_residual_energy

  def DEE(self, availableAAs, availableRots, doPruning=True):
    isRotPruned = {}
    # Initialize the pruning matrix to false.
    countTotalRots = 0
    
    for res_id_tup in self.ematrix.getSortedResIds():
      isRotPruned[res_id_tup] = {}
      for aa in availableAAs[res_id_tup]:
        isRotPruned[res_id_tup][aa] = {}
        for rot in availableRots[res_id_tup][aa]:
          isRotPruned[res_id_tup][aa][rot] = False
          countTotalRots +=1 
    
    rotamersPrunedThisRound = 1
    while doPruning and rotamersPrunedThisRound > 0:
      rotamersPrunedThisRound = 0
      
      # Candidate is the rotamer to be pruned; witness is the rotamer that will prune it.
      for res_id_tup in self.ematrix.getSortedResIds():
        for aa2 in availableAAs[res_id_tup]:
          for witness_rot in availableRots[res_id_tup][aa2]:
            if not isRotPruned[res_id_tup][aa2][witness_rot]:
              for aa1 in availableAAs[res_id_tup]:
                for candidate_rot in availableRots[res_id_tup][aa1]:
                  if not isRotPruned[res_id_tup][aa1][candidate_rot]:
                    if witness_rot != candidate_rot or aa1 != aa2:
                      if self.can_prune_rot(res_id_tup, aa1, candidate_rot, aa2, witness_rot, availableAAs, availableRots, isRotPruned):
                        isRotPruned[res_id_tup][aa1][candidate_rot] = True
                        #print "Pruned rotamer "+`res_id_tup`+" "+aa1+" "+`candidate_rot`
                        rotamersPrunedThisRound += 1
#      print "Rotamers pruned this round "+`rotamersPrunedThisRound`
    return isRotPruned

  # DEE: Can prune a candidate rotamer with respect to a witness?
  def can_prune_rot(self, res_id_tup, aa1, candidate_rot, aa2, witness_rot, availableAAs, availableRots, isRotPruned):
    pruning_interval = self.ematrix.getIntraE(res_id_tup,aa1,candidate_rot) - self.ematrix.getIntraE(res_id_tup,aa2,witness_rot)
    if self.subtrahend != {}:
      # account for the intervals in each negative state.
      for statename in self.subtrahend.keys():
        if res_id_tup in self.subtrahend[statename]['intra'].keys():
          subtrahend_interval1 = 0.0
          subtrahend_interval2 = 0.0
          if aa1 in self.subtrahend[statename]['intra'][res_id_tup].keys():
            subtrahend_interval1 = -self.subtrahend[statename]['intra'][res_id_tup][aa1]
            subtrahend_interval1 -= self.getSubtrahendResidual(statename, res_id_tup, aa1)
          if aa2 in self.subtrahend[statename]['intra'][res_id_tup].keys():
            subtrahend_interval2 = -self.subtrahend[statename]['intra'][res_id_tup][aa2]
            subtrahend_interval2 -= self.getSubtrahendResidual(statename, res_id_tup, aa2)
          pruning_interval = pruning_interval + subtrahend_interval1 - subtrahend_interval2

    for otherRes_id_tup in self.ematrix.getNeighbors(res_id_tup):
      minimumForThisResidue = float("inf")
      # go through each rotamer in otherRes for the allowed aminoacids in otherRes
      for otherAA in availableAAs[otherRes_id_tup]:
        for otherRot in availableRots[otherRes_id_tup][otherAA]:
          if not isRotPruned[otherRes_id_tup][otherAA][otherRot]:
            candidateEnergyWithOtherRot = self.ematrix.getPairE(res_id_tup, aa1, candidate_rot, otherRes_id_tup,otherAA,otherRot)
            witnessEnergyWithOtherRot = self.ematrix.getPairE(res_id_tup, aa2, witness_rot, otherRes_id_tup, otherAA, otherRot)
            if self.subtrahend != {}:
              for statename in self.subtrahend.keys():
                if res_id_tup in self.subtrahend[statename]['intra'] and otherRes_id_tup in self.subtrahend[statename]['intra'] \
                    and aa1 in self.subtrahend[statename]['intra'][res_id_tup] and otherAA in self.subtrahend[statename]['intra'][otherRes_id_tup]:
                  candidateEnergyWithOtherRot -= self.subtrahend[statename]['pair'][res_id_tup][otherRes_id_tup][aa1][otherAA]
                if res_id_tup in self.subtrahend[statename]['intra'] and otherRes_id_tup in self.subtrahend[statename]['intra'] \
                    and aa2 in self.subtrahend[statename]['intra'][res_id_tup] and otherAA in self.subtrahend[statename]['intra'][otherRes_id_tup]:
                  witnessEnergyWithOtherRot -= self.subtrahend[statename]['pair'][res_id_tup][otherRes_id_tup][aa2][otherAA]
            diff = candidateEnergyWithOtherRot - witnessEnergyWithOtherRot
            minimumForThisResidue = min(minimumForThisResidue, diff)
      pruning_interval += minimumForThisResidue
    if pruning_interval >= 0.0:
      return True
    else:
      return False
