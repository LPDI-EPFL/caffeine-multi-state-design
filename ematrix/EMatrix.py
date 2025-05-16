from ematrix.parseRosettaMatrixDirectory import *
from ematrix.ReducedEmatrix import *
from util.createMatrix import * 
from config import config
class EMatrix:
  doPruning=True
  # Initialize energy matrix to null or by merging two energy matrices.
  def __init__(self, mergingMatX=None, mergingMatY=None):
    self.singleBody = {}
    self.pairwiseE = {}
    if mergingMatX != None and mergingMatY != None: 
      self.initializeEmatByMergingMatrices(mergingMatX, mergingMatY)
      # Now prune
      self.tdDEE(config.DO_TD_DEE)

  def parseMatrix(self, rosetta_matrix_directory):
    self.singleBody, self.pairwiseE = parseRosettaMatrixDirectory(rosetta_matrix_directory)
    self.numRes = len(self.singleBody.keys())
    for res in self.singleBody: 
      assert('_order' in self.singleBody[res])
      assert('_map_res_rot_to_aa_rot' in self.singleBody[res])
    self.tdDEE(config.DO_TD_DEE)

  # Merge matrices X and Y by making rotamer tuples out of all rotamers that belong to the same amino acid/residue 
  def initializeEmatByMergingMatrices(self, mergingMatX, mergingMatY):
    # First merge the single body terms.
    for res_id_tup1X in mergingMatX.getSortedResIds():
      self.singleBody[res_id_tup1X] = {}
      self.singleBody[res_id_tup1X]['_order'] = []
      for aa1X in mergingMatX.getAllowedAAs(res_id_tup1X):
        self.singleBody[res_id_tup1X]['_order'].append(aa1X)
        self.singleBody[res_id_tup1X][aa1X] = {}
        # Check if the res,aa appears in both matrices.
        if res_id_tup1X in mergingMatY.getSortedResIds() and aa1X in mergingMatY.getAllowedAAs(res_id_tup1X):
          # Then we must merge rotamers
          for rot1X in mergingMatX.getRots(res_id_tup1X, aa1X):
            for rot1Y in mergingMatY.getRots(res_id_tup1X, aa1X):
              self.singleBody[res_id_tup1X][aa1X][(rot1X, rot1Y)] = {'E': (mergingMatX.getIntraE(res_id_tup1X, aa1X, rot1X) + mergingMatY.getIntraE(res_id_tup1X, aa1X, rot1Y)), 'source' : 'XY'}
        else:
          for rot1X in mergingMatX.getRots(res_id_tup1X, aa1X):
            self.singleBody[res_id_tup1X][aa1X][rot1X] = {'E': mergingMatX.getIntraE(res_id_tup1X, aa1X, rot1X), 'source': 'X'}
    # Repeat the process with matrix Y, but don't merge rotamers as they have already been merged. 
    for res_id_tup1Y in mergingMatY.getSortedResIds():
      if res_id_tup1Y not in self.getSortedResIds():
        self.singleBody[res_id_tup1Y] = {}
        self.singleBody[res_id_tup1Y]['_order'] = []
      for aa1Y in mergingMatY.getAllowedAAs(res_id_tup1Y):
        # Only add rotamers for AA not in X
        if aa1Y not in self.getAllowedAAs(res_id_tup1Y):
          self.singleBody[res_id_tup1Y]['_order'].append(aa1Y)
          self.singleBody[res_id_tup1Y][aa1Y] = {}
        for rot1Y in mergingMatY.getRots(res_id_tup1Y, aa1Y):
          self.singleBody[res_id_tup1Y][aa1Y][rot1Y] = {'E': mergingMatY.getIntraE(res_id_tup1Y, aa1Y, rot1Y), 'source': 'Y'}
    # Merge the pairwise terms. 
    for res_id_tup1 in self.getSortedResIds():
      self.pairwiseE[res_id_tup1] = {}
      for res_id_tup2 in self.getSortedResIds():
        # res_id_tup2 has to be a neighbor of res_id_tup1 in either X or Y
        if res_id_tup1 in mergingMatX.getSortedResIds() and res_id_tup2 in mergingMatX.getNeighbors(res_id_tup1) or \
            res_id_tup1 in mergingMatY.getSortedResIds() and res_id_tup2 in mergingMatY.getNeighbors(res_id_tup1):
          self.pairwiseE[res_id_tup1][res_id_tup2] = {}
          for aa1 in self.getAllowedAAs(res_id_tup1):
            self.pairwiseE[res_id_tup1][res_id_tup2][aa1] = {}
            for aa2 in self.getAllowedAAs(res_id_tup2):
              self.pairwiseE[res_id_tup1][res_id_tup2][aa1][aa2] = {}
              for rot1 in self.singleBody[res_id_tup1][aa1].keys():
                self.pairwiseE[res_id_tup1][res_id_tup2][aa1][aa2][rot1] = {}
                for rot2 in self.singleBody[res_id_tup2][aa2].keys():
                  energy = 0.0
                  source1 = self.singleBody[res_id_tup1][aa1][rot1]['source']
                  rot1x = rot1y = rot1
                  if source1 == 'XY':
                    rot1x, rot1y = rot1
                  source2 = self.singleBody[res_id_tup2][aa2][rot2]['source']
                  rot2x = rot2y = rot2
                  if source2 == 'XY':
                    rot2x, rot2y = rot2
                  
                  if (source1 == 'X' and source2 == 'Y') or (source1 == 'Y' and source2 =='X'):
                    energy = 0.0
                  elif (source1 == 'X' and source2 == 'X') or (source1 == 'XY' and source2 =='X') or (source1 == 'X' and source2 =='XY'):  
                    energy = mergingMatX.getPairE(res_id_tup1, aa1, rot1x, res_id_tup2, aa2, rot2x)
                  elif (source1 == 'Y' and source2 == 'Y') or (source1 == 'XY' and source2 =='Y') or (source1 == 'Y' and source2 =='XY'):  
                    energy = mergingMatY.getPairE(res_id_tup1, aa1, rot1y, res_id_tup2, aa2, rot2y)
                  elif (source1 == 'XY' and source2 == 'XY'):
                    energy = mergingMatX.getPairE(res_id_tup1, aa1, rot1x, res_id_tup2, aa2, rot2x) + mergingMatY.getPairE(res_id_tup1, aa1, rot1y, res_id_tup2, aa2, rot2y)

                  self.pairwiseE[res_id_tup1][res_id_tup2][aa1][aa2][rot1][rot2] = energy
            
  def getSortedResIds(self):
    return sorted(self.singleBody.keys())

  def getAllowedAAs(self, res_id_tup=None):
    if res_id_tup == None:
      allowedAAs = {}
      for res_id_tup in self.getSortedResIds():
        allowedAAs[res_id_tup] = list (self.singleBody[res_id_tup]['_order'])
      return allowedAAs
    return list(self.singleBody[res_id_tup]['_order'])

  def getIntraE(self, res_id_tup, aa, rot):
    return self.singleBody[res_id_tup][aa][rot]['E']

  def getPairE(self, res_id_tup1, aa1, rot1, res_id_tup2, aa2, rot2):
    ## If two residues are not neighbors, return 0. 
    if res_id_tup2 in self.pairwiseE[res_id_tup1]:
      return self.pairwiseE[res_id_tup1][res_id_tup2][aa1][aa2][rot1][rot2]
    else:
      return 0.0
  def getDihedrals (self, res_id_tup, aa, rot):
    return self.singleBody[res_id_tup][aa][rot]['Dih']
  
  # Returns the number of rotamers for a specific amino acid or all by default.
  def getRots(self, res_id_tup, aa):
    nonPrunedRots = []
    for rot in self.singleBody[res_id_tup][aa].keys():
      if not self.isRotPruned[res_id_tup][aa][rot]:
        nonPrunedRots.append(rot)
    return nonPrunedRots

  def getNeighbors(self, res_id_tup):
    if res_id_tup in self.pairwiseE:
      return sorted(self.pairwiseE[res_id_tup].keys())
    else:
      return []

  def areNeighbors(self, res_id_tup1, res_id_tup2):
    if res_id_tup2 in self.pairwiseE[res_id_tup1]:
      return True
    else:
      return False

  def getRotsPerRes(self, seq=[]):
    rots_per_res_2D = [0]*self.numRes
    map_resNums_to_residue_ids = self.getSortedResIds()
    for res_ix in range(self.numRes):
      res_id_tup = map_resNums_to_residue_ids[res_ix]
      AAs = []
      if seq[res_ix] == 'XXX':
        AAs = self.singleBody[res_id_tup]['_order']
      else:
        AAs = [seq[res_ix]]
      for aa in AAs:
        for rot in range(len(self.singleBody[res_id_tup][aa])):
          if not self.pruneMatrix_rotamers[res_id_tup][aa][rot]:
            rots_per_res_2D[res_ix] +=1
      assert(rots_per_res_2D[res_ix] > 0)
    return rots_per_res_2D

  # Type-dependent DEE: Can prune a rotamer? 
  def can_prune_rot(self, res_id_tup, aa1, candidate_rot, aa2, witness_rot):
    e1 = self.singleBody[res_id_tup][aa1][candidate_rot]['E']
    e2 = self.singleBody[res_id_tup][aa2][witness_rot]['E']
    self.getNeighbors(res_id_tup)
    for otherRes_id_tup in self.getNeighbors(res_id_tup):
      minimumForThisResidue = float("inf")
      e1_add = 0.0
      e2_add =0.0
      # go through each rotamer in otherRes for the allowed aminoacids in otherRes
      for otherAA in self.getAllowedAAs(otherRes_id_tup):
        for otherRot in self.getRots(otherRes_id_tup, otherAA):
          candidateEnergyWithOtherRot = self.pairwiseE[res_id_tup][otherRes_id_tup][aa1][otherAA][candidate_rot][otherRot]
          witnessEnergyWithOtherRot = self.pairwiseE[res_id_tup][otherRes_id_tup][aa2][otherAA][witness_rot][otherRot]
          diff = candidateEnergyWithOtherRot - witnessEnergyWithOtherRot
          if diff < minimumForThisResidue: 
            e1_add = candidateEnergyWithOtherRot
            e2_add = witnessEnergyWithOtherRot
            minimumForThisResidue = diff
#      pruning_interval += minimumForThisResidue
      e1 += e1_add
      e2 += e2_add
    if e1 - e2 >= 0.0:
      return True
    else:
      return False

  # Call type dependent DEE on the matrix.
  def tdDEE(self, doPruning=True):
    self.isRotPruned = {}
    # Initialize the pruning matrix to false.
    print "Pruning search space using type dependent DEE "
    countTotalRots = 0
    
    for res_id_tup in self.getSortedResIds():
      self.isRotPruned[res_id_tup] = {}
      for aa in self.getAllowedAAs(res_id_tup):
        self.isRotPruned[res_id_tup][aa] = {}
        for rot in self.singleBody[res_id_tup][aa].keys():
          self.isRotPruned[res_id_tup][aa][rot] = False
          countTotalRots +=1 
    
    print "Total number of rotamers: "+`countTotalRots`

    rotamersPrunedThisRound = 1
    while doPruning and rotamersPrunedThisRound > 0:
      rotamersPrunedThisRound = 0
      # Candidate is the rotamer to be pruned; witness is the rotamer that will prune it.
      for res_id_tup in self.getSortedResIds():
        for aa in self.getAllowedAAs(res_id_tup):
          for candidate_rot in self.getRots(res_id_tup, aa):
              # In Type dependent DEE the candidate rotamer and the witness rotamer are of the same type. 
              for witness_rot in self.getRots(res_id_tup, aa):
                if witness_rot != candidate_rot and not self.isRotPruned[res_id_tup][aa][candidate_rot] \
                    and not self.isRotPruned[res_id_tup][aa][witness_rot]:
                  if self.can_prune_rot(res_id_tup, aa, candidate_rot, aa, witness_rot):
                    self.isRotPruned[res_id_tup][aa][candidate_rot] = True
                    #print "Pruned rotamer "+`res_id_tup`+" "+aa+" "+`candidate_rot`
                    rotamersPrunedThisRound += 1
      print "Rotamers pruned this round "+`rotamersPrunedThisRound`


