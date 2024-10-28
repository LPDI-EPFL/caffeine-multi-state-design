from PGQueueNode import *
import os
from PGExpansionQueue import *
from sys import *
from mplp.MPLP import *
from ematrix.EMatrix import *
import config
import time
import pdb
class DynamicAS:

  def __init__(self, ematrix, assignedAAs={}):
    # Log levels: 0, superSmall, 1 small, 2 full
    self.logLevel = 0
    #self.logFile = open("output/"+str(os.getpid()) + ".out", "w", 256)
    self.ematrix = ematrix
    self.allowedAAs = {}
    ## Compute the amino acids that will be allowed during the search.
    for res_id_tup in ematrix.getSortedResIds():
      self.allowedAAs[res_id_tup] = {}
      # If there is only one amino acid allowed for this position in the ematrix
      if len(ematrix.getAllowedAAs(res_id_tup)) == 1:
        # Then assignedAAs is ignored. 
        self.allowedAAs[res_id_tup] = ematrix.getAllowedAAs(res_id_tup)
      elif res_id_tup in assignedAAs:  # AA has been assigned.
        self.allowedAAs[res_id_tup] = assignedAAs[res_id_tup]
      else: # Otherwise assign all others.
        self.allowedAAs[res_id_tup] = ematrix.getAllowedAAs(res_id_tup)

    self.relaxationMaster = MPLP(ematrix)
    print self.allowedAAs

  def doAS(self):
    curExpansion = PGExpansionQueue()
    rootNode = PGQueueNode(self.allowedAAs)
    rootNode.fScore = (self.relaxationMaster.optimizeEMPLP(self.allowedAAs))[0]
    rootNode.my_next_res_id_tup = self.relaxationMaster.getHeuristicNextResidueToExpandBasedOnBeliefs()
    curExpansion.insert(rootNode)

    while(not curExpansion.is_empty()):
      expNode = curExpansion.pop()[1]
      self.log(1,"Expanding node, res "+`expNode.myResId`+" aa "+`expNode.myAA`+" rot " +`expNode.myRot`+" fScore "+`expNode.fScore`)
      self.log(1,"Assigned sequence space: "+ `expNode.myConf`)

      if(len(expNode.myConf.keys()) == (len(self.ematrix.getSortedResIds()))):
        # We are done!
        self.log(1,"Done!! energy: "+`expNode.fScore`)
        #self.logFile.close()
        return expNode
      # Expand next residue
      next_res_id_tup = expNode.my_next_res_id_tup
      for aa in expNode.myAllowedAAs[next_res_id_tup]:
        for rot in self.ematrix.getRots(next_res_id_tup, aa):
          newNode = PGQueueNode(expNode.myAllowedAAs, expNode.myConf, next_res_id_tup, aa, rot)
          newNode.fScore = (self.relaxationMaster.optimizeEMPLP(newNode.myAllowedAAs, newNode.myConf))[0]
          newNode.my_next_res_id_tup = self.relaxationMaster.getHeuristicNextResidueToExpandBasedOnBeliefs(newNode.myConf)
          curExpansion.insert(newNode)
    print "Error with dynamic A*"
    return None

  def log(self, loglevel, string):
    if loglevel <= self.logLevel:
      self.logFile.write(string+"\n")
      print string
