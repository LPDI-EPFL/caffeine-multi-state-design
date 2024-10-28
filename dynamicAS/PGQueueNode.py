from heapq import *
import pdb
import random
def CompareNodes(nodeA, nodeB): 
  if nodeA.fScore > nodeB.fScore:
    return 1
  elif nodeA.fScore == nodeB.fScore:
    return 0
  else:
    return -1
class PGQueueNode:
  my_next_res_id_tup = ""
  fScore = float("-inf")
  def __init__(self, parentAllowedAAs, parentConf={}, myResId="", myAA="", myRot=""):
    # We assign the new rotamers at this position.  
    self.myResId = myResId
    self.myRot = myRot
    self.myAA = myAA
    self.myConf = dict(parentConf)
    self.myAllowedAAs = dict(parentAllowedAAs)
    if myResId != "":
      self.myConf[myResId] = {}
      self.myConf[myResId][myAA] = myRot
      self.myAllowedAAs[myResId] = [myAA]
