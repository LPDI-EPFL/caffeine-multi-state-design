from heapq import *
import ipdb
import random
class SequenceNode:
  def __init__(self, parentPartialSequence={}, new_res_id_tup=None, newAA=None):
    self.isExactScore = False
    self.fScore = float("-inf")
    # AvailableRots is a list of lists saying which rotamers are available for this node.
    self.partialSequence = dict(parentPartialSequence)
    self.myResIdTup = new_res_id_tup
    self.myAAtype = newAA
    if self.myAAtype != None:
      self.partialSequence[self.myResIdTup] = [self.myAAtype]
    self.next_res_id_tup = ""
    self.isExactScore = False

  def printPartialSequence(self):
    for res_id_tup in self.partialSequence[self.myResIdTup].keys():
      print `res_id_tup`+self.partialSequence[res_id_tup]

  # Returns a partial sequence where the amino acids of chain 
  # with label chainA are also assigned 
  # to chain with label chainB
  def getHomoDimericSequence(self, chainA, chainB):
    new_part_seq = {}
    for key in self.partialSequence:
      chain = key[0]
      resi = key[1]
      aatype = self.partialSequence[key]
      if chain == chainA:
        new_part_seq[(chainA, resi)] = aatype
        new_part_seq[(chainB, resi)] = aatype
    return new_part_seq

