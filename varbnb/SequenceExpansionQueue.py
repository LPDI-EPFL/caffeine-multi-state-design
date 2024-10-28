from heapq import *
from varbnb.SequenceNode import * 
# A priority queue for SequenceNode
class SequenceExpansionQueue:
  def __init__(self):
    self.thequeue = []
  def insert(self, node):
    mytuple = (node.fScore, node)
    heappush(self.thequeue,mytuple)
  def pop(self):
    return heappop(self.thequeue)
  def is_empty(self):
    if len(self.thequeue) > 0:
      return False
    else:
      return True
