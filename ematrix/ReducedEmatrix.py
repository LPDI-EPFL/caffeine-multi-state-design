class ReducedEMatrix:
  def __init__(self, numResidues, rotsPerRes, eIntraMat, ePairMat, map_rots_to_seq):
    self.numResidues = numResidues
    self.rotsPerRes = rotsPerRes
    self.eIntraMat = eIntraMat
    self.ePairMat = ePairMat
    self.map_rots_to_seq = map_rots_to_seq
