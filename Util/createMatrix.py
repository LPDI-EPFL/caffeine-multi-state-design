# Creates a 4d rotamer matrix of the form [resI][rotR][resJ][rotS]
# Initialize value could be:
#       "None" if the matrix has to be initialized to None
#       "0.0" for a float initialization
#       "0" for an integer
from collections import defaultdict
from collections import Counter
def create4DRotMatrix(numRes, rotsPerRes, initializeValue):
  rot4Dmat = [None]*numRes
  for resI in range(numRes):
    rot4Dmat[resI] = [None]*rotsPerRes[resI]
    for rotR in range(rotsPerRes[resI]):
      rot4Dmat[resI][rotR] = [None]*numRes
      for resJ in range(numRes):
        rot4Dmat[resI][rotR][resJ] = [None]*rotsPerRes[resJ] 
        for rotS in range(rotsPerRes[resJ]):
          rot4Dmat[resI][rotR][resJ][rotS] = initializeValue
  return rot4Dmat

# Creates a 3d rotamer matrix of the form [senderRes][recvRes][recvRot]
# Initialize value could be:
#       "None" if the matrix has to be initialized to None
#       "0.0" for a float initialization
#       "0" for an integer
def create3DMsgMat(numRes, rotsPerRes, initializeValue):
  msg3Dmat = [None]*numRes
  for sendRes in range(numRes):
    msg3Dmat[sendRes] = [None]*numRes
    for recvRes in range(numRes):
      msg3Dmat[sendRes][recvRes] = [None]*rotsPerRes[recvRes] 
      for recvRot in range(rotsPerRes[recvRes]):
        msg3Dmat[sendRes][recvRes][recvRot] = initializeValue
  return msg3Dmat

# Creates a 2d rotamer matrix of the form [resI][rotR]
# Initialize value could be:
#       "None" if the matrix has to be initialized to None
#       "0.0" for a float initialization
#       "0" for an integer
def create2DRotMatrix(numRes, rotsPerRes, initializeValue):
  rot2Dmat = [None]*numRes
  for resI in range(numRes):
    rot2Dmat[resI] = [None]*rotsPerRes[resI]
    for rotR in range(rotsPerRes[resI]):
      rot2Dmat[resI][rotR] = initializeValue
  return rot2Dmat

# Creates a 2d residue interaction matrix of the form [resI][resJ]
# Initialize value could be:
#       "None" if the matrix has to be initialized to None
#       "0.0" for a float initialization
#       "0" for an integer
def create2DResMatrix(numRes, initializeValue):
  res2Dmat = [None]*numRes
  for resI in range(numRes):
    res2Dmat[resI] = [None]*numRes
    for resJ in range(numRes):
      res2Dmat[resI][resJ] = initializeValue
  return res2Dmat

def multi_dimensions(n, type):
  """ Creates an n-dimension dictionary where the n-th dimension is of type 'type'
  """
  if n<=1:
    return type(0.0)
  return defaultdict(lambda:multi_dimensions(n-1, type))
