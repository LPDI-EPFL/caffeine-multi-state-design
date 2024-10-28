import pdb
import math

SMALL_FLOAT = 10e-8

# Returns the root mean squares of a list of scores 
def rootMeanSquares(scores):
  denom = (float(len(scores)))
  mysum = 0.0
  numerator = mysum + (x*x for x in scores)
  return math.sqrt(numerator/denom)

def harmonicMean(scores):
  sumInverted = 0.0
  for i in range(len(scores)):
    if scores[i] == 0.0:
      return 0.0
    sumInverted += 1.0/scores[i]
  return len(scores)/sumInverted

def eq( a, b, eps=0.0001 ):
  # Division by zero
  if max(abs(a), abs(b)) <= eps:
    return True
  else:
    return abs( a - b ) / max(abs(a),abs(b)) <= eps

def listContains(myList, myIntElement):
  for elem in range(len(myList)):
    if elem == myIntElement:
      return True
  return False

# Tell me the index of element with the lowest value in the array (for positive
# floats only, for now)
def argMin(myList):
  minIndex = -1
  minValue = float('inf')
  for index in range(len(myList)):
    if myList[index] < minValue:
      minIndex = index
      minValue = myList[index]
  return minIndex

# This harmonic mean ignores the lowest element in the array.
def harmonicMean2ndBest(scores):
  minElemIndex = argMin(scores)
  sumInverted = 0.0
  for i in range(len(scores)):
    if i != minElemIndex:
      if scores[i] == 0.0:
        return 0.0
      sumInverted += 1.0/scores[i]
  if sumInverted == 0.0:
    return 0.0
  return (len(scores)-1)/sumInverted

# returns the index of the minimum element of an array
def argmin(array): 
  minIndex = -1
  minValue = float("inf")
  for elemIndex in range(len(array)):
    if array[elemIndex] < minValue:
      minValue = array[elemIndex]
      minIndex = elemIndex
  return minIndex
