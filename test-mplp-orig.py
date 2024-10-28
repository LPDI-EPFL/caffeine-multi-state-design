#!/usr/local/bin/ipython
from util.createMatrix import *
from mplp.MPLPorig import *

numRes = 5
rotsPerRes = [1]*5
intraEmat = create2DRotMatrix(numRes, rotsPerRes, 0.0)
intraEmat[0][0] = 0.554935
intraEmat[1][0] = 0.577375
intraEmat[2][0] = -6.18313
intraEmat[3][0] = -4.17289
intraEmat[4][0] = -3.58061
pairEmat = create4DRotMatrix(numRes, rotsPerRes, 0.0)
pairEmat[0][0][1][0] = 0.880953
pairEmat[0][0][2][0] = 0.0396635
pairEmat[0][0][3][0] = 0.0772491
pairEmat[1][0][0][0] = 0.880953
pairEmat[1][0][2][0] = 0.196624
pairEmat[1][0][3][0] = 0.088757
pairEmat[1][0][4][0] = 0.0
pairEmat[2][0][0][0] = 0.0396635
pairEmat[2][0][1][0] = 0.196624
pairEmat[2][0][3][0] = -0.198208
pairEmat[2][0][4][0] = -4.1227
pairEmat[3][0][0][0] = 0.0772491
pairEmat[3][0][1][0] = 0.088757
pairEmat[3][0][2][0] = -0.198208
pairEmat[3][0][4][0] = -1.08624
pairEmat[4][0][1][0] = 0.0
pairEmat[4][0][2][0] = -4.1227
pairEmat[4][0][3][0] = -1.08624

mplporig = MPLPorig(numRes, rotsPerRes, intraEmat, pairEmat)
availableRots =[]
for i in range(numRes):
	availableRots.append([0])
ebound , sol = mplporig.optimizeEMPLP(availableRots, 1)
print "ebound ="+`ebound`

