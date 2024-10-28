# When we read in the energy matrix, anything below lowerCutoff is set to zero
lowerCutoff = 0.001
# Similarly, any energy above 1000, is set to 1000.
upperCutoff = 1000
# If MAX_PRODUCT is set to true, max product belief propagation is used, otherwise sum-product is used
MAX_PRODUCT = False
MPLP_ITERATIONS = 3
BP_ITERATIONS = 3
# Set any energy above a threshold to float("inf")
USE_STERIC_E = False
STERIC_E = 10
# Do Type dependent DEE when reading an eenergy matrix.
DO_TD_DEE = True
# Constraint thresholds for constraints 
CONSTRAINT_THRESHOLD_E = 5.0
CONSTRAINT_THRESHOLD_DELTA_E = 5.0
# Signs for constraint Delta E
SIGN = { 'complex': 1.0, 'ligand':-1.0, 'protein':-1.0, 'nodeltaE': 0.0}
