#!/usr/local/bin/ipython
from rosetta import *
rosetta.init()
pdb_structure_negative = "/Users/pablo/cur/msd_data/01-inputStructs/bclxl_bh3/3fdl_bclxl_bh3_min_res-renum.pdb"

pose=pose_from_pdb(pdb_structure_negative)
design

