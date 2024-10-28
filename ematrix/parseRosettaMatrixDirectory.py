from os import listdir
from os.path import isfile, join
from config import config
import sys
import gzip
import pdb
def parseRosettaMatrixDirectory(matrix_directory):
  onlyfiles = [f for f in listdir(matrix_directory) if isfile(join(matrix_directory, f))]
  ## We first read all the single body energy terms because these contain rotamer and aa type information. 
  singleBody = {}
  for matrix_filename in onlyfiles:
    if '-' not in matrix_filename: 
      ## Residue id and chain number are only in the filename.
      ## TODO: For now the sequence number has to be an integer. If it has a kabat numbering it is not supported.
      try:
        res_seq_num = int(matrix_filename.split("_")[0])
      except:
        print "Error: non-numeric residue numbering is not supported."
        sys.exit(1)

      chain = matrix_filename.split("_")[1].split(".")[0].rstrip()
      res_id_tup = (chain, res_seq_num)
      ## Open file as gzip file. 
      f=gzip.open(matrix_directory+"/"+matrix_filename,'rb')
      ## Second dimension is the AA id.
      singleBody[res_id_tup] = {}
      singleBody[res_id_tup]['_order'] = [] # It is important to maintain track on the order in which AA are ordered. 
      aa_count = 0
      rot_counter = 0
      for line in f:
        fields = line.split()
        energy = float (fields[0])
        if config.USE_STERIC_E:
          if energy > config.STERIC_E:
            energy = float("inf")
        aatype = fields[1]
        dihedrals = [float (x) for x in fields [2:] ]
        if aatype not in singleBody[res_id_tup]:
          singleBody[res_id_tup][aatype] = {}
          singleBody[res_id_tup]['_order'].append(aatype)
          rot_counter = 0 
        singleBody[res_id_tup][aatype][rot_counter] = ({'E': energy, 'Dih': dihedrals})
        rot_counter += 1

  # Create a map between the global number of rotamers at a residue and the rotamer index for each amino acid. 
  for res_id_tup in singleBody.keys():
    singleBody[res_id_tup]['_map_res_rot_to_aa_rot'] = []
    for aa in singleBody[res_id_tup]['_order']:
      for rot_for_aa in range(len(singleBody[res_id_tup][aa])):
        singleBody[res_id_tup]['_map_res_rot_to_aa_rot'].append((aa, rot_for_aa))

  ## We now fill all the pairwise energies. 
  pairwiseE = {}
  for matrix_filename in onlyfiles:
    if '-' in matrix_filename: 
      ## Residue id and chain number are only in the filename.
      res_seq_num1 = int(matrix_filename.split("-")[0].split("_")[0])
      chain1 = matrix_filename.split("-")[0].split("_")[1]
      res_id_tup1 = (chain1, res_seq_num1)
      res_seq_num2 = int(matrix_filename.split("-")[1].split("_")[0])
      chain2 = matrix_filename.split("-")[1].split("_")[1].split(".")[0]
      res_id_tup2 = (chain2, res_seq_num2)
      ## Initialize the pairwise energy matrices
      if res_id_tup1 not in pairwiseE: 
        pairwiseE [res_id_tup1] = {}
      pairwiseE[res_id_tup1][res_id_tup2] = {}
      for aa1 in singleBody[res_id_tup1]['_order']:
        pairwiseE[res_id_tup1][res_id_tup2][aa1] = {}
        for aa2 in singleBody[res_id_tup2]['_order']:
          pairwiseE[res_id_tup1][res_id_tup2][aa1][aa2] = {} # [[] for x in range(len(singleBody[res_id_tup1][aa1]))]
          for rot1_ix in singleBody[res_id_tup1][aa1].keys():
            pairwiseE[res_id_tup1][res_id_tup2][aa1][aa2][rot1_ix] = {} # [0.0]*len(singleBody[res_id_tup2][aa2])
      # Matrices are symmetrical. Do everything again but in reverse.
      #  NOTE: Technically one could use a triangular matrix. However I prefer to have them symmetrical even if it takes twice as much space.
      if res_id_tup2 not in pairwiseE:
        pairwiseE[res_id_tup2] = {}
      pairwiseE[res_id_tup2][res_id_tup1] = {}
      for aa2 in singleBody[res_id_tup2]['_order']:
        pairwiseE[res_id_tup2][res_id_tup1][aa2] = {}
        for aa1 in singleBody[res_id_tup1]['_order']:
          pairwiseE[res_id_tup2][res_id_tup1][aa2][aa1] = {}  #[[] for x in range(len(singleBody[res_id_tup2][aa2]))]
          for rot2_ix in singleBody[res_id_tup2][aa2].keys():
            pairwiseE[res_id_tup2][res_id_tup1][aa2][aa1][rot2_ix] = {}  #[0.0]*len(singleBody[res_id_tup1][aa1])

      # Open the pairwise file as gzip file
      f=gzip.open(matrix_directory+"/"+matrix_filename,'rb')
      # Now read the actual energies
      line_counter = 0
      for line in f:
        energies = [float(x) for x in line.split()]
        aa1, rot_ix_for_aa1 = singleBody[res_id_tup1]['_map_res_rot_to_aa_rot'][line_counter]
        for pairE_ix in range(len(energies)):
          try:
            aa2, rot_ix_for_aa2 = singleBody[res_id_tup2]['_map_res_rot_to_aa_rot'][pairE_ix]
            energy = energies[pairE_ix]
            if config.USE_STERIC_E:
              if energy > config.STERIC_E:
                energy = float("inf")
            pairwiseE[res_id_tup1][res_id_tup2][aa1][aa2][rot_ix_for_aa1][rot_ix_for_aa2] = energy
            pairwiseE[res_id_tup2][res_id_tup1][aa2][aa1][rot_ix_for_aa2][rot_ix_for_aa1] = energy
          except:
            pdb.set_trace()
        line_counter += 1
      
  return singleBody, pairwiseE
  
