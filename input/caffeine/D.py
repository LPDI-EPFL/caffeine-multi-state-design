base_data_directory = "/Users/pablo/cur/caf/"
base_ematrix_directory = {}
base_ematrix_directory = base_data_directory+"02d-energyMatrices_groupD/"
base_structure_directory = base_data_directory+"01-inputStructs/"

stateDefs = []
stateDefs.append({'type':'positive','name':'02-synCafB:CFF:synCafH', 'id': 'posBH', 'color':'magenta'})
#stateDefs.append({'type':'positive','name':'08-synCafH:CFF', 'id': 'posH', 'color':'magenta'})
## Negative states:
# Homo: assign to those residues in the second chain the same amino acid as the corresponding in the first chain.
stateDefs.append({'type':'negative','name':'03-synCafB:CFF:synCafB', 'id': 'negBB', 'color':'orange', 'homo': ['B', 'H']})
stateDefs.append({'type':'negative','name':'04-synCafH:CFF:synCafH', 'id': 'negHH',  'color':'blue', 'homo': ['H', 'B']})
#stateDefs.append({'type':'negative','name':'06-synCafB:CFF', 'id': 'negB', 'color':'magenta'})
#stateDefs.append({'type':'negative','name':'08-synCafH:CFF', 'id': 'negH', 'color':'magenta'})
# Reference/constraint states
# Complex energy reference
stateDefs.append({'type':'reference','name':'01-wtCafB:CFF:wtCafH', 'id': 'complex', 'color':'magenta'})
stateDefs.append({'type':'constraint','name':'02-synCafB:CFF:synCafH', 'id': 'complex', 'color':'magenta'})
# B chain reference
#stateDefs.append({'type':'reference','name':'05-wtCafB:CFF', 'id': 'protein', 'color':'magenta'})
#stateDefs.append({'type':'constraint','name':'06-synCafB:CFF', 'id': 'protein', 'color':'magenta'})
# H chain reference
stateDefs.append({'type':'reference','name':'07-wtCafH:CFF', 'id': 'ligand', 'color':'magenta'})
stateDefs.append({'type':'constraint','name':'08-synCafH:CFF', 'id': 'ligand', 'color':'magenta'})
