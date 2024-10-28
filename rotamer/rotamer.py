# Atoms for each side-chain angle for each residue
import sys
from pymol import cmd 
from pymol import editing
CHIS = {}
CHIS["ARG"] = [ ["N","CA","CB","CG" ],
                ["CA","CB","CG","CD" ],
                ["CB","CG","CD","NE" ],
                ["CG","CD","NE","CZ" ]
              ]
 
CHIS["ASN"] = [ ["N","CA","CB","CG" ],
                ["CA","CB","CG","OD1" ]
              ]
 
CHIS["ASP"] = [ ["N","CA","CB","CG" ],
                ["CA","CB","CG","OD1" ]
              ]
CHIS["CYS"] = [ ["N","CA","CB","SG" ]
              ]
 
CHIS["GLN"] = [ ["N","CA","CB","CG" ],
                ["CA","CB","CG","CD" ],
                ["CB","CG","CD","OE1"]
              ]
 
CHIS["GLU"] = [ ["N","CA","CB","CG" ],
                ["CA","CB","CG","CD" ],
                ["CB","CG","CD","OE1"]
              ]
 
CHIS["HIS"] = [ ["N","CA","CB","CG" ],
                ["CA","CB","CG","ND1"]
              ]
 
CHIS["ILE"] = [ ["N","CA","CB","CG1" ],
                ["CA","CB","CG1","CD1" ]
              ]
 
CHIS["LEU"] = [ ["N","CA","CB","CG" ],
                ["CA","CB","CG","CD1" ]
              ]
 
CHIS["LYS"] = [ ["N","CA","CB","CG" ],
                ["CA","CB","CG","CD" ],
                ["CB","CG","CD","CE"],
                ["CG","CD","CE","NZ"]
              ]
 
CHIS["MET"] = [ ["N","CA","CB","CG" ],
                ["CA","CB","CG","SD" ],
                ["CB","CG","SD","CE"]
              ]
 
CHIS["PHE"] = [ ["N","CA","CB","CG" ],
                ["CA","CB","CG","CD1" ]
              ]
 
CHIS["PRO"] = [ ["N","CA","CB","CG" ],
                ["CA","CB","CG","CD" ]
              ]
 
CHIS["SER"] = [ ["N","CA","CB","OG" ]
              ]
 
CHIS["THR"] = [ ["N","CA","CB","OG1" ]
              ]
 
CHIS["TRP"] = [ ["N","CA","CB","CG" ],
                ["CA","CB","CG","CD1"]
              ]
 
CHIS["TYR"] = [ ["N","CA","CB","CG" ],
                ["CA","CB","CG","CD1" ]
              ]
 
CHIS["VAL"] = [ ["N","CA","CB","CG1" ]
              ]
 
# Set a rotamer, based on a selection, a restype and chi angles
# chi is a list of angles.
def set_rotamer(sel, chi):
  at = cmd.get_model("byres ("+sel+")").atom[0]
 
  print "My residue name: "+at.resn

  list = chi
  for i in range(len(CHIS[at.resn])):
      print "Setting Chi"+str(i+1)+" to "+str(list[i])
      print(sel + ' and name '+CHIS[at.resn][i][0])
      editing.set_dihedral(sel + ' and name '+CHIS[at.resn][i][0],
                         sel + ' and name '+CHIS[at.resn][i][1],
                         sel + ' and name '+CHIS[at.resn][i][2],
                         sel + ' and name '+CHIS[at.resn][i][3], str(list[i]))
    # Remove some objects that got created
  cmd.delete("pk1")
  cmd.delete("pk2")
  cmd.delete("pkmol")
