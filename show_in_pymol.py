#!/usr/local/bin/python
### NOTE: Replace this with the location of your pymol installation.
import sys
#sys.path.append("/usr/local/Cellar/pymol/1.7.6.0_1/libexec/lib/python2.7/site-packages/")
#sys.path.append("/Applications/PyMOL.app/Contents/lib/python2.7/site-packages")
import pymol
from pymol import cmd
from rotamer.rotamer import *
import json

# Load json file 
jsonfile = open(sys.argv[1])
state = json.load(jsonfile)
pymol.finish_launching()
#for state in topNodes.keys():
ix = 0

name_dict = {}
name_dict['02-synCafB:CFF:synCafH'] = 'CIH1_CIH2'
name_dict['03-synCafB:CFF:synCafB'] = 'CIH1_CIH1'
name_dict['04-synCafH:CFF:synCafH'] = 'CIH2_CIH2'
name_dict['06-synCafB:CFF'] = 'CIH1'
name_dict['08-synCafH:CFF'] = 'CIH2'
chain_color = {}
chain_color['02-synCafB:CFF:synCafH'] = {'B': 'paleyellow', 'H': 'lightblue'}
chain_color['03-synCafB:CFF:synCafB'] = {'B': 'paleyellow', 'H': 'paleyellow'}
chain_color['04-synCafH:CFF:synCafH'] = {'B': 'lightblue', 'H': 'lightblue'}
chain_color['06-synCafB:CFF'] = {'B': 'paleyellow', 'H': 'paleyellow'}
chain_color['08-synCafH:CFF'] = {'B': 'lightblue', 'H': 'lightblue'}
for mytype in state: # 'positive' or 'negative'
    if isinstance(state[mytype], dict):
        for mystatename in state[mytype]:
            ix = ix+1
            objname = name_dict[mystatename]
            #objname = mytype+"_"+mystatename
#            mycolor = state[mytype][mystatename]['color']
            print(state[mytype][mystatename]['pdbfile'])

            cmd.load(state[mytype][mystatename]['pdbfile'], objname)
            cmd.wizard("mutagenesis")
            cmd.show("cartoon", objname)
            for chain in chain_color[mystatename].keys():
                cmd.color(chain_color[mystatename][chain], objname+" and elem c and chain "+chain)
            cmd.hide("lines", objname+" and elem h")

            for res_id_tup in state[mytype][mystatename]['residues'].keys():
                aa = state[mytype][mystatename]['residues'][res_id_tup]['aa']
                print "Mutating residue "+res_id_tup+" to "+aa
                chain=res_id_tup[2]
                res=res_id_tup[5:-1]
                selection=objname+" and chain "+chain+" and resi "+res
                print (selection)
                pymol.cmd.select("s_"+objname, selection)
                #selection2 = "(s_"+objname+")"
                selection2 = "s_"+objname
                pymol.cmd.get_wizard().set_mode(aa)
                pymol.cmd.get_wizard().do_select(selection2)
                pymol.cmd.get_wizard().apply()
                dih = state[mytype][mystatename]['residues'][res_id_tup]['Dih']
                selection=objname+" and chain "+chain+" and resi "+res
                pymol.cmd.select("s_"+objname, selection)
                if aa != 'ALA':
                    print('Setting dihedrals {} for object {}'.format("s_"+objname, dih))
                    set_rotamer("s_"+objname, dih)
                pymol.cmd.show("sticks", selection)
                pymol.cmd.delete("s_"+objname)

# Draw wildtype 
ix = ix+1
mytype = 'positive'
mystatename = '02-synCafB:CFF:synCafH'
objname = 'wtVH_wtVH'

cmd.load(state[mytype][mystatename]['pdbfile'], objname)
cmd.show("cartoon", objname)
for chain in chain_color[mystatename].keys():
    cmd.color('palegreen', objname+" and elem c and chain "+chain)
cmd.hide("lines", objname+" and elem h")

for res_id_tup in state[mytype][mystatename]['residues'].keys():
    aa = state[mytype][mystatename]['residues'][res_id_tup]['aa']
    print "Mutating residue "+res_id_tup+" to "+aa
    chain=res_id_tup[2]
    res=res_id_tup[5:-1]
    selection=objname+" and chain "+chain+" and resi "+res
    #selection2 = "(s_"+objname+")"
    selection=objname+" and chain "+chain+" and resi "+res
    pymol.cmd.select("s_"+objname, selection)
    pymol.cmd.show("sticks", selection)
    pymol.cmd.delete("s_"+objname)



