#!/usr/local/bin/ipython
from ematrix.EMatrix import *
from varbnb.VarbnbMSD import *
from config.config import *
import sys
import imp 


drug = sys.argv[1]
group = sys.argv[2]

# read data.
data = imp.load_source("data","input/"+drug+"/"+group+".py")

state={}
for statedef in data.stateDefs:
  mytype = statedef['type']
  myname = statedef['name']
  if mytype not in state:
    state[mytype] = {}
  state[mytype][myname] = {}
  state[mytype][myname]['pdbfile'] = data.base_structure_directory+myname+"/structure.pdb"
  state[mytype][myname]['color'] = statedef['color']
  state[mytype][myname]['ematrix_directory'] = data.base_ematrix_directory+myname+"/structure_0001/"
  print "Parsing "+state[mytype][myname]['ematrix_directory']
  state[mytype][myname]['ematrix'] = EMatrix()
  state[mytype][myname]['ematrix'].parseMatrix(state[mytype][myname]['ematrix_directory'])
  if 'id' in statedef:
    state[mytype][myname]['id'] = statedef['id']
  if 'homo' in statedef: 
    state[mytype][myname]['homo'] = statedef['homo'] 


if not os.path.exists("output/"+drug):
    os.makedirs("output/"+drug)

logfile= "output/"+drug+"/"+drug+"_"+group+".log"
print "Running, logging in "+logfile
varbnb_run = VarbnbMSD(state, logfile)
# Results are stored in state.
multistate_results = varbnb_run.optimize()
#multistate_results = varbnb_run.optimize({('B', 91): ['GLU']})
#multistate_results = varbnb_run.optimize({('B', 46): ['ARG']})
#multistate_results = varbnb_run.optimize({('B', 46): ['ARG'], ('H', 40): ['GLU'], ('H', 48): ['GLU']})
#multistate_results = varbnb_run.optimize({('B', 64): ['ARG'], ('B', 47): ['ARG'], ('H', 64): ['GLU'], ('H', 47): ['ASP']})


## Output results to output/ directory
import json 
## First convert tuples to string.

jsonarray = json.dumps(multistate_results, sort_keys=True, indent=4, separators=(',', ': '))
print jsonarray

outfile = open("output/"+drug+"/"+drug+"_"+group+".json", 'w')
outfile.write(jsonarray)
outfile.close()
