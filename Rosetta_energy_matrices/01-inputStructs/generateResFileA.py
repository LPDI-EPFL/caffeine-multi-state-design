#!/usr/local/bin/ipython
import os
defaultMuts = 'QNREDKS'
SynCafB = [('B', 45)]
SynCafH = [('H', 45)]

states = {
    '01-wtCafB:CFF:wtCafH': {'flex': SynCafB+SynCafH, 'mut': []},
    '02-synCafB:CFF:synCafH': {'flex':[], 'mut': SynCafB+SynCafH},
    '03-synCafB:CFF:synCafB': {'flex':[],  'mut': SynCafB+SynCafH},
    '04-synCafH:CFF:synCafH': {'flex': [],  'mut': SynCafB+SynCafH},
    '05-wtCafB:CFF': {'flex': SynCafB,  'mut': []},
    '06-synCafB:CFF': {'flex': [],  'mut': SynCafB},
    '07-wtCafH:CFF': {'flex': SynCafH,  'mut': []},
    '08-synCafH:CFF': {'flex': [],  'mut': SynCafH},
}

header = 'NATRO\nstart\n'
for statekey in states: 

  if not os.path.exists(statekey):
    os.makedirs(statekey)
  resfile = open (statekey+'/'+'resfileA', 'w')
  resfile.write(header)
  for flexres in states[statekey]['flex']:
    chain = flexres[0]
    res = flexres[1]
    outline = `res`+' '+chain+' NATAA\n'
    resfile.write(outline)
  for mutres in states[statekey]['mut']:
    chain = mutres[0]
    res = mutres[1]
    if res == 113 and chain =='A':
      outline = `res`+' '+chain+' PIKAA '+'FYNQDE'+'\n'
    elif res == 59 and chain =='A':
      outline = `res`+' '+chain+' PIKAA '+"SNTV"+'\n'
    elif res == 11 and chain =='A':
      outline = `res`+' '+chain+' PIKAA '+"ASTN"+'\n'
    elif res == 14 and chain =='A':
      outline = `res`+' '+chain+' PIKAA '+"EQNDE"+'\n'
    else:
      outline = `res`+' '+chain+' PIKAA '+defaultMuts+'\n'
    
    resfile.write(outline)



