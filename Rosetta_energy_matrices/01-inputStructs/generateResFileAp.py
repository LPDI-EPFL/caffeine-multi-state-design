#!/usr/local/bin/ipython
import os
defaultMuts = 'QNREDKS'
SynCafB = [('B', 45), ('B', 91) ]
SynCafH = [('H', 45)]
SynCafB_H = [('H', 91)]
#SynCafH = [('H', 45), ('H', 91) ]

states = {
    #'01-wtCafB:CFF:wtCafH': {'flex': SynCafB+SynCafH, 'mut': []},
    '01-wtCafB:CFF:wtCafH': {'flex': SynCafB_H+SynCafB+SynCafH, 'mut': []},
    '02-synCafB:CFF:synCafH': {'flex':SynCafB_H, 'mut': SynCafB+SynCafH},
    '03-synCafB:CFF:synCafB': {'flex':[],  'mut': [('B', 45), ('B', 91), ('H', 91), ('H', 45)]},
    '04-synCafH:CFF:synCafH': {'flex': [('H', 91), ('B', 91)],  'mut': [('H', 45), ('B', 45)]},
    '05-wtCafB:CFF': {'flex': SynCafB,  'mut': []},
    '06-synCafB:CFF': {'flex': [],  'mut': SynCafB},
    '07-wtCafH:CFF': {'flex': SynCafH+SynCafB_H,  'mut': []},
    '08-synCafH:CFF': {'flex': [],  'mut': SynCafH+SynCafB_H},
}

header = 'NATRO\nstart\n'
for statekey in states: 

  if not os.path.exists(statekey):
    os.makedirs(statekey)
  resfile = open (statekey+'/'+'resfileAp', 'w')
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



