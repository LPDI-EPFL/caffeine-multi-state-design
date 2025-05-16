#!/usr/local/bin/ipython
import os
defaultMuts = 'FWYLIAVHGDQNSRKE'
#SynCafB = [('B', 47), ('B', 64) ]
SynCafB = [('B', 49), ('B', 108)]
SynCafH = [('H', 49), ('H', 108)]
SynCafH_B = []
SynCafB_H = []
flexB = []
flexH = []

states = {
    '01-wtCafB:CFF:wtCafH': {'flex': SynCafB+SynCafH+flexB+flexH, 'mut': []},
    '02-synCafB:CFF:synCafH': {'flex':flexB+flexH+SynCafH_B+SynCafB_H, 'mut': SynCafB+SynCafH},
    '03-synCafB:CFF:synCafB': {'flex':flexB+flexH+SynCafH+SynCafH_B,  'mut': SynCafB+SynCafB_H},
    '04-synCafH:CFF:synCafH': {'flex': flexB+flexH+SynCafB+SynCafB_H,  'mut': SynCafH_B+SynCafH},
    '05-wtCafB:CFF': {'flex': SynCafB+flexB+SynCafH_B,  'mut': []},
    '06-synCafB:CFF': {'flex': flexB,  'mut': SynCafB+SynCafH_B},
    '07-wtCafH:CFF': {'flex': SynCafH+flexH+SynCafB_H,  'mut': []},
    '08-synCafH:CFF': {'flex': flexH,  'mut': SynCafH+SynCafB_H},
}

header = 'NATRO\nstart\n'
for statekey in states: 

  if not os.path.exists(statekey):
    os.makedirs(statekey)
  resfile = open (statekey+'/'+'resfileC', 'w')
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



