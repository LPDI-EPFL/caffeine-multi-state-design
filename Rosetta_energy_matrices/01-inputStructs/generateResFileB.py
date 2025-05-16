#!/usr/local/bin/ipython
import os
defaultMuts = 'QNREDK'
#SynCafB = [('B', 47), ('B', 64) ]
SynCafB = [('B', 41)]
SynCafH = [('H', 64)]
SynCafH_B = [('B', 64)]
SynCafB_H = [('H', 47), ('H', 41)]
flexB = [('B', 47), ('B', 67)]
flexH = [('H', 67), ('H', 41)]
flex_wt_01 = [('B', 41), ('B', 47), ('B', 64), ('B', 67), ('H', 41), ('H', 47), ('H', 64), ('H', 67)]

mut_syn_02 = [('B', 47), ('H', 64)]
flex_syn_02 = [('B', 41), ('B', 64), ('B', 67), ('H', 41), ('H', 47), ('H', 67)]

mut_syn_03 = [('B', 47), ('H', 47)]
flex_syn_03 = [('B', 41), ('B', 64), ('B', 67), ('H', 41), ('H', 64), ('H', 67)]

mut_syn_04 = [('B', 64), ('H', 64)]
flex_syn_04 = [('B', 41), ('B', 47), ('B', 67), ('H', 41), ('H', 47), ('H', 67)]

flex_wt_05 = [('B', 41), ('B', 47), ('B', 64), ('B', 67)]

flex_syn_06 = [('B', 41), ('B', 64), ('B', 67)]
mut_syn_06 = [('B', 47)]

flex_wt_07 = [('H', 41), ('H', 47), ('H', 64), ('H', 67)]

flex_syn_08 = [('H', 41), ('H', 47), ('H', 67)]
mut_syn_08 = [('H', 64)]

states = {
    '01-wtCafB:CFF:wtCafH': {'flex': flex_wt_01, 'mut': []},
    '02-synCafB:CFF:synCafH': {'flex': flex_syn_02, 'mut': mut_syn_02},
    '03-synCafB:CFF:synCafB': {'flex': flex_syn_03,  'mut': mut_syn_03},
    '04-synCafH:CFF:synCafH': {'flex': flex_syn_04,  'mut': mut_syn_04},
    '05-wtCafB:CFF': {'flex': flex_wt_05,  'mut': []},
    '06-synCafB:CFF': {'flex': flex_syn_06,  'mut': mut_syn_06},
    '07-wtCafH:CFF': {'flex': flex_wt_07,  'mut': []},
    '08-synCafH:CFF': {'flex': flex_syn_08,  'mut': mut_syn_08},
}

header = 'NATRO\nstart\n'
for statekey in states: 

  if not os.path.exists(statekey):
    os.makedirs(statekey)
  resfile = open (statekey+'/'+'resfileB', 'w')
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



