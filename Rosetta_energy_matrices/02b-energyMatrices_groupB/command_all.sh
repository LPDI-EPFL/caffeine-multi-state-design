#!/bin/bash
BASE_DIRECTORY="/Users/pablo/cur/caf/"
#for state in 01-wtCafD:CFF:wtCafG 02-synCafD:CFF:synCafG 03-synCafD:CFF:synCafD 04-synCafG:CFF:synCafG
for state in 01-wtCafB:CFF:wtCafH 02-synCafB:CFF:synCafH 03-synCafB:CFF:synCafB 04-synCafH:CFF:synCafH 05-wtCafB:CFF 06-synCafB:CFF 07-wtCafH:CFF 08-synCafH:CFF
do
  mkdir -p $state
  cd $state
  INPUT_STRUCT="$BASE_DIRECTORY/01-inputStructs/$state/structure.pdb"
  RESFILE="$BASE_DIRECTORY/01-inputStructs/$state/resfileB"
  rm structure_0001/*
  PARAMS_FILE="$BASE_DIRECTORY/01-inputStructs/$state/CFF.params"
  /Users/pablo/local_lpdi/main/source/bin/ig_dump.macosclangdebug -s $INPUT_STRUCT -resfile $RESFILE -use_input_sc -extra_res_fa $PARAMS_FILE 
  cd ..
done
