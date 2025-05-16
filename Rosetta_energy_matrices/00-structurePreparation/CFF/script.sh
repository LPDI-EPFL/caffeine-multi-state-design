#babel -h CFF.sdf CFF_withH.sdf
#ROSETTA_BIN=/home/gainza/local_lpdi/Rosetta/main/source/scripts/python/public/molfile_to_params.py
MOLFILE_TO_PARAMS=/home/gainza/local_lpdi/Rosetta/main/source/scripts/python/public/molfile_to_params.py
$MOLFILE_TO_PARAMS -n CFF -p CFF --conformers-in-one-file CFF_withH.sdf
