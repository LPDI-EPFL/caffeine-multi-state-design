#!/bin/bash
for mydir in 0*
do
    cd $mydir
    ln -s ../CAF.pdb structure.pdb
    ln -s ../CFF.params .
    ln -s ../CFF_conformers.pdb .
    cd ..
done
