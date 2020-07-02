#!/bin/bash

#export EIC_LEVEL=dev
#source /afs/rhic.bnl.gov/eic/restructured/etc/eic_cshrc.csh

# Other arguments for script to run
breit=0


for i in {0..99}
    do
    ./EventLoop ../MCData/100Minclusive/truth_pE275_pE18_minqsq9_$i.root ../MCData/100Minclusive/smeared_pE275_pE18_minqsq9_$i.root dataFiles/100Minclusive/pE275_pE18_minqsq9_frame${breit}_d0kpi_$i.root $breit
    done


