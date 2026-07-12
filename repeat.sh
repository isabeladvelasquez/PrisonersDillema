#!/bin/bash

reps=100
folder1="SPD_2C_L100_T10000_TIMEvsPCP"
folder2="SPD_2C_L100_T10000_2CPROB"
folder3="SPD_2C_L100_T10000_SNAPSHOT"

gcc '2C_L100_RHOvar_Mvar.c' lat2eps_lib.c -o '2C_L100_RHOvar_Mvar' -O3 -lgsl -lgslcblas -lm

for i in $(seq 1 $reps); do
    ./2C_L100_RHOvar_Mvar 
done

#mkdir -p "$folder1"
#mkdir -p "$folder2"
#mkdir -p "$folder3"

#mv *SAVETIME*.dat "$folder1/"/ 2>/dev/null
#mv *SAVE2CPROB*.dat "$folder2/"/ 2>/dev/null
#mv *ATVS*.eps "$folder3/"/ 2>/dev/null
