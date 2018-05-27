#!/bin/bash
species=32

for ((mastpos=1;mastpos<=$[$species];mastpos++))
do

export t=$mastpos;
sbatch MASTIF.pbs 
sleep 60.0s

done
