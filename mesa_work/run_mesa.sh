#!/bin/bash

date
hostname
# module load Anaconda3-5.1.0
source ~/.bash_profile

# cd "/scratch/ht06/yl7032/numax-sc-metallicity/hpc/coarse_v5/coarse$1/"
python3 "driver_prergbtip.py" $1

date
exit

###PBS -l nodes=node21:ppn=12
###PBS -l nodes=node43:ppn=12
###PBS -q physics
