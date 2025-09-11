#!/bin/bash
#SBATCH --job-name=IQTREE
#SBATCH --mem=20G
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1

module load iq-tree/2.2.2.7 

RG=$1

iqtree -s ${RG} -m MFP -bb 2000 -nt AUTO -lmap 20000

