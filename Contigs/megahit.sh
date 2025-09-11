#!/bin/bash
#SBATCH --job-name=megahit

# Record the execution date, hostname, and current directory.
date
hostname
pwd

# Making sure logs directory exists.
mkdir -p logs

# Load the megahit module.
module load megahit/1.2.9

# Assign input parameters (assuming you will pass them as arguments or otherwise specify).
R1=$1
R2=$2

# Use basename to construct output prefix
output_prefix=$(basename "${R1}" .fwd_p.fq)

# Run megahit with specified parameters.
megahit -1 "${R1}" -2 "${R2}" -o "out_${output_prefix}"

# Print the completion date.
date
