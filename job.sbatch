#!/bin/bash
#SBATCH --mail-user=emailaddress
#SBATCH -t 120:0:0
#SBATCH -N 1
#SBATCH -n 1

./a.out $niter