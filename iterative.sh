#!/bin/bash

g++ -O3 -std=c++11 ../pin_buckling_new_3d/*.cpp

for niter in $(seq 1 500);do	
export niter
echo "submitting $niter"	
sbatch -o sim_${niter} -e err_${niter} --job-name sim_${niter} job.sbatch
sleep 1
done


