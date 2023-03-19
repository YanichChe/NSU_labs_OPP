#!/bin/sh
#PBS -l walltime=00:00:50
#PBS -l select=2:ncpus=12:ompthreads=24

cd $PBS_O_WORKDIR
echo "OMP_NUM_THREADS = $OMP_NUM_THREADS"
echo
./omp_slae.out