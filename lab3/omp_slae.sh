#!/bin/bash
#PBS -l walltime=00:30:00
#PBS -l select=1:ncpus=4:ompthreads=4
cd $PBS_O_WORKDIR


OMP_NUM_THREADS="1" ./omp_slae_schedule.out
OMP_NUM_THREADS="2" ./omp_slae_schedule.out
OMP_NUM_THREADS="4" ./omp_slae_schedule.out
OMP_NUM_THREADS="8" ./omp_slae_schedule.out
OMP_NUM_THREADS="16" ./omp_slae_schedule.out

