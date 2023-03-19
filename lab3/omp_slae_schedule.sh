#!/bin/bash
#PBS -l walltime=00:30:00
#PBS -l select=1:ncpus=4:ompthreads=4
cd $PBS_O_WORKDIR

echo "OMP_NUM_THREADS = $OMP_NUM_THREADS"

echo "Static"
echo

for (( i = 1; i <= 10; i++ ))
do
OMP_SCHEDULE="static, $i" ./omp_slae_schedule.out
done

for (( i = 50; i <= 100; i+=50))
do
OMP_SCHEDULE="static, $i" ./omp_slae_schedule.out
done

for (( i = 500; i <= 1500; i+=500 ))
do
OMP_SCHEDULE="static,$i" ./omp_slae_schedule.out
done
echo

echo "Dynamic"
echo

for (( i = 1; i <= 10; i++ ))
do
OMP_SCHEDULE="dynamic,$i" ./omp_slae_schedule.out
done

for (( i = 50; i <= 100; i+=50 ))
do
OMP_SCHEDULE="dynamic,$i" ./omp_slae_schedule.out
done

for (( i = 500; i <= 1500; i+=500 ))
do
OMP_SCHEDULE="dynamic,$i" ./omp_slae_schedule.out
done
echo

echo "Guided"
echo

for (( i = 1; i <= 10; i++ ))
do
OMP_SCHEDULE="guided,$i" ./omp_slae_schedule.out
done

for (( i = 50; i <= 100; i+=50 ))
do
OMP_SCHEDULE="guided,$i" ./omp_slae_schedule.out
done

for (( i = 500; i <= 1500; i+=500 ))
do
OMP_SCHEDULE="guided,$i" ./omp_slae_schedule.out
done
echo