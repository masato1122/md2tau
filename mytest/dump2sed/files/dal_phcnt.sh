#!/bin/sh
#PBS -q default
#PBS -l nodes=1:ppn=16:cnode011
#PBS -j oe
#PBS -N cnt_band

export OMP_NUM_THREADS=1
cd $PBS_O_WORKDIR
nprocs=16
MPIRUN=/opt/intel/compilers_and_libraries_2017.4.196/linux/mpi/intel64/bin/mpirun

#----- RUN the JOB -----#
$MPIRUN -np ${nprocs} lmp_intel < band_cnt.in

