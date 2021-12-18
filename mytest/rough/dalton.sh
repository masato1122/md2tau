#!/bin/sh
#PBS -q default
#PBS -l nodes=1:ppn=2:cnode001
#PBS -j oe
#PBS -N test_rough

export OMP_NUM_THREADS=1
cd $PBS_O_WORKDIR
nprocs=2
#MPIRUN=/opt/intel/compilers_and_libraries_2017.4.196/linux/mpi/intel64/bin/mpirun

#----- RUN the JOB -----#
python ../../peak_detecter4SED.py -f SED_rtz1.txt -x 3 -y 4 --peakmin 10.0 \
     --calculation 1 -N $nprocs --nkskip 12 > log1.txt

for i in plot fitting; do
python ../../peak_analyzer.py --sfile SED_rtz1.txt -x 3 -y 4 --wunit Hz \
    --pfile SED_rtz1_PEAKS.txt --JTYPE $i
done

