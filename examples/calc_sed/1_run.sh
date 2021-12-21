#!/bin/sh
WORKDIR=$(cd $(dirname $0); pwd)

cd $WORKDIR
BIN=../../dump2sed/dump2sed
mpirun -np 1 $BIN sed.in

#python ../../Util/sed_sum_xyz.py -1 SED_X.txt -2 SED_Y.txt -3 SED_Z.txt

#for i in 0; do
#python ${SEDROOT}/dump2band_cylind_rotate.py --file md4band.dump --pfile data.lammps \
#    --nuat 40 --nrot 10 --Tcoord $i \
#    --CALTYPE serial
#done

#sh ${SEDROOT}/draw_sed.sh 1
#for itype in 1 2 3; do
#    sh ${SEDROOT}/draw_sed_detail.sh 0 $itype
#done

