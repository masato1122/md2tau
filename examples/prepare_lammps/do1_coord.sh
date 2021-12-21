#!/bin/sh
WORKDIR=$(cd $(dirname $0); pwd)
#UDIR=/home/ohnishi/work/Peapod/4_Again/UTIL
UDIR=./UTIL

if [ $# -ge 1 ]; then
    ical=$1
else
    ical=new
fi

ODIR=./CAL${ical}
if [ -e $ODIR ]; then
    echo " $ODIR already exists."
    exit
fi

# ---- parameters ----------------------------
# Number of unit cells of the SWNT
Ncell=32
#Ncell=1024
Next=0        # adjusting # of C60s

# --- for test
#tnpt=0; tpre=10     # ps
# --- for accurate calculation
tnpt=0; tpre=500     # ps

nstructures=260
##nstructures=16400     # 2^14
##nstructures=32780    # 2^15

# --------------------------------------------
Mfull=`expr $Ncell / 4 - $Next`
# --- make a SWNT
n=10; m=10
LFILE=log_tmp.txt
mkcnt $n $m > $LFILE
mkc60 > $LFILE
mv ${n}${m}CNT.xyz cnt.xyz

# --- get length of CNT and decide the unit size
line=`sed -n 2p cnt.xyz`; data=(`echo $line`)
Lunit=${data[5]}
Llong=`echo "$Lunit * $Ncell" | bc -l`
Lxy=`echo "$Llong * 0.5" | bc -l`

# --- output
echo ""
echo " # of unit cells : $Ncell"
echo " # of C60s       : $Mfull"
echo " box size : $Llong x $Lxy x $Lxy"

# --- make a peapod
python $UDIR/mkpeapod.py  \
    --Fswnt cnt.xyz --Mswnt $Ncell \
    --Ffull C60.xyz --Mfull $Mfull \
    --Lxy $Lxy > $LFILE

# --- make lammps input
sh $UDIR/mk_lmpin4sed.sh $tnpt $tpre $nstructures

# --- set directory for a SWNT and peapod
for m in swnt pea; do
    XFILE=${m}.xyz
    DIR1=$WORKDIR/CAL${ical}/${m}
    xyztolmposi $XFILE > $LFILE
    if [ ! -e $DIR1 ]; then
        mkdir -p $DIR1
    fi
    mv data.lammps $XFILE band_${m}.in $DIR1
    cp $UDIR/opt.tersoff $DIR1
    cp $UDIR/sed_all.in $DIR1
    
done
rm temp.xyz cnt.xyz C60.xyz

# ---- make input scripts
cd $WORKDIR/CAL${ical}
sh ../$UDIR/mkscripts.sh $ical


