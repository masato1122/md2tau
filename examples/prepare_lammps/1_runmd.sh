#!/bin/sh
WORKDIR=$(cd $(dirname $0); pwd)

for ical in `seq 1 5`; do
    dir1=$WORKDIR/CAL${ical}
    cd $dir1
    sh ../mkseki1.sh ${ical}
    #qsub 1_seki.sh
done

