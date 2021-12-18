#!/usr/bin/env bash
WORKDIR=$(cd $(dirname $0); pwd)

FILE="check.txt"

lw=1.0; ps=0.8;
gnuplot<<EOF
set xlabel ""
set xrange [0:55]
#set xtics $xtics; set format x ""

set ylabel ""
set yrange [*:1000]
#set ytics $ytics; set format y ""
set logscale y

set out "out.eps"
plot \
"$FILE" u 1:2 w lp ls 2 lw $lw ps $ps dt (0,0) title "",\
"$FILE" u 1:3 w l  ls 3 lw $lw dt (0,0) title ""
EOF

EFIG=fig.eps
ps2eps -l -f -B out.eps; mv out.eps.eps $EFIG; epstopdf $EFIG
echo $EFIG
rm $EFIG

