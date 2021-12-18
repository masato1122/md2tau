#!/bin/bash
if [ $# -ge 1 ]; then
    FILE=$1
else
    echo -n "Input file name: "
    read FILE
fi

if [ ! -e $FILE ]; then
    echo "cannot find $FILE"
    exit
fi

y0=1e-3
y1=10


lw=1.0; ps=0.8;
gnuplot<<EOF
set size ratio 0.3 0.7

set xlabel ""
set xrange [*:*]
#set xtics $xtics; set format x ""

set ylabel ""
set yrange [$y0:$y1]
set logscale y; set format y "10^{%L}"

set out "out.eps"
plot \
for [i=2:11] "$FILE" u 1:i w l ls i lw 0.7 dt (0,0) title "",\
"$FILE" u 1:2 w lp ls 2 lw 1.0 ps 0.3 dt (0,0) title ""
EOF

EFIG=data.eps
ps2eps -l -f -B out.eps; mv out.eps.eps $EFIG; epstopdf $EFIG
echo $EFIG

#PFIG=fig.png
#convert -trim -transparent white -density 400 out.eps $PFIG
#echo $PFIG

