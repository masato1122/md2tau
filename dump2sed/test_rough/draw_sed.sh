wmax=54
lw1=4.0; lw2=3.0; ps=0.8;
f1=18; f2=18; f3=16; f4=12

#--- directory
FDIR=./FIGSED
if [ ! -e $FDIR ]; then
    mkdir $FDIR
fi

#--- variable
if [ $# -ge 1 ]; then
    FILE=$1
else
    echo -n "Input file name: "
    read FILE
fi

#- TYPE: 0.{1,2,3} 1.{1,3,4}
if [ $# -ge 2 ]; then
    TYPE=$2
else
    TYPE=1
fi

#--- 2.band
cb0=*; cb1=*

if [ $TYPE == "1" ]; then
    l1=1; l2=2; l3=3; cb1=1e5
else
    l1=1; l2=3; l3=4; cb1=1e6
fi

OUT=out_${FILE}.eps
gnuplot<<EOF
set terminal postscript eps color enhanced "Helvetica" $f1
set key font "Helvetica,$f2"
set tics font "Helvetica,$f3"
set tics scale 0.5
set size ratio 2.5 0.8

set label 1 at graph -0.1,1.04 right "$LAB"

#--- for pm3d
#set ylabel "{/Symbol w} (THz)" offset 1.0,0.0
#set ytics 10 offset 0.8,0.0
#set xtics 0.5 offset 0.0,0.8
#--- for plot
set ylabel "{/Symbol w} (THz)" offset 2.2,0.0
set ytics 10 offset 0.6,0.0
set xtics 0.5 offset 0.0,0.2

set xrange [0.0:1]
set mxtics 5
set xtics ('0' 0, '1' 1)

#set yrange [0:${wmax}]
set mytics 5
set border lw 0.5

#set pm3d map
#set palette defined(1 "#000066", 2 "#1e90ff", 10 "#FFCC00")
set palette rgbformula 22,13,-31
set logscale cb

#set cbrange [$cb0:$cb1]
#unset colorbox

#-- for pm3d
#set colorbox user horizontal size 0.1405,0.03 origin 0.4055,0.785
#set colorbox user horizontal size 0.1310,0.03 origin 0.3842,0.755
#-- for plot
set colorbox user horizontal size 0.174,0.03 origin 0.3505,0.735

#set format cb "10^{%L}"
set format cb "%L"
set cbtics offset 0.0,2.333 font "Helvetica,$f4"
set mcbtics 10

set out "${OUT}"
plot "$FILE" every 1 u ${l1}:${l2}:${l3} w p pt 5 ps 0.5 lc pal notitle
EOF
#splot "$FILE" every 100 u ${l1}:${l2}:${l3} w pm3d lc pal notitle
#EOF

ps2eps -l -f -B ${OUT}; mv ${OUT}.eps ${OUT}
PFIG=$FDIR/${FILE%.txt}.png
convert -transparent white -density 300 ${OUT} $PFIG
echo $PFIG
rm ${OUT}

