wmax=55
lw1=4.0; lw2=3.0; ps=0.8;
f1=18; f2=18; f3=16; f4=12

FILE=SED_rtz1.txt

cb0=1e-1
cb1=100

gnuplot<<EOF
#set terminal postscript eps color enhanced "Helvetica" $f1
#set key font "Helvetica,$f2"
#set tics font "Helvetica,$f3"
#set ticscale 0.5

set size ratio 3.0 0.9

set xrange [0.0:1]
set xtics 0.5 offset 0.0,0.8
set mxtics 5
set xtics ('0' 0, '1' 1)

set ylabel "{/Symbol w} (THz)" offset 1.4,0.0
set yrange [0:${wmax}]
set ytics 10 offset 0.8,0.0
set mytics 5
set border lw 0.5

set pm3d map
set logscale cb
set palette defined(1 "#000066", 2 "#1e90ff", 10 "#FFCC00")
set palette rgbformula 22,13,-31
set palette maxcolors 10

set cbrange [$cb0:$cb1]
#unset colorbox

#set colorbox user horizontal size 0.1405,0.03 origin 0.4055,0.785
set colorbox user horizontal size 0.1310,0.03 origin 0.3842,0.755
set format cb "10^{%L}"
set cbtics offset 0.0,2.4 font "Helvetica,$f4"
set mcbtics 10

set out "out.eps"
splot "$FILE" u 1:(\$3*1e-12):4 w pm3d lc pal notitle
EOF

#splot "$FILE" u 1:(\$3*1e-12):4 w pm3d lc pal notitle
#EOF

ps2eps -l -f -B out.eps; mv out.eps.eps out.eps
PFIG=sed.png
convert -density 400 out.eps $PFIG
echo $PFIG
rm out.eps


