wmax=55
lw1=4.0; lw2=3.0; ps=0.8;
f1=18; f2=18; f3=16; f4=12

FILE=SED_rtz1_PEAKS.txt

cb0=0.1
cb1=10

gnuplot<<EOF
set size ratio 3.0 0.9

set label 1 at graph -0.1,1.04 right "$LAB"

set xrange [0.0:1]
set xtics 0.5 offset 0.0,0.8
set mxtics 5
set xtics ('0' 0, '1' 1)

set ylabel "{/Symbol w} (THz)" offset 1.4,0.0
set yrange [0:${wmax}]
set ytics 10 offset 0.8,0.0
set mytics 5
set border lw 0.5

set palette model HSV rgbformulae 3,2,2
set palette maxcolors 10
set cbrange [$cb0:$cb1]
set logscale cb

set out "out.eps"
plot for [i=0:640] "$FILE" every :::i::i u 3:5:6 w p ls 1 lc pal pt 6 ps 0.2 notitle
EOF

#splot "$FILE" u 1:(\$3*1e-12):4 w pm3d lc pal notitle
#EOF

ps2eps -l -f -B out.eps; mv out.eps.eps out.eps
PFIG=peak.png
convert -density 400 out.eps $PFIG
echo $PFIG
rm out.eps


