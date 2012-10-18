#!/bin/bash

[ $# -eq 0 ] && { echo "Usage: $0 coeff_vector pause_time "; exit 1; }

coe=$1
pause=$2

echo "#!/usr/bin/gnuplot -persist" > plot.vectff.gpi
echo "reset" >> plot.vectff.gpi
echo "a=a+1" >> plot.vectff.gpi
echo "set xr [-5:5]" >> plot.vectff.gpi
echo "set yr [-5:5]" >> plot.vectff.gpi
for ((i=1;i<=1500;i++))
do
echo "p 'VECTFF' index $i u 1:2:(\$3*$coe):(\$4*$coe) w vectors" >> plot.vectff.gpi
echo "pause $pause" >> plot.vectff.gpi
done
echo "if(a<2) reread" >> plot.vectff.gpi

echo "#!/usr/bin/gnuplot -persist" > plot
echo "reset" >> plot
echo "set term x11" >> plot
echo "unset key" >> plot
echo "a=0;" >> plot
echo "load 'plot.vectff.gpi'" >> plot

chmod u+x plot
./plot


