#!/bin/bash

[ $# -eq 0 ] && { echo "Usage: $0 coeff_vector pause_time "; exit 1; }

coe=$1
pause=$2
i=767
#first mode is 0

echo "#!/usr/bin/gnuplot -persist"                                > plot.vectff
echo "reset"                                                     >> plot.vectff
echo "set term x11"                                              >> plot.vectff
echo "set xr [-4:4]"                                             >> plot.vectff
echo "set yr [-4:4]"                                             >> plot.vectff
echo "unset key"                                                 >> plot.vectff
echo "set title 'highest energy mode = $((i+1))' "               >> plot.vectff
echo "p 'VECTFF' index $i u 1:2:(\$3*$coe):(\$4*$coe) w vectors" >> plot.vectff
chmod u+x plot.vectff
./plot.vectff


