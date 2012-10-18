#!/bin/bash

[ $# -eq 0 ] && { echo "Usage: $0 box ncell "; exit 1; }

pi=3.1415926535897932384626
box=$1
ncell=$2
kcell=`echo "scale=10;2*$pi*$ncell/$box" | bc`

echo "#!/usr/bin/gnuplot -persist" > plot.band
echo "set term postscript enhanced eps color 14" >> plot.band
echo 'set output "band.eps"' >> plot.band

#echo "set style line 1 lw 2 lt 1 lc 1" >> plot.band
#echo "set style line 2 lw 2 lt 1 lc 2" >> plot.band
#echo "set style line 3 lw 2 lt 1 lc 3" >> plot.band
echo "set style line 1 lw 2 lt 1 " >> plot.band
echo "set style line 2 lw 2 lt 1 " >> plot.band
echo "set style line 3 lw 2 lt 1 " >> plot.band

echo "set xr[0:2*$kcell+0.5*$kcell]" >> plot.band

echo "set xlabel 'energy {/Symbol w}'" >> plot.band
echo 'set format x ""' >> plot.band
echo "set xtics ('{/Symbol G}' 0,'{/Symbol G}' 2*$kcell,'X' $kcell,'L' 2*$kcell+0.5*$kcell)" >> plot.band

echo "p 'DOSKFF.vib+band_GtoX' u 2:5 w l title 'T1' ls 1 ,\
'DOSKFF.vib+band_GtoX' u 2:6 w l title 'T2' ls 2 ,\
'DOSKFF.vib+band_GtoX' u 2:7 w l title 'L'  ls 3 ,\
'DOSKFF.vib+band_GtoK' u (-\$2+2*$kcell):5 w l title '' ls 1 ,\
'DOSKFF.vib+band_GtoK' u (-\$2+2*$kcell):6 w l title '' ls 2 ,\
'DOSKFF.vib+band_GtoK' u (-\$2+2*$kcell):7 w l title '' ls 3 ,\
'DOSKFF.vib+band_GtoL' u ( \$2+2*$kcell):5 w l title '' ls 1 ,\
'DOSKFF.vib+band_GtoL' u ( \$2+2*$kcell):6 w l title '' ls 2 ,\
'DOSKFF.vib+band_GtoL' u ( \$2+2*$kcell):7 w l title '' ls 3" >> plot.band
echo "  #    EOF" >> plot.band
chmod u+x plot.band
./plot.band
gv band.eps

