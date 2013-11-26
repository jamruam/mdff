#!/bin/bash

EXEDLPOLY=/home/filipe/dl_poly2/src/dl_class_1.9/execute/DLPOLY_PARA.X
EXEMDFF=/home/filipe/mdff/bin/mdff.x

cd dl_poly
rm STATIS HISTORY OUTPUT RDFDAT REVCON  REVIVE
mpirun -np 4 $EXEDLPOLY
tail -n +4 RDFDAT > tmp 
mv tmp RDFDAT 
tail -n +3 STATIS > tmp
mv tmp STATIS
../read_statis.sh > ene
cd ..

cd mdff
mpirun -np 4 $EXEMDFF control_md.F > stdout_md
grep tot OSZIFF | awk '{print $3,$15}' > ene
mpirun -np 4 $EXEMDFF control_gr.F
cd ..

cat > plot.ene << eof
#!/usr/bin/gnuplot -persist
reset
#-13921.0267036681
p 'dl_poly/ene','mdff/ene' w l
eof
chmod u+x plot.ene
./plot.ene

cat > plot.gr << eof
#!/usr/bin/gnuplot -persist
reset
p 'dl_poly/RDFDAT' w l ,'mdff/GRTFF' w l
eof
chmod u+x plot.gr
./plot.gr

