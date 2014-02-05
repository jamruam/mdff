#!/bin/bash

EXEDLPOLY=/home/filipe/recherche/dl_class_1.9/execute/DLPOLY.Y
EXEMDFF=/home/filipe/recherche/mdff/bin/mdff.x
EXEMOLDY=/home/filipe/recherche/moldy/moldy

dlpoly=false
moldy=false
mdff=true

if $dlpoly; then
	echo "running DL_POLY test (ex7)"
	cd dl_poly
	rm STATIS OUTPUT RDFDAT REVCON  REVIVE
	$EXEDLPOLY
	tail -n +4 RDFDAT > tmp 
	mv tmp RDFDAT 
	tail -n +3 STATIS > tmp
	mv tmp STATIS
	./read_statis.sh > ene
	cd ..
fi

if $moldy; then
	echo "running moldy test (ex7)"
	cd moldy
	$EXEMOLDY control.argon > output.argon
	./plotrdf output.argon > gr.moldy
	awk 'BEGIN {lr=3} {if ((NF==13)&&($1!="Dip")&&(lr%3==0)) {print $0 } if ((NF==13)&&($1!="Dip")) {lr=lr+1} }' output.argon > ene
	cd ..
fi

if $mdff; then
	echo "running mdff test (ex7)"
	cd mdff
	$EXEMDFF control_md.F > stdout_md
	grep tot OSZIFF | awk '{print $3,$15}' > ene
	$EXEMDFF control_gr.F > stdout_gr
	cd ..
fi


cat > plot.ene << eof
#!/usr/bin/gnuplot -persist
reset
#-13921.0267036681
p 'dl_poly/ene' w l,'mdff/ene' u 1:(\$2-13921.0267036681) w l,'moldy/ene' u 0:(\$3*120) w l
eof
chmod u+x plot.ene
./plot.ene

cat > plot.gr << eof
#!/usr/bin/gnuplot -persist
reset
p 'dl_poly/RDFDAT' w l ,'mdff/GRTFF' w l,'moldy/gr.moldy' w l
eof
chmod u+x plot.gr
./plot.gr

