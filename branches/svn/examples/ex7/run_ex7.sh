#!/bin/bash

EXEDLPOLY=
EXEMDFF=mdff.x

dlpoly=false
mdff=true
ensemble=nve

if $dlpoly; then
	echo "running DL_POLY test (ex7)"
	cd dl_poly
		rm STATIS OUTPUT RDFDAT REVCON  REVIVE
		$EXEDLPOLY
		tail -n +4 RDFDAT > RDFDAT.$ensemble
	cd ..
fi

if $mdff; then
	echo "running mdff test (ex7)"
	cd mdff
	$EXEMDFF control_md.F > stdout_md.$ensemble
	$EXEMDFF control_gr.F > stdout_gr.$ensemble
	mv GRTFF GRTFF.$ensemble
	cd ..
fi

cat > plot.gr << eof
#!/usr/bin/gnuplot -persist
reset
p 'dl_poly/RDFDAT.$ensemble' w l ,'mdff/GRTFF.$ensemble' w l
eof
chmod u+x plot.gr
./plot.gr

