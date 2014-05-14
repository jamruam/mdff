#!/bin/bash

for ((tgt=400;tgt<=450;tgt++)) 
do

cat > control.F << EOF
&controltag
        calc='stochio'
&end

&stochiotag
        geo2=100.0
        target_nions = $tgt
        density = 3.66d0
&end
EOF

~/mdff/src/mdff.x control.F > log
grep EXACT log > exact
grep APPROXIMATED log > appro
if [ -s exact ]; then
	echo "found exact stochiometry for $tgt ions"
	mv log log.$tgt	
	ntype=`grep ntype log.$tgt | awk '{print $NF}'`
	total=`grep total log.$tgt | awk '{print $NF}'`
	box=`grep "cell param. a" log.$tgt | awk '{print $(NF-1)}'`
	echo "cell parameter of cubic box : " $box
	
	awk -v kel=$ntype '{   { if ($1=="#elem") {lr=NR} }
                               { for (k=1;k<=kel+1;k++) { 
                                                           if ((NR==lr+k)&&(lr!=0)) { printf("%4s",$1) } 
                                                        }
                               }
                           }' log.$tgt > atype 
	awk -v kel=$ntype '{   { if ($1=="#elem") {lr=NR} }
                               { for (k=1;k<=kel+1;k++) { 
                                                           if ((NR==lr+k)&&(lr!=0)) { printf("%4s",$2) } 
                                                        }
                               }
                           }' log.$tgt > itype 
	t=`cat atype`
	i=`cat itype`
        echo "type info : "
	echo "$t"
	echo "$i"
	echo " generate random structure for $tgt ions"
	random_struct.py -n $total -i "$i" -t "$t" -a "$box 0.0 0.0" -b "0.0 $box 0.0" -c "0.0 0.0 $box" > log.random
	mv POSFF.randomPY POSFF.randomPY.$tgt
	fftocar -i POSFF.randomPY.$tgt -o POSFF.randomPY.$tgt.vasp 
fi
if [ -s appro ]; then
	echo "found approximated stochiometry for $tgt ions"
fi

done
rm control.F


