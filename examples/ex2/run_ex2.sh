#!/bin/bash

EXE=../../bin/./mdff.x

echo "EXAMPLE 2 : calculation of total energy of LJ clusters"
echo "reference structures and energies from :"
echo "The Cambridge Cluster Database"
echo "http://www-wales.ch.cam.ac.uk/CCD.html" 
echo ""

echo "N (atoms)      mdff.x       CCD          diff"
for (( cluster=3;cluster<150;cluster++)) 
do

	echo "CLUSTER_LJ" > tmp.file
	echo "10000.000   1" >> tmp.file
	echo "A" >> tmp.file
	echo $cluster >> tmp.file
	awk '{print $1,$2,$3,"0.0000 0.0000 0.0000"}' clusters/$cluster >> tmp.file 
	mv tmp.file POSFF


	$EXE control.F > stdout
	echo "$cluster `grep "Etot" OUTFF | awk '{print $14}'` `grep " $cluster " REFERENCE| awk '{print $2}'` " | awk '{printf(" %4i %16.5f %12.5f %12.5f\n",$1,$2,$3,$2-$3)}'

done

echo " "
echo "Lowest energy icosahedral minima at sizes with non-icosahedral global minima. " 
echo " " 
echo "N (atoms)      mdff.x       CCD          diff"
for cluster in 38 75 76 77 98 102  103  104
do
        echo "CLUSTER_LJ" > tmp.file
        echo "10000.000   1" >> tmp.file
        echo "A" >> tmp.file
        echo $cluster >> tmp.file
        awk '{print $1,$2,$3,"0.0000 0.0000 0.0000"}' clusters/$((cluster))i >> tmp.file
        mv tmp.file POSFF

        $EXE control.F > stdout
	echo " "$((cluster))i" `grep "Etot" OUTFF | awk '{print $14}'` `grep " $((cluster))i " REFERENCE| awk '{print $2}'` " | awk '{printf(" %4s %16.5f %12.5f %12.5f\n",$1,$2,$3,$2-$3)}'
done

