#!/bin/bash
EXE=../../bin/mdff.x
sep="=================================================="
echo $sep
echo "#Example 0: Minimum settings with FCC structure"
echo "#just look at control.F"
echo $sep
echo "FCC LJ a=1.5414"
echo $sep
$EXE control.F > stdout 
grep natm stdout
etot=`grep Etot stdout | awk '{print $NF}'` 
etot2=`echo ${etot} | sed -e 's/[eE]+*/\\*10\\^/'`
natm=`grep natm stdout | awk '{print $NF}'`
lr=`grep "long range correction :" stdout |  awk '{print $NF}'`
etotr=`echo "scale=8;$etot2/$natm" | bc -l`
etotlr=`echo "scale=8;($etot2+($lr))/$natm " | bc -l`
echo "Etot                  = " $etot
echo "Etot/atom             = " $etotr
echo "Etot/atom + lr corre. = " $etotlr
exit 0
