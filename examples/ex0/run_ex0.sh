#!/bin/bash
EXE=../../src/mdff.x
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
etot3=`echo "scale=8;$etot2*1.0" | bc -l`
etotlr=`echo "scale=8;($etot2)+($lr)" | bc -l`
etotr=`echo "scale=8;$etot2/$natm" | bc -l`
etotrlr=`echo "scale=8;($etot2+($lr))/$natm " | bc -l`
echo "Etot                  = " $etot3
echo "Etot      + lr corre. = " $etotlr
echo ""
etotdl=`head -n 4 dl_poly/STATIS | tail -n 1 | awk '{print $1}'`
echo "dl_poly               = " $etotdl 
exit 0
