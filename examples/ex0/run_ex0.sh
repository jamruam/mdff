#!/bin/bash
EXE=../../bin/mdff.x
sep="=================================================="
echo $sep
echo "#Example 0: Minimum settings with FCC structure"
echo "#just look at control.F"
echo $sep
echo "4x4x4 FCC LJ a=6.0000000000"
echo $sep
$EXE control.F > stdout.4x4x4 
grep Etot stdout.4x4x4 
exit 0
