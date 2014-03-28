#!/bin/bash

echo "#Example 1: LJ (argon) fcc structure at low temperature"
echo "The configuration is readed in POSFF"
echo "more info in control.F"
EXE=../../bin/mdff.x 
mpirun -n 2 $EXE control.F | tee stdout
