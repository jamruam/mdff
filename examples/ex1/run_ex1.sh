#!/bin/bash

echo "#Example 1: LJ fcc structure at low temperature"
echo "The configuration is readed in POSFF"
EXE=../../bin/mdff.x 
/usr/local/bin/mpirun -n 2 $EXE control.F
