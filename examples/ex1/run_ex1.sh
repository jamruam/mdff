#!/bin/bash

echo "#Example 1: LJ fcc structure at low temperature"
echo "The configuration is readed in POSFF"

EXE=../../src/mdff.x 
mpirun -np 2 $EXE control.F
