#!/bin/bash

EXE=../../bin/./mdff.x


echo "================================================================================================================================="
echo "EXAMPLE 5 : Density of states of FCC Lennard-Jones"
echo "================================================================================================================================="
echo ""
echo "generate the fcc structure"
$EXE control_md.F > stdout.md
mv OUTFF OUTFF.md
mv TRAJFF TRAJFF.fcc
cp TRAJFF.fcc TRAJFF

echo "optimisation of the structure"
echo "optalgo=sastry"
$EXE control_sastry.F > stdout.sastry
mv ISTHFF ISTHFF.sastry
mv ISCFF ISCFF.sastry
tail -n2 ISTHFF.sastry
echo "optalgo=lbfgs"
$EXE control_lbfgs.F > stdout.lbfgs
mv ISTHFF ISTHFF.lbfgs
mv ISCFF ISCFF.lbfgs
tail -n2 ISTHFF.lbfgs
echo "optalgo=m1qn3"
$EXE control_m1qn3.F > stdout.m1qn3
mv ISTHFF ISTHFF.m1qn3
mv ISCFF ISCFF.m1qn3
tail -n2 ISTHFF.m1qn3
echo ""
echo "the structural config ISCFF.* should be the same"
cp ISCFF.sastry ISCFF
echo ""
echo "hessian matrix and dos calculation"
$EXE control_vib.F > stdout.vib
mv DOSFF DOSFF.vib
./plot_vectff.sh 15.0 0.001
echo ""
echo "complete DOS from kpoints in IBZKPTFF"
echo "(nkphon+1) x (nkphon+1) x (nkphon+1) k-points"
$EXE control_vib+dos.F > stdout.vib+dos
mv DKFF DKFF.vib+dos
mv DOSKFF DOSKFF.vib+dos
echo ""
echo "plot dos"
./plot_dos_ex5.sh 
echo ""
echo "Band calculation"
echo "Gamma to X"
$EXE control_vib+band_GtoX.F > stdout.vib+band_GtoX
mv DOSKFF DOSKFF.vib+band_GtoX
echo "Gamma to K"
$EXE control_vib+band_GtoK.F > stdout.vib+band_GtoK
mv DOSKFF DOSKFF.vib+band_GtoK
echo "Gamma to L"
$EXE control_vib+band_GtoL.F > stdout.vib+band_GtoL
mv DOSKFF DOSKFF.vib+band_GtoL
echo ""
echo "plot band stucture"
./plot_band_ex5.sh 
box=`grep "cell parameter" OUTFF.md | awk '{print $NF}'`
ncell=`grep Face-centered OUTFF.md | awk '{print $NF}'`
./plot_band_ex5.sh $box $ncell

echo "end of example 5"
