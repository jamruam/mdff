#!/bin/bash

EXE=../../bin/./mdff.x


echo "================================================================================================================================="
echo "================================================================================================================================="
echo ""
echo "EXAMPLE 4 : testing electric field gradient"
echo ""
echo "================================================================================================================================="
echo "================================================================================================================================="
echo ""
echo "First example: Two charges in a box=10"
echo " " 
cp TRAJFF.twocharges TRAJFF
echo "&fieldtag"
echo "      qch(1) =  1.0"
echo "      qch(2) = -1.0"
echo "&end"
echo " "
echo "================================================================================================================================="
echo ""
echo "Direct summation"
echo "&eftag (main parameters): "
echo "ncelldirect = 40"
echo "cutefg = 1000.0d0"
echo "&end"
echo ""
cat control_direct.F > control.F
echo "&fieldtag" >> control.F
echo "      qch(1) =  1.0" >> control.F
echo "      qch(2) = -1.0" >> control.F
echo "&end" >> control.F
$EXE control.F > stdout_direct_twocharges
cat EFGALL > EFGALL.direct_twocharges
tail -n3 EFGALL 
echo ""


echo "================================================================================================================================="
echo ""
echo "Ewald summation"
echo "&eftag (main parameters): "
echo "ncellewald = 11"
echo "alphaES = 0.8"
echo "&end"
echo ""

cat control_ewald.F > control.F
echo "&fieldtag" >> control.F
echo "      qch(1) =  1.0" >> control.F
echo "      qch(2) = -1.0" >> control.F
echo "&end" >> control.F
$EXE control.F > stdout_ewald_twocharges
cat EFGALL > EFGALL.ewald_twocharges
tail -n3 EFGALL 
echo ""
echo "================================================================================================================================="
echo ""
echo "From Nymand and Linse, J. Chem. Phys. 112, 6152 (2000)"
echo ""
echo "Method           V1,xx                 V1,yy"
echo "ES             2.0003749            -1.0001874" 
echo "DS             2.0003749            -1.0001874"
echo ""
echo ""
echo "================================================================================================================================="
echo "================================================================================================================================="
echo ""
echo "Second example: cluster (9 atoms)"
echo ""
cp TRAJFF.cluster TRAJFF
echo ""
echo "&fieldtag"
echo "      qch(1) =  0.1  ! A particules"
echo "      qch(2) = -0.8  ! B particules"
echo "&end"
echo " "
echo "================================================================================================================================="
echo ""
echo "Direct summation"
echo "&eftag (main parameters): "
echo "ncelldirect = 40"
echo "cutefg = 1000.0d0"
echo "&end"
echo ""
cat control_direct.F > control.F
echo "&fieldtag" >> control.F
echo "      qch(1) =  0.1" >> control.F
echo "      qch(2) = -0.8" >> control.F
echo "&end" >> control.F
$EXE control.F > stdout_direct_cluster
cat NMRFF  > NMRFF.direct_cluster
tail -n10 NMRFF 
echo ""


echo "================================================================================================================================="
echo ""
echo "Ewald summation"
echo "&eftag (main parameters): "
echo "ncellewald = 11"
echo "alphaES = 0.8"
echo "&end"
echo ""

cat control_ewald.F > control.F
echo "&fieldtag" >> control.F
echo "      qch(1) =  0.1" >> control.F
echo "      qch(2) = -0.8" >> control.F
echo "&end" >> control.F
$EXE control.F > stdout_ewald_cluster
cat NMRFF  > NMRFF.ewald_cluster
tail -n10 NMRFF 
echo ""
echo "================================================================================================================================="
echo ""
echo "GULP output (note : e²=14.3998 eV/A)"
echo "EFG Tensor properties :"
echo ""
echo "  -------------------------------------------------------------------------------"
echo "     Site no.    Diagonalised EFGs (V/Angs**2)          eVzz/h        Asymmetry"
echo "                   xx        yy        zz             (MHz/barn)      Parameter"
echo "  -------------------------------------------------------------------------------"
echo "         1       15.2454   16.4848  -31.7301            76.7225        0.0391"
echo "         2       13.1974   15.5994  -28.7968            69.6297        0.0834"
echo "         3       20.5185   23.5380  -44.0565           106.5272        0.0685"
echo "         4       12.0418   14.6251  -26.6669            64.4796        0.0969"
echo "         5        8.6231   10.8748  -19.4979            47.1452        0.1155"
echo "         6       22.1626   24.1320  -46.2946           111.9390        0.0425"
echo "         7        7.7041    9.9389  -17.6430            42.6603        0.1267"
echo "         8       16.7300   17.0914  -33.8214            81.7789        0.0107"
echo "         9       -0.8228   -3.4544    4.2772            10.3421        0.6153"
echo "  -------------------------------------------------------------------------------"

echo "                   xx        yy        zz             (MHz/barn)      Parameter "  >  gulp.out
echo "         1       15.2454   16.4848  -31.7301            76.7225        0.0391"     >> gulp.out
echo "         2       13.1974   15.5994  -28.7968            69.6297        0.0834"     >> gulp.out
echo "         3       20.5185   23.5380  -44.0565           106.5272        0.0685"     >> gulp.out
echo "         4       12.0418   14.6251  -26.6669            64.4796        0.0969"     >> gulp.out
echo "         5        8.6231   10.8748  -19.4979            47.1452        0.1155"     >> gulp.out
echo "         6       22.1626   24.1320  -46.2946           111.9390        0.0425"     >> gulp.out
echo "         7        7.7041    9.9389  -17.6430            42.6603        0.1267"     >> gulp.out
echo "         8       16.7300   17.0914  -33.8214            81.7789        0.0107"     >> gulp.out
echo "         9       -0.8228   -3.4544    4.2772            10.3421        0.6153"     >> gulp.out

echo "                   EWALD(mdff)                      |                    DIRECT (mdff)                  |                       GULP" > COMP
echo "       vxx          vyy          vzz          eta   |      vxx          vyy          vzz          eta   |      vxx          vyy          vzz          eta" >> COMP
echo "---------------------------------------------------------------------------------------------------------------------------------------------------------------">> COMP
paste NMRFF.ewald_cluster NMRFF.direct_cluster gulp.out | tail -n9 | awk '{printf("%12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s\n",$3,$4,$5,$6,$9,$10,$11,$12,$14,$15,$16,$18)}' >> COMP

echo ""
echo "look at COMP file for direct comparison"
