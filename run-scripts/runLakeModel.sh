#!/bin/bash
#=================================================================

#setenv WKDIR `pwd`

cd Parameter-files
gfortran -o parameter writeParameter.f90
cd ..

j=1
while [ $j -le 1000 ]
do
   cd Parameter-files
   echo "$j" > row.txt
   ./parameter.exe
   cat lake.inc.save fort.23 > lake.inc
   mv lake.inc ..
   cd ..
   gfortran -c *.f90
   gfortran -o lake *.o
   ./lake.exe
   mv surface.dat surface-"$j".txt
   rm lake.inc
   ((j++)) 
done