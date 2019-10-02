#!/bin/bash
#=================================================================

#setenv WKDIR `pwd`

f95 -o parameter writeParameter.f90

j=1
while [ $j -le 1000 ]
do
   echo "$j" > row.txt
   ./parameter.exe
   cat lake.inc.save fort.23 > lake.inc
   f95 -o lake env_sub.f90
   ./lake.exe
   mv surface.dat OUTPUT/surface-"$j".txt
   ((j++)) 
done
