#Powershell script
#=================================================================

cd Parameter-files
gfortran -o parameter writeParameter.f90
cd ..

$j=1
while($j -le 2){
   cd Parameter-files
   Set-Content .\row.txt "$j"
   ./parameter.exe
   Get-Content lake.inc.save,fort.23 | Set-Content lake.inc
   mv lake.inc ..
   cd ..
   gfortran -c *.f90
   gfortran -o lake *.o
   ./lake.exe
   mv surface.dat surface-"$j".txt
   rm lake.inc
   $j++
}   