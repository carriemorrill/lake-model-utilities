#Powershell script
#=================================================================

gfortran -o parameter writeParameter.f90

$j=1
while($j -le 1000){
   Set-Content .\row.txt "$j"
   ./parameter.exe
   Get-Content lake.inc.save,fort.23 | Set-Content lake.inc
   gfortran -o lake HB-windy-sed-varD2.f90
   ./lake.exe
   mv surface.dat surface-"$j".txt
   $j++
}   
