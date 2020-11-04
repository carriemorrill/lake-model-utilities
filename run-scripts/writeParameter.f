      program writeParameter

      open(unit=20,file="lake-params.txt",status="old")  ! latin hypercube scalings
      open(unit=22,file="row.txt",status="old")     ! row of LHS file to use, from shell script

      read(22,*) inum

      do j=1,inum
        read(20,*) v1, v2, v3, v4, v5, v6, v7  ! CHANGE this line to match number of parameters in lake-params.txt
!  80    format(a7,1x,f12.10,1x,f12.10,1x,f12.10)
      end do

! CHANGE following lines to specify ranges for parameters you wish to vary
      v1 = 2.e-3*v1 + 1.e-3   ! CDRN from 1 to 3 e-3
      v2 = 0.3*v2 + 0.2       ! eta from 0.2 to 0.5
      v3 = 0.2*v3 + 0.7       ! albsnow from 0.7 to 0.9
      v4 = 0.3*v4 + 0.4       ! albslush from 0.4 to 0.7
      v5 = 2.e6*v5 + 2.e6     ! csed from 2e6 to 4e6
      v6 = 2.0*v6 + 0.5       ! condsed from 0.5 to 2.5
      v7 = 0.15*v7 + 0.05     ! albsed from 0.05 to 0.2

! CHANGE following lines to write lake.inc parameter statements for parameters you wish to vary
      write(23,*) "      parameter (cdrn = ", v1,")"
      write(23,*) "      parameter (eta = ", v2,")"
      write(23,*) "      parameter (alb_snow = ", v3,")"
      write(23,*) "      parameter (alb_slush = ", v4,")"
      write(23,*) "      parameter (csed = ", v5,")"
      write(23,*) "      parameter (condsed = ", v6,")"
      write(23,*) "      parameter (albsd = ", v7,")"

      end
