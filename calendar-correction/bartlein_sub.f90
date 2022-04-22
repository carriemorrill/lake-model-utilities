subroutine cal_adjust_PMIP (variable, vfill, endageBP, begageBP, begyrCE, agestep, nsimyrs, nt, calendar_type, var3d_in, var3d_out)

!! subroutines from Bartlein and Shafer 2019 for paleo calendar corrections.
!! Options for daily data and 4D data are currently commented out, such that
!! input values must be a monthly timeseries from one location/grid cell.
!! For transient time series (as opposed to equilibrium, where year is repeated)
!! make sure that orb params and month lengths are being updated each year.

implicit none

character(64), intent(in)    :: variable             ! CMIP/PMIP variable name (e.g. "tas", "pr")
real(8), intent(in)          :: vfill                ! fill value
integer(4), intent(in)       :: begageBP             ! beginning year (BP) (negative, e.g. 10 ka = -10000 BP)
integer(4), intent(in)       :: endageBP             ! ending year (BP)
integer(4), intent(in)       :: begyrCE              ! beginning (pseudo-) year of individual model simulation
integer(4), intent(in)       :: agestep              ! age step size (interval between age calculations)
integer(4), intent(in)       :: nsimyrs              ! number of years of simulation for each age
integer(4), intent(in)       :: nt                   ! total nmber of months (nt = ny*nm)
character(32), intent(in)    :: calendar_type        ! calendar type
!real(4), allocatable, intent(in)    :: var3d_in(:,:,:)      ! 3-D (e.g. nlon,nlat,ndtot) input data 
!real(4), allocatable, intent(out)    :: var3d_out(:,:,:)    ! 3-D (e.g. nlon,nlat,nt) output adjusted data 
real(4), intent(in)    :: var3d_in(1,1,nt)      ! 3-D (e.g. nlon,nlat,ndtot) input data 
real(4), intent(out)   :: var3d_out(1,1,nt)    ! 3-D (e.g. nlon,nlat,nt) output adjusted data 

external get_month_lengths, mon_to_day_ts, day_to_mon_ts

! past ages are negative, e.g. 21 ka = 21,000 cal yr BP = -21000, and 1950 CE = 0 cal yr BP = 0
! simulation age-related variables (controls of orbital parameters)
integer(4)              :: nages                ! number of simulation ages
integer(4), allocatable :: iageBP(:)            ! (ny) year BP 1950 (negative, e.g. 1900 CE = -50 years BP 1950)

! simulation year-related variables (controls of vernal equinox day and leap-year status)
integer(4), allocatable :: iyearCE(:)           ! (ny) yearCE simulation year (e.g. 1850CE, 850CE, etc.)

! month-length variables
integer(4), allocatable :: imonlen_0ka(:,:)     ! (ny,nm) integer-value month lengths -- 0ka
integer(4), allocatable :: imonmid_0ka(:,:)     ! (ny,nm) integer-value mid days -- 0ka
integer(4), allocatable :: imonbeg_0ka(:,:)     ! (ny,nm) integer-value month beginning days -- 0ka
integer(4), allocatable :: imonend_0ka(:,:)     ! (ny,nm) integer-value month ending days -- 0ka
integer(4), allocatable :: imonlen(:,:)         ! (ny,nm) integer-value month lengths (paleo)
integer(4), allocatable :: imonmid(:,:)         ! (ny,nm) integer-value mid days (paleo)
integer(4), allocatable :: imonbeg(:,:)         ! (ny,nm) integer-value month beginning (paleo)
integer(4), allocatable :: imonend(:,:)         ! (ny,nm) integer-value month ending days (paleo)
real(8), allocatable    :: rmonlen(:,:)         ! (ny,nm) real-value month lengths (paleo)
real(8), allocatable    :: rmonmid(:,:)         ! (ny,nm) real-value mid days (paleo)
real(8), allocatable    :: rmonbeg(:,:)         ! (ny,nm) real-value month beginning (paleo)
real(8), allocatable    :: rmonend(:,:)         ! (ny,nm) real-value month ending days (paleo)
real(8), allocatable    :: VE_day(:)            ! (ny) vernal equinox day in simulation year
real(8), allocatable    :: SS_day(:)            ! (ny) (northern) summer solstice day in simulation year
integer(4), allocatable :: ndays(:)             ! (ny) number of days in year
integer(4), allocatable :: imonlen_0ka_ts(:)    ! (nt) integer-value month lengths at present as time series
real(8), allocatable    :: rmonmid_ts(:)        ! (nt) real-value paleo month mid days as time series
real(8), allocatable    :: rmonbeg_ts(:)        ! (nt) real-value paleo month beginning days as time series
real(8), allocatable    :: rmonend_ts(:)        ! (nt) real-value paleo month ending as time series
integer(4), allocatable :: imonmid_ts(:)        ! (nt) integer-value paleo month mid days as time series
integer(4), allocatable :: imonbeg_ts(:)        ! (nt) integer-value paleo month beginning days as time series
integer(4), allocatable :: imonend_ts(:)        ! (nt) integer-value paleo month ending as time series
!integer(4), allocatable :: ndays_ts(:)          ! (nt) integer-value times series of paleo year lengths
real(8), allocatable    :: mon_time(:)          ! (nt) new monthly time values for daily-input files
real(8), allocatable    :: mon_time_bnds(:,:)   ! (2,nt) new monthly time-bounds values for daily input files

!character(8)            :: time_freq            ! type of CMIP/PMIP time frequency (e.g. Aclim, Amon, day, etc.)

! data
integer(4)              :: ny                   ! number of years and total nmber of months (nt = ny*nm)
integer(4), parameter   :: nm=12                ! number of months per year
integer(4)              :: ndtot,ndtot_0ka      ! total number of days
real(8), allocatable    :: xdh(:,:)             ! 3-D (e.g. ivar_dimlen(1),ndtot) pseudo- or actual daily data
real(8), allocatable    :: var_adj(:,:)         ! 3-D (e.g. ivar_dimlen(1),nt) adjusted data 
!real(8), allocatable    :: var4d_in(:,:,:,:)    ! 4-D (e.g. nlon,nlat,nlev,nt) input data
!real(4), allocatable    :: var4d_out(:,:,:,:)   ! 4-D (e.g. nlon,nlat,nlev,nt) output adjusted data

integer, parameter      :: maxdims = 5
!integer(4)              :: invar_dimlen(maxdims)

! smoothing parameters for multi-year pseudo-daily interpolation
integer(4)              :: nw_tmp=21, nsw_tmp=20    ! smoothing parameters
logical                 :: smooth=.true., restore=.true.
logical                 :: no_negatives = .false.   ! restrict pseudo-daily interpolated values to positive values?

integer(4)              :: n,m,j,k,i          ! indices

!*************************************************************************
! Step 1:  Allocate month-length arrays
!*************************************************************************

!*******************************************************
! these will get read in from python
!    variable = 'tas'
!    endageBP = -6000.
!    begageBP = -6000.
!    agestep = 1.   ! use 10 for decadal averages, etc.
!    nsimyrs = 200.   ! number of sim yrs per time, 1 for transient, >1 for equilibrium
!    calendar_type = '365_day'
!************************************************************

    no_negatives = .false.
    select case (trim(variable))
    case ('pr','clt','sic')
        no_negatives = .true.
    case default
        continue
    end select

! allocate month-length arrays
    nages = (endageBP - begageBP)/agestep + 1
    ny = nages * nsimyrs
!    nt = ny * nm
    allocate (iageBP(ny), iyearCE(ny))
    allocate (imonlen_0ka(ny,nm),imonmid_0ka(ny,nm),imonbeg_0ka(ny,nm),imonend_0ka(ny,nm))
    allocate (imonlen(ny,nm), imonmid(ny,nm), imonbeg(ny,nm), imonend(ny,nm))
    allocate (rmonlen(ny,nm), rmonmid(ny,nm), rmonbeg(ny,nm), rmonend(ny,nm))
    allocate (VE_day(ny), SS_day(ny), ndays(ny), imonlen_0ka_ts(nt))
    allocate (imonmid_ts(nt), imonbeg_ts(nt),imonend_ts(nt), rmonmid_ts(nt), rmonbeg_ts(nt), rmonend_ts(nt))
    allocate (mon_time(nt), mon_time_bnds(2,nt))

!*************************************************************************
! Step 2:  get month lengths
!*************************************************************************
! 0 ka month lengths (used in pseudo-daily interpolation)
!    if (trim(time_freq) .ne. 'day') then
        call get_month_lengths(calendar_type, 0, agestep, nages, begyrCE, nsimyrs, &
            iageBP, iyearCE, imonlen_0ka, imonmid_0ka, imonbeg_0ka, imonend_0ka, rmonlen, rmonmid, rmonbeg, rmonend, &
                VE_day, SS_day, ndays)
!    end if

! paleo month lengths
    call get_month_lengths(calendar_type, begageBP, agestep, nages, begyrCE, nsimyrs, &
        iageBP, iyearCE, imonlen, imonmid, imonbeg, imonend, rmonlen, rmonmid, rmonbeg, rmonend, VE_day, SS_day, ndays)

!   write(*,*) "Mod_len  Mod_beg  Mod_end  Pal_len  Pal_beg  Pal_end"
!do j=1,12
!   write(*,*) imonlen_0ka(1,j),imonbeg_0ka(1,j),imonend_0ka(1,j),imonlen(1,j),imonbeg(1,j),imonend(1,j)
!end do

! reshape month lengths into time series, and get cumulative number of days in years
    ndtot_0ka = 0
    do n=1,ny
        ndtot_0ka=ndtot_0ka + ndays(n)
!        ndyr = ndyr + ndays(n)
        do m=1,nm
            i = (n-1)*nm + m
            imonlen_0ka_ts(i) = imonlen_0ka(n,m)
            rmonmid_ts(i) = rmonmid(n,m); rmonbeg_ts(i) = rmonbeg(n,m); rmonend_ts(i) = rmonend(n,m)
            imonmid_ts(i) = imonmid(n,m); imonbeg_ts(i) = imonbeg(n,m); imonend_ts(i) = imonend(n,m)
!            ndays_ts(i) = ndyr
        end do
    end do

! total number of days in simulation
    ndtot = 0
    do n=1,ny
        ndtot = ndtot + ndays(n)
    end do

!********************************************************************
! Step 3:  Get the input variable to be adjusted, allocate output arrays
!********************************************************************
! allocate variables
!    select case(invar_ndim)
!    case (3)
!        if (trim(time_freq) .eq. 'day') then
!            allocate(var3d_in(invar_dimlen(1),invar_dimlen(2),ndtot))
!        else
!!            allocate(var3d_in(invar_dimlen(1),invar_dimlen(2),nt))
!!            allocate(var3d_in(:,:,:))
!        end if
!!        allocate(var3d_out(invar_dimlen(1),invar_dimlen(2),nt))
!!          allocate(var3d_out(:,:,:))
!    case (4)
!        if (trim(time_freq) .eq. 'day') then
!            allocate(var4d_in(invar_dimlen(1),invar_dimlen(2),invar_dimlen(3),ndtot))
!        else
!            allocate(var4d_in(invar_dimlen(1),invar_dimlen(2),invar_dimlen(3),nt))
!        end if
!        allocate(var4d_out(invar_dimlen(1),invar_dimlen(2),invar_dimlen(3),nt))
!    case default
!        stop "allocating variables"
!    end select
!    allocate(xdh(invar_dimlen(1),ndtot), var_adj(invar_dimlen(1),nt))
    allocate(xdh(1,ndtot),var_adj(1,nt))

!************************************************************************
! Step 4:  Get calendar-adjusted values
!************************************************************************    
!    select case(invar_ndim)
!    case (3) 
!        ntotalpts = invar_dimlen(1) * invar_dimlen(2)
!    case (4)
!        ntotalpts = invar_dimlen(1) * invar_dimlen(2) * invar_dimlen(3)
!    end select
        
    ! Note: OpenMP seems to work best if only innermost loop is parallelized
!    select case (invar_ndim)
!    case (3)
!!        do k=1,invar_dimlen(2)
        do k=1,1
            ! if daily data, copy into xdh
!            if (trim(time_freq) .eq.'day') xdh(:,:) = var3d_in(:,k,:)
        
            !$omp parallel do
!!            do j=1,invar_dimlen(1)
              do j=1,1
!                if (trim(time_freq) .eq.'day') xdh(j,:) = var3d_in(j,k,:)
                ! unless the input data are daily, do pseudo-daily interpolation of the monthly input data
!                if (trim(time_freq) .ne. 'day') then
                    ! interpolate
                    call mon_to_day_ts(nt, imonlen_0ka_ts, dble(var3d_in(j,k,:)), dble(vfill), &
                        no_negatives, smooth, restore, ndtot, nw_tmp, nsw_tmp, xdh(j,:))
                    ! reaggregate daily data using correct calendar
                    call day_to_mon_ts(ny,ndays,rmonbeg,rmonend,ndtot,xdh(j,:),dble(vfill),var_adj(j,:))
!                else
                    ! input data are already daily, so just reaggregate using correct calendar
!                    call day_to_mon_ts(ny,ndays,rmonbeg,rmonend,ndtot,xdh(j,:),dble(vfill),var_adj(j,:))
!                end if
            
                var3d_out(j,k,:)=sngl(var_adj(j,:))

            end do
            !$omp end parallel do
        end do
        
!    case(4)
!        do l=1,invar_dimlen(3)
!            do k=1,invar_dimlen(2)      
!                ! if daily data, copy into xdh
!                if (trim(time_freq) .eq.'day') xdh(:,:) = var4d_in(:,k,l,:)
        
!                !$omp parallel do
!                do j=1,invar_dimlen(1)
!                    if (trim(time_freq) .eq.'day') xdh(j,:) = var4d_in(j,k,l,:)
!                    ! unless the input data are daily, do pseudo-daily interpolation of the monthly input data
!                    if (trim(time_freq) .ne. 'day') then
!                        ! interpolate
!                        call mon_to_day_ts(nt, imonlen_0ka_ts, dble(var4d_in(j,k,l,:)), dble(vfill), &
!                            no_negatives, smooth, restore, ndtot, nw_tmp, nsw_tmp, xdh(j,:))
!                        ! reaggregate daily data using correct calendar
!                        call day_to_mon_ts(ny,ndays,rmonbeg,rmonend,ndtot,xdh(j,:),dble(vfill),var_adj(j,:))
!                    else
!                        ! input data are already daily, so just reaggregate using correct calendar
!                        call day_to_mon_ts(ny,ndays,rmonbeg,rmonend,ndtot,xdh(j,:),dble(vfill),var_adj(j,:))
!                    end if
!                    var4d_out(j,k,l,:)=sngl(var_adj(j,:))
!                end do
!                !$omp end parallel do
!            end do
!        end do
!    case default
!        stop "adjusting"
!   end select

!************************************************************************
! Step 5:  Deallocate
!************************************************************************    

    deallocate (iageBP, iyearCE)
    deallocate (imonlen_0ka,imonmid_0ka,imonbeg_0ka,imonend_0ka)
    deallocate (imonlen, imonmid, imonbeg, imonend)
    deallocate (rmonlen, rmonmid, rmonbeg, rmonend)
    deallocate (VE_day, SS_day, ndays, imonlen_0ka_ts)
    deallocate (imonmid_ts, imonbeg_ts,imonend_ts, rmonmid_ts, rmonbeg_ts, rmonend_ts)
    deallocate (mon_time, mon_time_bnds)
!    select case (invar_ndim)
!    case (3)
!!        deallocate (var3d_in, var3d_out)
!    case (4)
!        deallocate (var4d_in, var4d_out)
!    case default
!        stop "deallocating"
!    end select
!!    deallocate (xdh, var_adj)
    
    end subroutine
    
!**********************************************************************************************
subroutine get_month_lengths(calendar_type, begageBP, agestep, nages, begyrCE, nsimyrs, &
    iageBP, iyearCE, imonlen, imonmid, imonbeg, imonend, rmonlen, rmonmid, rmonbeg, rmonend, VE_day, SS_day, ndays)
!***********************************************************************************************

    implicit none

    external GISS_orbpars, GISS_srevents, monlen, adjust_to_ref_length, adjust_to_ref_day, &
    adjust_to_yeartot, integer_monlen, imon_begmidend

    integer(4), parameter   :: nm = 12              ! number of months in the year
    integer(4), parameter   :: nd_360 = 360         ! number of days in a 360-day year
    integer(4), parameter   :: daysinmonth360 = 30  ! number of days in a month of a 360-day year
    integer(4), parameter   :: nd_365 = 365         ! number of days in a 365-day "noleaps" year
    integer(4), parameter   :: nd_366 = 366         ! number of days in a 366-day leap year

    ! other calendar-related variables
    real(8)     :: veqday_360 = 80.0d0              ! fixed vernal equinox day, 360-day year
    real(8)     :: veqday_365 = 80.5d0              ! fixed vernal equinox day, 365-day year
    real(8)     :: veqday_366 = 81.5d0              ! fixed vernal equinox day, 366-day year
    real(8)     :: ssday_360 = 170.5d0              ! fixed (northern) summer solstice day, 360-day year
    real(8)     :: ssday_365 = 173.0d0              ! fixed (northern) summer solstice day, 365-day year
    real(8)     :: ssday_366 = 173.5d0              ! fixed (northern) summer solstice day, 366-day year
!    real(8)     :: tropical_year = 365.24219876d0   ! length of a tropical year (days)
    real(8)     :: progreg_year = 365.2425d0        ! length of a Gregorian year (days)

    integer(4)  :: nd_progreg                       ! number of days in a 365 or 366-day year proleptic_gregorian calendar
!    real(8)     :: veqday_progreg                   ! vernal equinox day in a 365 or 366-day year proleptic_gregorian calendar
    real(8)     :: perihelion                       ! perihelion day (in a proleptic Gregorian calendar)
    real(8)     :: aphelion                         ! aphelion day (in a proleptic Gregorian calendar)

    ! month-length definitions
    real(8)     :: present_mon_360(nm) = 30.0d0     ! present-day month lengths in 360-day year
    real(8)     :: present_mon_noleap(nm) = &       ! present-day month lengths in 365-day (noleap) year
        (/ 31.0d0, 28.0d0, 31.0d0, 30.0d0, 31.0d0, 30.0d0, 31.0d0, 31.0d0, 30.0d0, 31.0d0, 30.0d0, 31.0d0 /)
    real(8)     :: present_mon_leap(nm) = &         ! present-day month lengths in 366-day (leap) year
        (/ 31.0d0, 29.0d0, 31.0d0, 30.0d0, 31.0d0, 30.0d0, 31.0d0, 31.0d0, 30.0d0, 31.0d0, 30.0d0, 31.0d0 /)
!    real(8)     :: present_mon_365_trop(nm) = &     ! present-day month lengths in a tropical year (note Feb.)
!        (/ 31.0d0, 28.24219876d0, 31.0d0, 30.0d0, 31.0d0, 30.0d0, 31.0d0, 31.0d0, 30.0d0, 31.0d0, 30.0d0, 31.0d0 /)
    real(8)     :: present_mon_365_progreg(nm) = &  ! present-day month lengths in a Gregorian year (note Feb.)
        (/ 31.0d0, 28.2425d0, 31.0d0, 30.0d0, 31.0d0, 30.0d0, 31.0d0, 31.0d0, 30.0d0, 31.0d0, 30.0d0, 31.0d0 /)
    
    real(8)     :: present_beg_360(nm) = &          ! present-day month beginning day in 360-day year
        (/  0.0d0, 30.0d0, 60.0d0, 90.0d0, 120.0d0, 150.0d0, 180.0d0, 210.0d0, 240.0d0, 270.0d0, 300.0d0, 330.0d0 /)
    real(8)     :: present_mid_360(nm) = &          ! present-day month middle day in 360-day year
        (/ 15.0d0, 45.0d0, 75.0d0, 105.0d0, 135.0d0, 165.0d0, 195.0d0, 225.0d0, 255.0d0, 285.0d0, 315.0d0, 345.0d0 /)
    real(8)     :: present_end_360(nm) = &          ! present-day month ending day in 360-day year
        (/ 30.0d0, 60.0d0, 90.0d0, 120.0d0, 150.0d0, 180.0d0, 210.0d0, 240.0d0, 270.0d0, 300.0d0, 330.0d0, 360.0d0 /)
    real(8)     :: present_beg_365(nm) = &          ! present-day month beginning day in 365-day (noleap) year
        (/  0.0d0, 31.0d0, 59.0d0, 90.0d0, 120.0d0, 151.0d0, 181.0d0, 212.0d0, 243.0d0, 273.0d0, 304.0d0, 334.0d0 /)
    real(8)     :: present_mid_365(nm) = &          ! present-day month middle day in 365-day (noleap) year
        (/ 15.5d0, 45.0d0, 74.5d0, 105.0d0, 135.5d0, 166.0d0, 196.5d0, 227.5d0, 258.0d0, 288.5d0, 319.0d0, 349.5d0 /)
    real(8)     :: present_end_365(nm) = &          ! present-day month ending day in 365-day (noleap) year
        (/ 31.0d0, 59.0d0, 90.0d0, 120.0d0, 151.0d0, 181.0d0, 212.0d0, 243.0d0, 273.0d0, 304.0d0, 334.0d0, 365.0d0 /)
    real(8)     :: present_beg_366(nm) = &          ! present-day month beginning day in 366-day (leap) year
        (/  0.0d0, 31.0d0, 60.0d0, 91.0d0, 121.0d0, 152.0d0, 182.0d0, 213.0d0, 244.0d0, 274.0d0, 305.0d0, 335.0d0 /)
    real(8)     :: present_mid_366(nm) = &          ! present-day month beginning day in 366-day (leap) year
        (/ 15.5d0, 45.5d0, 75.5d0, 106.0d0, 136.5d0, 167.0d0, 197.5d0, 228.5d0, 259.0d0, 289.5d0, 320.0d0, 350.5d0 /)
    real(8)     :: present_end_366(nm) = &          ! present-day month beginning day in 366-day (leap) year
        (/ 31.0d0, 60.0d0, 91.0d0, 121.0d0, 152.0d0, 182.0d0, 213.0d0, 244.0d0, 274.0d0, 305.0d0, 335.0d0, 366.0d0 /)
!    real(8)     :: present_beg_365_trop(nm) = &     ! present-day month beginning day in a tropical year
!        (/  0.0000d0, 31.0000d0, 59.2422d0, 90.2422d0, 120.2422d0, 151.2422d0, 181.2422d0, 212.2422d0, 243.2422d0, &
!            273.2422d0, 304.2422d0, 334.2422d0 /)
!    real(8)     :: present_mid_365_trop(nm) = &     ! present-day month middle day in a tropical year
!        (/ 15.5000d0, 45.1211d0, 74.7422d0, 105.2422d0, 135.7422d0, 166.2422d0, 196.7422d0, 227.7422d0, 258.2422d0, &
!            288.7422d0, 319.2422d0, 349.7422d0 /)
!    real(8)     :: present_end_365_trop(nm) = &     ! present-day month ending day in a tropical year
!        (/ 31.0000d0, 59.2422d0, 90.2422d0, 120.2422d0, 151.2422d0, 181.2422d0, 212.2422d0, 243.2422d0, 273.2422d0, &
!            304.2422d0, 334.2422d0, 365.2422d0 /)
    real(8)     :: present_beg_365_progreg(nm) = &  ! present-day month beginning day in a Gregorian year
        (/  0.0000d0, 31.0000d0, 59.2425d0, 90.2425d0, 120.2425d0, 151.2425d0, 181.242d0, 212.2425d0, 243.2425d0, &
            273.2425d0, 304.2425d0, 334.2425d0 /)
    real(8)     :: present_mid_365_progreg(nm) = &  ! present-day month beginning day in a Gregorian year
        (/ 15.5000d0, 45.1213d0, 74.7425d0, 105.2425d0, 135.7425d0, 166.2425d0, 196.7425d0, 227.7425d0, 258.2425d0, &
        288.7425d0, 319.2425d0, 349.7425d0 /)
    real(8)     :: present_end_365_progreg(nm) = &  ! present-day month beginning day in a Gregorian year
        (/ 31.0000d0, 59.2425d0, 90.2425d0, 120.2425d0, 151.2425d0, 181.2425d0, 212.2425d0, 243.2425d0, 273.2425d0, &
        304.2425d0, 334.2425d0, 365.2425d0 /)

    ! calendar type
    character(23), intent(in)   :: calendar_type

    ! simulation age-related variables (controls orbital elements)
    integer(4), intent(in)  :: begageBP                     ! beginning year (BP) (negative, e.g. 10 ka = -10000 BP)
    integer(4), intent(in)  :: agestep                      ! age step size
    integer(4), intent(in)  :: nages                        ! number of simulation ages

    ! individual model simulation year-related variables (controls equinox and solstice days and leap-year status)
    integer(4), intent(in)  :: begyrCE                      ! beginning (pseudo-) year of individual model simulation
    integer(4), intent(in)  :: nsimyrs                      ! number of years of simulation

    ! (output) month-length variables
    integer(4), intent(out) :: iageBP(nages*nsimyrs)        ! year BP 1950 (negative, e.g. 1900 CE = -50.0d0 BP)
    integer(4), intent(out) :: iyearCE(nages*nsimyrs)       ! yearCE simulation year
    integer(4), intent(out) :: imonlen(nages*nsimyrs,nm)    ! integer-value month lengths
    integer(4), intent(out) :: imonmid(nages*nsimyrs,nm)    ! integer-value mid-month days
    integer(4), intent(out) :: imonbeg(nages*nsimyrs,nm)    ! integer-value beginning days
    integer(4), intent(out) :: imonend(nages*nsimyrs,nm)    ! integer-value ending days
    real(8), intent(out)    :: rmonlen(nages*nsimyrs,nm)    ! real-value month lengths
    real(8), intent(out)    :: rmonmid(nages*nsimyrs,nm)    ! real-value mid-month days
    real(8), intent(out)    :: rmonbeg(nages*nsimyrs,nm)    ! real-value month beginning day
    real(8), intent(out)    :: rmonend(nages*nsimyrs,nm)    ! real-value month ending days
    real(8), intent(out)    :: VE_day(nages*nsimyrs)        ! real-value vernal equinox day in simulation year
    real(8), intent(out)    :: SS_day(nages*nsimyrs)        ! real-value (northern) summer solstice in simulation year
    integer(4), intent(out) :: ndays(nages*nsimyrs)         ! integer number of days in simulation year

    ! subroutine GISS_orbpars() and GISS_srevents() input and output arguments
    character(2)            :: year_type = 'BP'             ! AD (AD/BC), CE (CE/BCE), BP (before 1950)
    real(8)                 :: AgeBP                        ! age (BP 1950) (input)
    real(8)                 :: eccen                        ! eccentricity of orbital ellipse
    real(8)                 :: obliq_deg                    ! obliquity (degrees)
    real(8)                 :: perih_deg                    ! longitude of perihelion (degrees)
    real(8)                 :: precc                        ! climatological precession parameter = eccen * sin(omegvp)
    real(8)                 :: veqday                       ! (real) day of vernal equinox

    ! monlen() subroutine arguments
    real(8)                 :: yrlen                        ! real number of days in year (year length)
    integer(4)              :: ndyr                         ! integer number of days in year

    ! present-day month-length variables
    real(8)                 :: rmonlen_0ka(nm)              ! real calculated month lengths at 0 ka
    real(8)                 :: rmonlen_0ka_leap(nm)         ! real calculated month lengths at 0 ka in a leap year
    real(8)                 :: present_monlen(nm)           ! "present day" month lengths
    real(8)                 :: rmonbeg_0ka(nm)              ! real calculated month beginning day at 0ka
    real(8)                 :: rmonmid_0ka(nm)              ! real calculated month mid-month day at 0ka
    real(8)                 :: rmonend_0ka(nm)              ! real calculated month ending day at 0ka
    real(8)                 :: rmonbeg_0ka_leap(nm)         ! real calculated month beginning day at 0ka in a leap year
    real(8)                 :: rmonmid_0ka_leap(nm)         ! real calculated month mid-month day at 0ka in a leap year
    real(8)                 :: rmonend_0ka_leap(nm)         ! real calculated month ending day at 0ka in a leap year

    ! arrays for calculating various month-length statistics
    real(8)                 :: rmonlen_rel(nm)              ! difference between real month lengths and present
    real(8)                 :: rmonlen_ratio(nm)            ! ratio of real month lengths and present
    real(8)                 :: ryeartot(nages*nsimyrs)      ! real total number of days in year
    integer(4)              :: iyeartot(nages*nsimyrs)      ! integer total number of days in year
    
!    integer(4)              :: yearlen_CE
!    integer(4)              :: yearlen_BP
    
    ! indices
    integer(4)              :: n        ! simulation-age index
    integer(4)              :: i        ! simulation-year index
    integer(4)              :: ii       ! age and year index
    integer(4)              :: m        ! month index

    logical                 :: adjust_to_0ka = .true.   ! adjust month-lengths and start days to common values at 0 ka/1950 CE

    ! check for supported calendar types
    select case (trim(calendar_type))
    case ('360_day','noleap','365_day','365.2425','proleptic_gregorian','progreg','gregorian','standard')
        continue
    case default
        stop "Calendar type not supported"
    end select
    
! ===================================================================================================================
! Step 1:  generate target years -- experiment ages (in yrs BP) and simulation years (in yrs CE)

    ii=0
    do n = 1, nages
        do i = 1, nsimyrs
            ii = ii+1
            iageBP(ii) = begageBP + (n - 1) * agestep
            iyearCE(ii) = begyrCE + (i - 1)
        end do
    end do

    ! initialize arrays
    imonlen = 0; imonmid = 0; imonbeg = 0; imonend = 0
    rmonlen = 0.0d0; rmonmid = 0.0d0; rmonbeg = 0.0d0; rmonend = 0.0d0
    VE_day = 0.0d0; ndays = 0
    ryeartot = 0.0d0; iyeartot = 0

! ===================================================================================================================
! Step 2:  get 0 ka (1950CE) calculated month lengths and beginning, middle and end days, which can be used to adjust 
! all other calculated month lengths to nominal "present-day" values

! orbital elements for 0 ka
    ageBP = 0.0d0
    call GISS_orbpars('BP', ageBP, eccen, obliq_deg, perih_deg, precc)

! ===================================================================================================================
! Step 3:  0 ka calculated month lengths, also set subroutine arguments

    select case (trim(calendar_type))
    case ('360_day')
        yrlen = dble(nd_360); ndyr = nd_360; veqday = veqday_360; present_monlen = present_mon_360
        call monlen(yrlen, veqday, int(present_monlen), eccen, perih_deg, rmonlen_0ka, rmonbeg_0ka, rmonmid_0ka, rmonend_0ka)
    case ('noleap', '365_day')
        yrlen = dble(nd_365); ndyr = nd_365; veqday = veqday_365; present_monlen = present_mon_noleap
        call monlen(yrlen, veqday, int(present_monlen), eccen, perih_deg, rmonlen_0ka, rmonbeg_0ka, rmonmid_0ka, rmonend_0ka)
    case ('366_day')
        yrlen = dble(nd_366); ndyr = nd_366; veqday = veqday_366; present_monlen = present_mon_leap
        call monlen(yrlen, veqday,  int(present_monlen), eccen, perih_deg, rmonlen_0ka, rmonbeg_0ka, rmonmid_0ka, rmonend_0ka)
    case ('365.2425')
        yrlen = dble(progreg_year); ndyr = nd_365; veqday = veqday_365; present_monlen = present_mon_365_progreg
        call monlen(yrlen, veqday, int(present_monlen), eccen, perih_deg, rmonlen_0ka, rmonbeg_0ka, rmonmid_0ka, rmonend_0ka)
    case ('proleptic_gregorian','progreg','gregorian','standard')
        ! 365-day year
        yrlen = dble(nd_365); ndyr = nd_365; veqday = veqday_365; present_monlen = present_mon_noleap
        call monlen(yrlen, veqday, int(present_monlen), eccen, perih_deg, rmonlen_0ka, rmonbeg_0ka, rmonmid_0ka, rmonend_0ka)
        ! 366-day year
        yrlen = dble(nd_366); ndyr = nd_366; veqday = veqday_366; present_monlen = present_mon_leap
        call monlen(yrlen, veqday, int(present_monlen), eccen, perih_deg, rmonlen_0ka_leap, & 
            rmonbeg_0ka_leap, rmonmid_0ka_leap, rmonend_0ka_leap)
    case default
        stop "calendar type not supported"
    end select

! loop over simulation ages and years

    ii = 0
    do n = 1, nages
        do i = 1, nsimyrs
            ii = ii + 1

! ===================================================================================================================
! Step 4:  orbital elements for simulation age (e.g. 6 ka)

            call GISS_orbpars('BP', dble(iageBP(n)), eccen, obliq_deg, perih_deg, precc)
            
! ===================================================================================================================
! Steps 5&6:  get real-valued month lengths for different calendars (step 5) and adjust values to 0 ka (step 6)

            select case (trim(calendar_type))
                
            ! non time-varying calendars
            case ('360_day')
                yrlen = dble(nd_360); ndyr = nd_360; veqday = veqday_360
                call monlen(yrlen, veqday, int(present_mon_360), eccen, perih_deg, rmonlen(ii,:), &
                    rmonbeg(ii,:), rmonmid(ii,:), rmonend(ii,:))
                ! adjust values so that 0 ka (1950 CE) will have nominal present-day month lengths
                if (adjust_to_0ka) call adjust_to_ref_length(rmonlen(ii,:), rmonlen_0ka, present_mon_360)
                if (adjust_to_0ka) call adjust_to_ref_day(rmonbeg(ii,:), rmonbeg_0ka, present_beg_360)
                if (adjust_to_0ka) call adjust_to_ref_day(rmonmid(ii,:), rmonmid_0ka, present_mid_360)
                if (adjust_to_0ka) call adjust_to_ref_day(rmonend(ii,:), rmonend_0ka, present_end_360)
            case ('noleap', '365_day')
                yrlen = dble(nd_365); ndyr = nd_365; veqday = veqday_365
                call monlen(yrlen, veqday, int(present_mon_noleap), eccen, perih_deg, rmonlen(ii,:), &
                    rmonbeg(ii,:), rmonmid(ii,:), rmonend(ii,:))
                ! adjust values so that 0 ka (1950 CE) will have nominal present-day month lengths
                if (adjust_to_0ka) call adjust_to_ref_length(rmonlen(ii,:), rmonlen_0ka, present_mon_noleap)
                if (adjust_to_0ka) call adjust_to_ref_day(rmonbeg(ii,:), rmonbeg_0ka, present_beg_365)
                if (adjust_to_0ka) call adjust_to_ref_day(rmonmid(ii,:), rmonmid_0ka, present_mid_365)
                if (adjust_to_0ka) call adjust_to_ref_day(rmonend(ii,:), rmonend_0ka, present_end_365)
            case ('366_day')
                yrlen = dble(nd_366); ndyr = nd_366; veqday = veqday_366
                call monlen(yrlen, veqday, int(present_mon_leap), eccen, perih_deg, rmonlen(ii,:), &
                    rmonbeg(ii,:), rmonmid(ii,:), rmonend(ii,:))
                ! adjust values so that 0 ka (1950 CE) will have nominal present-day month lengths
                if (adjust_to_0ka) call adjust_to_ref_length(rmonlen(ii,:), rmonlen_0ka, present_mon_leap)
                if (adjust_to_0ka) call adjust_to_ref_day(rmonbeg(ii,:), rmonbeg_0ka, present_beg_366)
                if (adjust_to_0ka) call adjust_to_ref_day(rmonmid(ii,:), rmonmid_0ka, present_mid_366)
                if (adjust_to_0ka) call adjust_to_ref_day(rmonend(ii,:), rmonend_0ka, present_end_366)                
            case ('365.2425')
                yrlen = dble(progreg_year); ndyr = nd_365; veqday = veqday_365
                call monlen(yrlen, veqday, int(present_mon_365_progreg), eccen, perih_deg, rmonlen(ii,:), &
                    rmonbeg(ii,:), rmonmid(ii,:), rmonend(ii,:))
                ! adjust values so that 0 ka (1950 CE) will have nominal present-day month lengths
                if (adjust_to_0ka) call adjust_to_ref_length(rmonlen(ii,:), rmonlen_0ka, present_mon_365_progreg)
                if (adjust_to_0ka) call adjust_to_ref_day(rmonbeg(ii,:), rmonbeg_0ka, present_beg_365_progreg)
                if (adjust_to_0ka) call adjust_to_ref_day(rmonmid(ii,:), rmonmid_0ka, present_mid_365_progreg)
                if (adjust_to_0ka) call adjust_to_ref_day(rmonend(ii,:), rmonend_0ka, present_end_365_progreg)                
                
            ! proleptic_gregorian-like calendars
            case ('proleptic_gregorian','progreg','gregorian','standard')
                ! check for leap year
                nd_progreg = yearlen_CE(iyearCE(ii))
                if (nd_progreg .eq. 365) then ! (non leap year)
                    yrlen = dble(nd_365); ndyr = nd_365; veqday = veqday_365
                    ! get real-valued month lengths
                    call monlen(yrlen, veqday, int(present_mon_noleap), eccen, perih_deg, rmonlen(ii,:), &
                        rmonbeg(ii,:), rmonmid(ii,:), rmonend(ii,:))
                    ! adjust values so that 0 ka (1950 CE) will have nominal present-day month lengths
                    if (adjust_to_0ka) call adjust_to_ref_length(rmonlen(ii,:), rmonlen_0ka, present_mon_noleap)
                    if (adjust_to_0ka) call adjust_to_ref_day(rmonbeg(ii,:), rmonbeg_0ka, present_beg_365)
                    if (adjust_to_0ka) call adjust_to_ref_day(rmonmid(ii,:), rmonmid_0ka, present_mid_365)
                    if (adjust_to_0ka) call adjust_to_ref_day(rmonend(ii,:), rmonend_0ka, present_end_365)
                else ! ndprogreg = 366 (leap year)
                    yrlen = dble(nd_366); ndyr = nd_366; veqday = veqday_366
                    call monlen(yrlen, veqday, int(present_mon_leap), eccen, perih_deg, rmonlen(ii,:), &
                        rmonbeg(ii,:), rmonmid(ii,:), rmonend(ii,:))
                    ! adjust values so that 0 ka (1950 CE) will have nominal present-day month lengths
                    if (adjust_to_0ka) call adjust_to_ref_length(rmonlen(ii,:), rmonlen_0ka_leap, present_mon_leap)
                    if (adjust_to_0ka) call adjust_to_ref_day(rmonbeg(ii,:), rmonbeg_0ka_leap, present_beg_366)
                    if (adjust_to_0ka) call adjust_to_ref_day(rmonmid(ii,:), rmonmid_0ka_leap, present_mid_366)
                    if (adjust_to_0ka) call adjust_to_ref_day(rmonend(ii,:), rmonend_0ka_leap, present_end_366)
                end if
            case default ! 
                stop "calendar_type"
            end select

! ===================================================================================================================
! Step 7:  require the sum of month lengths each year to exactly equal the year length

            if (adjust_to_0ka) call adjust_to_yeartot(rmonlen(ii,:), yrlen, ryeartot(ii))

! ===================================================================================================================
! Step 8:  integer month lengths

            call integer_monlen(rmonlen(ii,:), ndyr, imonlen(ii,:), iyeartot(ii))

! ===================================================================================================================
! Step 9: integer beginning, middle and ending days

            call imon_begmidend(imonlen(ii,:), rmonbeg(ii,:), imonbeg(ii,:), imonmid(ii,:), imonend(ii,:))       
            
! ===================================================================================================================
! Step 10: get vernal equinox and summer solstice days for ii-th simulation year

            select case (trim(calendar_type))
            case ('360_day')
                VE_day(ii) = veqday_360; SS_day(ii) = ssday_360; ndays(ii) = nd_360
            case ('noleap', '365_day', '365.2425')
                VE_day(ii) = veqday_365; SS_day(ii) = ssday_365; ndays(ii) = nd_365
            case ('366_day')
                VE_day(ii) = veqday_366; SS_day(ii) = ssday_366; ndays(ii) = nd_366
            case ('proleptic_gregorian','progreg','gregorian','standard')
                ! get vernal equinox and summer solstice days for model simulation year (not age)
                year_type = 'CE'
                call GISS_srevents(year_type, iyearCE(ii), progreg_year, VE_day(ii), SS_day(ii), perihelion, aphelion, ndays(ii))
            case default
                stop "paleo calendar type not defined"
            end select
            
            ! various month-length statistics
            do m=1,nm
                rmonlen_rel(m) = rmonlen(ii,m)-present_monlen(m)
                rmonlen_ratio(m) = rmonlen(ii,m)/present_monlen(m)
            end do
            
        end do
    end do
    
    contains
    
    integer(4) function yearlen_CE(yearCE)
! gets number of days in a CE year -- no year zero or Gregorian-Julian adjustment

    integer(4)  :: yearCE   ! YearCE

    yearlen_CE = 365
    if (mod(yearCE, 4) .eq. 0) yearlen_CE = 366
    if (mod(yearCE, 100) .eq. 0) yearlen_CE = 365
    if (mod(yearCE, 400) .eq. 0) yearlen_CE = 366

    end function yearlen_CE

end subroutine get_month_lengths

!************************************************************************************************************
subroutine monlen(yrlen, veqday, imonlen, eccen, perih, rmonlen, rmonbeg, rmonmid, rmonend)
! calculate month lengths, and beginning, middle and ending day times for a particular orbital configuration, 
! calendar and vernal equinox day using a "traverse-time/time-of-flight" expression based on Kepler's equation:
! (Curtis, H.D. (2014) Orbital Mechanics for Engineering Students, Elsevier, Ch. 3)
!************************************************************************************************************

    implicit none
    
    external kepler_time

    integer(4), parameter   :: nm=12               ! number of months per year
    real(8), intent(in)     :: yrlen                ! total length of year (days)
!    integer(4), intent(in)  :: ndyr                 ! number of days in year
    real(8), intent(in)     :: veqday               ! real vernal equinox day
    integer(4), intent(in)  :: imonlen(nm)          ! number of days in each month at present (calendar dependent)
    real(8), intent(in)     :: eccen, perih         ! orbital parameters

    real(8), intent(out)    :: rmonlen(nm)          ! real month lengths
    real(8), intent(out)    :: rmonbeg(nm)          ! real month beginning day
    real(8), intent(out)    :: rmonmid(nm)          ! real month middle day
    real(8), intent(out)    :: rmonend(nm)          ! real month ending day
    
    real(8)                 :: angle_veq_to_Jan1    ! angle between vernal equinox and Jan 1
    real(8)                 :: tt_perih_to_veq      ! traverse time, perihelion to vernal equinox
    real(8)                 :: angle_perih_to_Jan1  ! angle, perihelion to Jan 1
    real(8)                 :: angle_perih_to_veq   ! angle, perihelion to vernal equinox
    
    real(8)                 :: mon_angle(nm)        ! month angle (degrees)
    real(8)                 :: beg_angle(nm), mid_angle(nm), end_angle(nm) ! beginning, middle and end of month, relative to Jan 1
    
    ! angles and traverse times for beginning, middle and ending days of each month 
    real(8)                 :: perih_angle_month_beg(nm), tt_month_beg(nm), t_month_beg(nm)
    real(8)                 :: perih_angle_month_mid(nm), tt_month_mid(nm), t_month_mid(nm)
    real(8)                 :: perih_angle_month_end(nm), tt_month_end(nm), t_month_end(nm)

    real(8)     :: pi                               ! pi

    integer(4)  :: m                                ! indices

    pi=4.0d0*datan(1.0d0)

    rmonlen = 0.0d0; rmonbeg = 0.0d0; rmonmid = 0.0d0; rmonend = 0.0d0 
    
    ! angle and traverse time, perihelion to vernal equinox 
    angle_perih_to_veq = 360.00d0 - perih
    call kepler_time(eccen, yrlen, angle_perih_to_veq, tt_perih_to_veq)

    ! angle, perihelion to Jan1
    angle_perih_to_Jan1 = angle_perih_to_veq - veqday * (360.0d0/yrlen)

    ! angle, vernal equinox to Jan1 (Jan 1 begins at 0.0 (or 360.0 degrees))
    ! for checking, angle_perih_to_Jan1 * (yrlen/360.0d0) should equal 0.0
    angle_veq_to_Jan1 = 360.0d0 - veqday * (360.0d0/yrlen)

    ! angles, traverse times, month lengths, etc. for individual months
    ! month angle (and beginning, middle and end of month, relative to Jan1)
    mon_angle(1) = dble(imonlen(1) * (360.0d0/yrlen))
    beg_angle(1) = 0.0d0
    mid_angle(1) = beg_angle(1) + mon_angle(1)/2.0d0
    end_angle(1) = mon_angle(1)
    do m=2,nm
        mon_angle(m) = dble(imonlen(m)) * (360.0d0/yrlen)
        beg_angle(m) = beg_angle(m-1) + mon_angle(m-1)
        mid_angle(m) = beg_angle(m) + mon_angle(m) / 2.0d0
        end_angle(m) = beg_angle(m) + mon_angle(m)
    end do

    ! angles from perihelion etc. for month beginning, mid, and ending days
    do m = 1, nm
    
        ! angle from perihelion for each month's beginning, middle and ending days (on circular orbit)
        perih_angle_month_beg(m) = angle_perih_to_Jan1 + beg_angle(m) 
        perih_angle_month_mid(m) = angle_perih_to_Jan1 + mid_angle(m) 
        perih_angle_month_end(m) = angle_perih_to_Jan1 + end_angle(m)  
    
        ! traverse times from perihelion for each month's beginning, middle and ending days (on elliptical orbit)
        call kepler_time(eccen, yrlen, perih_angle_month_beg(m), tt_month_beg(m))
        call kepler_time(eccen, yrlen, perih_angle_month_mid(m), tt_month_mid(m))
        call kepler_time(eccen, yrlen, perih_angle_month_end(m), tt_month_end(m))
    
        ! traverse time from vernal equinox for each month's beginning, middle and ending days (on elliptical orbit)
        t_month_beg(m) = tt_month_beg(m) - tt_perih_to_veq
        t_month_mid(m) = tt_month_mid(m) - tt_perih_to_veq
        t_month_end(m) = tt_month_end(m) - tt_perih_to_veq
    
        ! beginning, middle and ending days of each month (relative to Jan1)
        rmonbeg(m) = dmod(t_month_beg(m) + veqday, yrlen)
        if (perih_angle_month_beg(m) .gt. 360.0d0) rmonbeg(m) = rmonbeg(m) + yrlen
        rmonmid(m) = dmod(t_month_mid(m) + veqday, yrlen)
        if (perih_angle_month_mid(m) .gt. 360.0d0) rmonmid(m) = rmonmid(m) + yrlen
        rmonend(m) = dmod(t_month_end(m) + veqday, yrlen) 
        if (perih_angle_month_end(m) .gt. 360.0d0) rmonend(m) = rmonend(m) + yrlen
    
    end do
    
    ! fixup for ending day of last month (when end of last month appers in beginning of year)
    if (rmonend(nm) .lt. 30.0d0) rmonend(nm) = rmonend(nm) + yrlen  

    ! month length (month beginning to next month beginning)
    do m = 1,nm-1
        rmonlen(m) = t_month_beg(m+1) - t_month_beg(m)
        if (rmonlen(m) .le. 0.0d0) rmonlen(m) = rmonlen(m) + yrlen
    end do
    rmonlen(nm) = t_month_beg(1) - t_month_beg(nm)
    if (rmonlen(nm) .le. 0.0d0) rmonlen(nm) = rmonlen(nm) + yrlen

end subroutine monlen

!************************************************************************************************************    
subroutine kepler_time(eccen, T, theta_deg, time)
! travel time along orbit relative to periapsis/perihelion (i.e. 0.0 at perihelion)
! input:  true anomaly (theta, degrees); output:  traverse time since perihelion (time, same units as T)
! (Curtis, H.D. (2014) Orbital Mechanics for Engineering Students, Elsevier, Ch. 3)
!************************************************************************************************************

    implicit none
    
    real(8), intent(in)     :: eccen
    real(8), intent(in)     :: T                ! orbital period (yrlen)
    real(8), intent(in)     :: theta_deg        ! "true" anomaly (angle from perihelion (degrees)
    real(8), intent(out)    :: time             ! traverse time from periapsis (e.g. perihelion)
    
!    real(8)                 :: radians
    real(8)                 :: theta_rad        ! theta_rad (radians)
    real(8)                 :: M                ! mean anomaly
    real(8)                 :: E                ! eccentric anomaly
    real(8)                 :: pi
    
    pi=4.0d0*datan(1.0d0)  
    
    theta_rad = radians(theta_deg)
    E = 2.0d0 * datan( dsqrt((1.0d0 - eccen) / (1.0d0 + eccen)) * dtan(theta_rad/ 2.0d0) )  ! eq 3.13b
    M = E - eccen*dsin(E)                                                                   ! eq 3.14
    time = (M / (2.0d0 * pi)) * T                                                           ! eq 3.15
    if (time .lt. 0.0d0) time = time + T

    contains
    real(8) function radians(d)
! decimal degrees to radians

    real(8) d, pi
    pi=4.0d0*datan(1.0d0)
    radians=d*((2.0d0*pi)/360.0d0)

    end function radians
end subroutine kepler_time    

!************************************************************************************************************
subroutine adjust_to_ref_length(rmonlen, rmonlenref, rmonlentarg)
!************************************************************************************************************

    implicit none

    integer(4), parameter   :: nm=12               ! number of months per year
    real(8), intent(inout)  :: rmonlen(nm)      ! real month lengths
    real(8), intent(in)     :: rmonlenref(nm)   ! reference month lengths (usually calculated 0 ka)
    real(8), intent(in)     :: rmonlentarg(nm)  ! target month lengths (usually conventional 0 ka)
    
    real(8)                 :: rmonlentot_in, rmonlentot_out

    integer(4)              :: m

    rmonlentot_in = 0.0d0; rmonlentot_out = 0.0d0
    ! adjust rmonlen to reference year
    do m=1,nm
        rmonlentot_in = rmonlentot_in + rmonlen(m)
        rmonlen(m) = rmonlen(m) * (rmonlentarg(m)/rmonlenref(m))
        rmonlentot_out = rmonlentot_out + rmonlen(m)
    end do
    
end subroutine adjust_to_ref_length

!************************************************************************************************************
subroutine adjust_to_ref_day(rmonday, rdayref, rdaytarg)
!************************************************************************************************************

    implicit none

    integer(4), parameter   :: nm=12               ! number of months per year
    real(8), intent(inout)  :: rmonday(nm)      ! real beginning, middle or ending days
    real(8), intent(in)     :: rdayref(nm)      ! reference day (usually calculated 0 ka)
    real(8), intent(in)     :: rdaytarg(nm)     ! target day (usually conventional 0 ka)

    integer(4)              :: m

    ! adjust rmonday to reference year 
    do m=1,nm
        rmonday(m) = rmonday(m) - (rdayref(m) - rdaytarg(m))
    end do

end subroutine adjust_to_ref_day
!************************************************************************************************************
subroutine adjust_to_yeartot(rmonlen, ryeartottarg, ryeartot)
!************************************************************************************************************
    implicit none

    integer(4), parameter   :: nm = 12              ! number of months in the year
    real(8), intent(inout)  :: rmonlen(nm)      ! real month lengths
    real(8), intent(in)     :: ryeartottarg     ! total annual month lengths target
    real(8), intent(out)    :: ryeartot         ! total annual month lengths

    integer(4)              :: m

    ! get ryeartot
    ryeartot=0.0d0
    do m=1,nm
        ryeartot = ryeartot + rmonlen(m)
    end do

    ! adjust rmonlen to ryeartot target
    do m=1,nm
        rmonlen(m) = rmonlen(m) * (ryeartottarg/ryeartot)
    end do

    ! recalc ryeartot
    ryeartot=0.0d0
    do m=1,nm
        ryeartot = ryeartot + rmonlen(m)
    end do

end subroutine adjust_to_yeartot

!************************************************************************************************************
subroutine integer_monlen(rmonlen, ndtarg, imonlen, iyeartot)
!************************************************************************************************************
    implicit none

    integer(4), parameter   :: nm = 12              ! number of months in the year
    real(8), intent(in)     :: rmonlen(nm)      ! real month lengths
    integer(4), intent(in)  :: ndtarg           ! total annual integer number of days target
    integer(4), intent(out) :: imonlen(nm)      ! integer month lengths
    integer(4), intent(out) :: iyeartot

    integer(4)              :: diff_sign
    real(8)                 :: ryeartot
    real(8)                 :: inc
    integer(4)              :: i,m

    ! integer month lengths
    iyeartot=0; ryeartot=0.0d0
    do  m=1,nm
        imonlen(m)=idint(dnint(rmonlen(m)))
        iyeartot=iyeartot+imonlen(m)
        ryeartot=ryeartot+rmonlen(m)
    end do

    ! force integer month lengths to sum to ndtarg
    if (iyeartot.ne.ndtarg) then

        ! add (diff_sign = 1) or subtract (diff_sign = -1) a small increment to each rmonlen value
        ! iterate until iyeartot = ndtarg, incrementing small increment via i (i*inc)
        diff_sign=1
        if ((iyeartot-ndtarg) .gt. 0) diff_sign=-1
        inc=0.000001; i=0
        do
            if (iyeartot.eq.ndtarg) exit
            i=i+1
            iyeartot=0
            do m=1,nm
                imonlen(m)=idint(dnint(rmonlen(m)+i*inc*diff_sign))
                iyeartot=iyeartot+imonlen(m)
            end do

            diff_sign=1
            if ((iyeartot-ndtarg) .gt. 0) diff_sign=-1

        end do

    end if

end subroutine integer_monlen
!************************************************************************************************************
subroutine imon_begmidend(imonlen, rmonbeg, imonbeg, imonmid, imonend)
!************************************************************************************************************
    implicit none

    integer(4), parameter   :: nm = 12              ! number of months in the year
    integer(4), intent(in)  :: imonlen(nm)  ! integer month length
    real(8), intent(in)     :: rmonbeg(nm)  ! real month beginning day
    integer(4), intent(out) :: imonbeg(nm)  ! integer beginning day of month
    integer(4), intent(out) :: imonmid(nm)  ! integer mid-month day
    integer(4), intent(out) :: imonend(nm)  ! integer ending day of month

    integer(4)              :: m

    imonbeg(1) = idint(rmonbeg(1))
    imonend(1) = imonbeg(1) + imonlen(1) - 1
    do m=2,nm
        imonbeg(m) = imonend(m-1) + 1
        imonend(m) = imonend(m-1) + imonlen(m)
    end do
    do m=1,nm
        imonmid(m) = imonbeg(m) + imonlen(m)/2
    end do

end subroutine imon_begmidend
!************************************************************************************************************
subroutine GISS_orbpars(year_type, year, eccen, obliq_deg, perih_deg, precc)
!************************************************************************************************************
    implicit none
    
    external orbpar
    
    ! past ages are negative, e.g. 21 ka = 21,000 cal yr BP = -21000, and 1950 CE = 0 cal yr BP = 0 here
    ! NOTE:  Year CE/AD = 0 is assumed to exist, and is equivalent to 1950 BP (-1950)
    
    character(2), intent(in)    :: year_type    ! year type (AD/BC, CE/BCE, BP (1950 CE))
    real(8), intent(in)         :: year        ! Year (negative, e.g. 1900 CE = -50 BP)
    real(8), intent(out)        :: eccen        ! eccentricity of orbital ellipse 
    real(8), intent(out)        :: obliq_deg    ! obliquity (degrees)
    real(8), intent(out)        :: perih_deg    ! longitude of perihelion (degrees)
    real(8), intent(out)        :: precc        ! climatological precession parameter = eccen * sin(omegvp)
 
    real(8)                     :: obliq        ! obliquity (latitude of Tropic of Cancer) (radians)
    real(8)                     :: omegvp       ! longitude of perihelion = spatial angle from vernal equinox to perihelion (radians)

    real(8)                     :: YearCE       ! YearCE -- for consistency with old code (input to ORBPAR() is YearCE)
    real(8)                     :: YearBP       ! YearBP
    real(8)                     :: pi           ! pi
    
    pi=4.0d0*datan(1.0d0)

    ! NOTE:  Year CE/AD = 0 is assumed to exist, and is equivalent to 1950 BP (-1950)
    ! subroutine orbpar() expects real-valued YearCE, but converts to YearBP for calculations
    select case (year_type)
    case ('CE', 'AD', 'ce', 'ad')
        YearCE = year
        YearBP = year - 1950.0d0
    case ('BP', 'bp')
        YearCE = year + 1950.0d0
        YearBP = year
    case default
        stop 'year_type'
    end select
    
    call orbpar(yearCE, eccen, obliq, omegvp)
    
    obliq_deg=obliq/((2.0d0*pi)/360.0d0)
    precc=eccen*dsin(omegvp)
    perih_deg=omegvp/((2.0d0*pi)/360.0d0)

end subroutine GISS_orbpars
!************************************************************************************************************    
subroutine orbpar(year, eccen, obliq, omegvp)   ! year should be YearCE
!************************************************************************************************************
! NOTE:  Year CE/AD = 0 is assumed to not exist
    
! .f90 version of the ORBPAR() subroutine contained in
!  https://data.giss.nasa.gov/ar5/SOLAR/ORBPAR.FOR downloaded 2017-09-04 17:17
!  also retrievable from https://web.archive.org/web/20150920211936/http://data.giss.nasa.gov/ar5/solar.html
! 
!  ORBPAR calculates the three orbital parameters as a function of
!  YEAR.  The source of these calculations is: Andre L. Berger,
!  1978, "Long-Term Variations of Daily Insolation and Quaternary
!  Climatic Changes", JAS, v.35, p.2362.  Also useful is: Andre L.
!  Berger, May 1978, "A Simple Algorithm to Compute Long Term
!  Variations of Daily Insolation", published by the Institut
!  d'Astronomie et de Geophysique, Universite Catholique de Louvain,
!  Louvain-la Neuve, No. 18.
! 
!  Tables and equations refer to the first reference (JAS).  The
!  corresponding table or equation in the second reference is
!  enclosed in parentheses.  The coefficients used in this
!  subroutine are slightly more precise than those used in either
!  of the references.  The generated orbital parameters are precise
!  within plus or minus 1000000 years from present.
! 
!  Input:  YEAR   = years A.D. are positive, B.C. are negative
!  Output: ECCEN  = eccentricity of orbital ellipse
!          OBLIQ  = latitude of Tropic of Cancer in radians
!          OMEGVP = longitude of perihelion 
!                 = spatial angle from vernal equinox to perihelion
!                   in radians with sun as angle vertex
! 
    implicit none
    
    real(8), intent(in)     :: year     ! Year CE 
    real(8), intent(out)    :: eccen, obliq, omegvp
    
    real(8)                 :: table1(3,47),table4(3,19),table5(3,78)
    real(8)                 :: pi, twopi, piz180
    real(8)                 :: ym1950, sumc, arg, obliqd, esinpi, ecospi, pie, fsinfd, psi
    
    integer(4)              :: i
    
! Table 1 (2).  Obliquity relative to mean ecliptic of date: OBLIQD
    data table1 / &
        -2462.2214466d0, 31.609974d0, 251.9025d0, &
         -857.3232075d0, 32.620504d0, 280.8325d0, &
         -629.3231835d0, 24.172203d0, 128.3057d0, &
         -414.2804924d0, 31.983787d0, 292.7252d0, &
         -311.7632587d0, 44.828336d0,  15.3747d0, &
          308.9408604d0, 30.973257d0, 263.7951d0, &
         -162.5533601d0, 43.668246d0, 308.4258d0, &
         -116.1077911d0, 32.246691d0, 240.0099d0, &
          101.1189923d0, 30.599444d0, 222.9725d0, &
          -67.6856209d0, 42.681324d0, 268.7809d0, &
           24.9079067d0, 43.836462d0, 316.7998d0, &
           22.5811241d0, 47.439436d0, 319.6024d0, &
          -21.1648355d0, 63.219948d0, 143.8050d0, &
          -15.6549876d0, 64.230478d0, 172.7351d0, &
           15.3936813d0,  1.010530d0,  28.9300d0, &
           14.6660938d0,  7.437771d0, 123.5968d0, &
          -11.7273029d0, 55.782177d0,  20.2082d0, &
           10.2742696d0,  0.373813d0,  40.8226d0, &
            6.4914588d0, 13.218362d0, 123.4722d0, &
            5.8539148d0, 62.583231d0, 155.6977d0, &
           -5.4872205d0, 63.593761d0, 184.6277d0, &
           -5.4290191d0, 76.438310d0, 267.2772d0, &
            5.1609570d0, 45.815258d0,  55.0196d0, &
            5.0786314d0,  8.448301d0, 152.5268d0, &
           -4.0735782d0, 56.792707d0,  49.1382d0, &
            3.7227167d0, 49.747842d0, 204.6609d0, &
            3.3971932d0, 12.058272d0,  56.5233d0, &
           -2.8347004d0, 75.278220d0, 200.3284d0, &
           -2.6550721d0, 65.241008d0, 201.6651d0, &
           -2.5717867d0, 64.604291d0, 213.5577d0, &
           -2.4712188d0,  1.647247d0,  17.0374d0, &
            2.4625410d0,  7.811584d0, 164.4194d0, &
            2.2464112d0, 12.207832d0,  94.5422d0, &
           -2.0755511d0, 63.856665d0, 131.9124d0, &
           -1.9713669d0, 56.155990d0,  61.0309d0, &
           -1.8813061d0, 77.448840d0, 296.2073d0, &
           -1.8468785d0,  6.801054d0, 135.4894d0, &
            1.8186742d0, 62.209418d0, 114.8750d0, &
            1.7601888d0, 20.656133d0, 247.0691d0, &
           -1.5428851d0, 48.344406d0, 256.6114d0, &
            1.4738838d0, 55.145460d0,  32.1008d0, &
           -1.4593669d0, 69.000539d0, 143.6804d0, &
            1.4192259d0, 11.071350d0,  16.8784d0, &
           -1.1818980d0, 74.291298d0, 160.6835d0, &
            1.1756474d0, 11.047742d0,  27.5932d0, &
           -1.1316126d0,  0.636717d0, 348.1074d0, &
            1.0896928d0, 12.844549d0,  82.6496d0/
    
! Table 4 (1).  Fundamental elements of the ecliptic: ECCEN sin(pi)
    data table4/ &
          .01860798d0,  4.207205d0,  28.620089d0, &
          .01627522d0,  7.346091d0, 193.788772d0, &
         -.01300660d0, 17.857263d0, 308.307024d0, &
          .00988829d0, 17.220546d0, 320.199637d0, &
         -.00336700d0, 16.846733d0, 279.376984d0, &
          .00333077d0,  5.199079d0,  87.195000d0, &
         -.00235400d0, 18.231076d0, 349.129677d0, &
          .00140015d0, 26.216758d0, 128.443387d0, &
          .00100700d0,  6.359169d0, 154.143880d0, &
          .00085700d0, 16.210016d0, 291.269597d0, &
          .00064990d0,  3.065181d0, 114.860583d0, &
          .00059900d0, 16.583829d0, 332.092251d0, &
          .00037800d0, 18.493980d0, 296.414411d0, &
         -.00033700d0,  6.190953d0, 145.769910d0, &
          .00027600d0, 18.867793d0, 337.237063d0, &
          .00018200d0, 17.425567d0, 152.092288d0, &
         -.00017400d0,  6.186001d0, 126.839891d0, &
         -.00012400d0, 18.417441d0, 210.667199d0, &
          .00001250d0,  0.667863d0,  72.108838d0/
    
!  Table 5 (3).  General precession in longitude: psi
    DATA TABLE5/ &
         7391.0225890d0, 31.609974d0, 251.9025d0, &
         2555.1526947d0, 32.620504d0, 280.8325d0, &
         2022.7629188d0, 24.172203d0, 128.3057d0, &
        -1973.6517951d0,  0.636717d0, 348.1074d0, &
         1240.2321818d0, 31.983787d0, 292.7252d0, &
          953.8679112d0,  3.138886d0, 165.1686d0, &
         -931.7537108d0, 30.973257d0, 263.7951d0, &
          872.3795383d0, 44.828336d0,  15.3747d0, &
          606.3544732d0,  0.991874d0,  58.5749d0, &
         -496.0274038d0,  0.373813d0,  40.8226d0, &
          456.9608039d0, 43.668246d0, 308.4258d0, &
          346.9462320d0, 32.246691d0, 240.0099d0, &
         -305.8412902d0, 30.599444d0, 222.9725d0, &
          249.6173246d0,  2.147012d0, 106.5937d0, &
         -199.1027200d0, 10.511172d0, 114.5182d0, &
          191.0560889d0, 42.681324d0, 268.7809d0, &
         -175.2936572d0, 13.650058d0, 279.6869d0, &
          165.9068833d0,  0.986922d0,  39.6448d0, &
          161.1285917d0,  9.874455d0, 126.4108d0, &
          139.7878093d0, 13.013341d0, 291.5795d0, &
         -133.5228399d0,  0.262904d0, 307.2848d0, &
          117.0673811d0,  0.004952d0,  18.9300d0, &
          104.6907281d0,  1.142024d0, 273.7596d0, &
           95.3227476d0, 63.219948d0, 143.8050d0, &
           86.7824524d0,  0.205021d0, 191.8927d0, &
           86.0857729d0,  2.151964d0, 125.5237d0, &
           70.5893698d0, 64.230478d0, 172.7351d0, &
          -69.9719343d0, 43.836462d0, 316.7998d0, &
          -62.5817473d0, 47.439436d0, 319.6024d0, &
           61.5450059d0,  1.384343d0,  69.7526d0, &
          -57.9364011d0,  7.437771d0, 123.5968d0, &
           57.1899832d0, 18.829299d0, 217.6432d0, &
          -57.0236109d0,  9.500642d0,  85.5882d0, &
          -54.2119253d0,  0.431696d0, 156.2147d0, &
           53.2834147d0,  1.160090d0,  66.9489d0, &
           52.1223575d0, 55.782177d0,  20.2082d0, &
          -49.0059908d0, 12.639528d0, 250.7568d0, &
          -48.3118757d0,  1.155138d0,  48.0188d0, &
          -45.4191685d0,  0.168216d0,   8.3739d0, &
          -42.2357920d0,  1.647247d0,  17.0374d0, &
          -34.7971099d0, 10.884985d0, 155.3409d0, &
           34.4623613d0,  5.610937d0,  94.1709d0, &
          -33.8356643d0, 12.658184d0, 221.1120d0, &
           33.6689362d0,  1.010530d0,  28.9300d0, &
          -31.2521586d0,  1.983748d0, 117.1498d0, &
          -30.8798701d0, 14.023871d0, 320.5095d0, &
           28.4640769d0,  0.560178d0, 262.3602d0, &
          -27.1960802d0,  1.273434d0, 336.2148d0, &
           27.0860736d0, 12.021467d0, 233.0046d0, &
          -26.3437456d0, 62.583231d0, 155.6977d0, &
           24.7253740d0, 63.593761d0, 184.6277d0, &
           24.6732126d0, 76.438310d0, 267.2772d0, &
           24.4272733d0,  4.280910d0,  78.9281d0, &
           24.0127327d0, 13.218362d0, 123.4722d0, &
           21.7150294d0, 17.818769d0, 188.7132d0, &
          -21.5375347d0,  8.359495d0, 180.1364d0, &
           18.1148363d0, 56.792707d0,  49.1382d0, &
          -16.9603104d0,  8.448301d0, 152.5268d0, &
          -16.1765215d0,  1.978796d0,  98.2198d0, &
           15.5567653d0,  8.863925d0,  97.4808d0, &
           15.4846529d0,  0.186365d0, 221.5376d0, &
           15.2150632d0,  8.996212d0, 168.2438d0, &
           14.5047426d0,  6.771027d0, 161.1199d0, &
          -14.3873316d0, 45.815258d0,  55.0196d0, &
           13.1351419d0, 12.002811d0, 262.6495d0, &
           12.8776311d0, 75.278220d0, 200.3284d0, &
           11.9867234d0, 65.241008d0, 201.6651d0, &
           11.9385578d0, 18.870667d0, 294.6547d0, &
           11.7030822d0, 22.009553d0,  99.8233d0, &
           11.6018181d0, 64.604291d0, 213.5577d0, &
          -11.2617293d0, 11.498094d0, 154.1631d0, &
          -10.4664199d0,  0.578834d0, 232.7153d0, &
           10.4333970d0,  9.237738d0, 138.3034d0, &
          -10.2377466d0, 49.747842d0, 204.6609d0, &
           10.1934446d0,  2.147012d0, 106.5938d0, &
          -10.1280191d0,  1.196895d0, 250.4676d0, &
           10.0289441d0,  2.133898d0, 332.3345d0, &
          -10.0034259d0,  0.173168d0,  27.3039d0/
!
    pi = 4.0d0*datan(1.0d0)
    twopi = 2.0d0*pi
    piz180 = twopi/360.0d0
    
    ! NOTE:  change Year CE to Year BP for calculations
    ym1950 = year-1950.0d0

    !  Obliquity from Table 1 (2):
    !    OBLIQ# = 23.320556 (degrees)             Equation 5.5 (15)
    !    OBLIQD = OBLIQ# + sum[A cos(ft+delta)]   Equation 1 (5)
    
    sumc = 0.
    do  i=1,47
        arg = piz180*(ym1950*table1(2,i)/3600.0d0+table1(3,i))
        sumc = sumc + table1(1,i)*dcos(arg)
    end do
    obliqd = 23.320556d0 + sumc/3600.0d0
    obliq = obliqd*piz180

    !  Eccentricity from Table 4 (1):
    !    ECCEN sin(pi) = sum[M sin(gt+beta)]           Equation 4 (1)
    !    ECCEN cos(pi) = sum[M cos(gt+beta)]           Equation 4 (1)
    !    ECCEN = ECCEN sqrt[sin(pi)^2 + cos(pi)^2]
    
    esinpi = 0.d0
    ecospi = 0.d0
    do i=1,19
        arg = piz180*(ym1950*table4(2,i)/3600.0d0+table4(3,i))
        esinpi = esinpi + table4(1,i)*dsin(arg)
        ecospi = ecospi + table4(1,i)*dcos(arg)
    end do
    eccen = sqrt(esinpi*esinpi+ecospi*ecospi)

    !  Perihelion from Equation 4,6,7 (9) and Table 4,5 (1,3):
    !    PSI# = 50.439273 (seconds of degree)         Equation 7.5 (16)
    !    ZETA =  3.392506 (degrees)                   Equation 7.5 (17)
    !    PSI = PSI# t + ZETA + sum[F sin(ft+delta)]   Equation 7 (9)
    !    PIE = atan[ECCEN sin(pi) / ECCEN cos(pi)]
    !    OMEGVP = PIE + PSI + 3.14159                 Equation 6 (4.5)
  
    pie = atan2(esinpi,ecospi)
    fsinfd = 0.d0
    do i=1,78
        arg = piz180*(ym1950*table5(2,i)/3600.+table5(3,i))
        fsinfd = fsinfd + table5(1,i)*dsin(arg)
    end do
    psi = piz180*(3.392506d0+(ym1950*50.439273d0+fsinfd)/3600.d0)
    omegvp = modulo(pie+psi+0.5d0*twopi, twopi)

end subroutine orbpar
!*****************************************************************************************************************
subroutine GISS_srevents(year_type, iyear, EDAYzY, veqday, ssday, perihelion, aphelion, ndays_in_year)
! subroutines based on
! SREVENTS.FOR    Solar EVENTS each year    2012/05/29
! https://data.giss.nasa.gov/ar5/SOLAR/SREVENTS.FOR downloaded 2017-09-12
!*****************************************************************************************************************

    implicit none

    external orbpar, DtoYMDHM

    integer(4), parameter       :: nm=12
    character(2), intent(in)    :: year_type
    integer(4), intent(in)      :: iyear
    real(8), intent(in)         :: EDAYzY
    real(8), intent(out)        :: veqday, ssday, perihelion, aphelion
    integer(4), intent(out)     :: ndays_in_year

    integer(4)                  :: YearCE, YearBP
        ! components of event dates
    integer(4)      :: JVEYR,JVEMON,JVEDAT,JVEHR,JVEMIN, JSSYR,JSSMON,JSSDAT,JSSHR,JSSMIN
    integer(4)      :: JAEYR,JAEMON,JAEDAT,JAEHR,JAEMIN, JWSYR,JWSMON,JWSDAT,JWSHR,JWSMIN
    integer(4)      :: JPRYR,JPRMON,JPRDAT,JPRHR,JPRMIN, JAPYR,JAPMON,JAPDAT,JAPHR,JAPMIN
    integer(4)      :: JEXYR,JEXMON,JEXDAT,JEXHR,JEXMIN
    integer(4)      :: KPERIH, KAPHEL

    real(8)         :: vereqx
!    real(8)         :: vernal
!    logical(4)      :: isleap

!  The unit (days) means days measured since 2000 January 1, hour 0

    real(8)         :: year, eccen, obliq, omegvp
    real(8)         :: bsemi
    real(8)         :: TAofVE, EAofVE, MAofVE
    real(8)         :: PERIH1, PERIH2, APHEL1, APHEL2
    real(8)         :: TAofSS, EAofSS, MAofSS
    real(8)         :: SUMSOL
    real(8)         :: TAofAE, EAofAE, MAofAE
    real(8)         :: AUTEQX
    real(8)         :: TAofWS, EAofWS, MAofWS
    real(8)         :: WINSOL

    real(8)         :: pi, twopi

    real(8)     :: nd_noleap(nm) = &       ! present-day month lengths in 365-day (noleap) year
        (/ 31.0d0, 28.0d0, 31.0d0, 30.0d0, 31.0d0, 30.0d0, 31.0d0, 31.0d0, 30.0d0, 31.0d0, 30.0d0, 31.0d0 /)
    real(8)     :: nd_leap(nm) = &         ! present-day month lengths in 366-day (leap) year
        (/ 31.0d0, 29.0d0, 31.0d0, 30.0d0, 31.0d0, 30.0d0, 31.0d0, 31.0d0, 30.0d0, 31.0d0, 30.0d0, 31.0d0 /)
    real(8)     :: accumday_noleap(nm) = &       ! present-day accumulated day lengths in 365-day (noleap) year
        (/  0.0d0, 31.0d0, 59.0d0, 90.0d0,120.0d0,151.0d0,181.0d0,212.0d0,243.0d0,273.0d0,304.0d0,334.0d0 /)
    real(8)     :: accumday_leap(nm) = &         ! present-day accumulated day lengths in 366-day (leap) year
        (/  0.0d0, 31.0d0, 60.0d0, 91.0d0,121.0d0,152.0d0,182.0d0,213.0d0,244.0d0,274.0d0,305.0d0,335.0d0 /)


    pi=4.0d0*datan(1.0d0)
    twopi = pi * 2.0d0

    ! NOTE:  Year CE/AD = 0 is assumed to exist, and is equivalent to 1950 BP (-1950)
    ! subroutine orbpar() expects real-valued YearCE, but converts to YearBP for calculations
    select case (year_type)
    case ('CE', 'AD', 'ce', 'ad')
        YearCE = iyear
        YearBP = iyear - 1950
    case ('BP', 'bp')
        YearCE = iyear + 1950
        YearBP = iyear
    case default
        stop 'year_type'
    end select

!  Determine orbital parameters
    YEAR = dble(YearCE)
    CALL ORBPAR (YEAR, ECCEN,OBLIQ,OMEGVP) ! orbpar() expects YearCE input
    BSEMI  = dSQRT (1.0d0-ECCEN*ECCEN)
!  Vernal Equinox
    VEREQX = VERNAL (YearCE, EDAYzY)! NOTE:  made EDAYzY an argument
    !VEREQX = dmod(VERNAL (YearCE, EDAYzY), EDAYzy) ! NOTE:  made EDAYzY an argument
    CALL DtoYMDHM (VEREQX, JVEYR,JVEMON,JVEDAT,JVEHR,JVEMIN)
    TAofVE = - OMEGVP
    EAofVE = dATAN2 (dSIN(TAofVE)*BSEMI, dCOS(TAofVE)+ECCEN)
    MAofVE = EAofVE - ECCEN*dSIN(EAofVE)
    IF(MAofVE.lt.0.0d0)  MAofVE = MAofVE + TWOPI
!  Perihelion
    KPERIH = 0
    PERIH1 = VEREQX - MAofVE*EDAYzY/TWOPI
    PERIH2 = PERIH1 + EDAYzY
    Call DtoYMDHM (PERIH1, JPRYR,JPRMON,JPRDAT,JPRHR,JPRMIN)
    If (JPRYR /= IYEAR)  GoTo 210
    KPERIH = 1
    Call DtoYMDHM (PERIH2, JEXYR,JEXMON,JEXDAT,JEXHR,JEXMIN)
    If (JEXYR == IYEAR)  KPERIH = 2
    GoTo 220
210 Call DtoYMDHM (PERIH2, JPRYR,JPRMON,JPRDAT,JPRHR,JPRMIN)
    If (JPRYR == IYEAR)  KPERIH = 1
!  Aphelion
220 KAPHEL = 0
    APHEL1 = PERIH2 - 0.5d0*EDAYzY
    APHEL2 = PERIH2 + 0.5d0*EDAYzY
    Call DtoYMDHM (APHEL1, JAPYR,JAPMON,JAPDAT,JAPHR,JAPMIN)
    If (JAPYR /= IYEAR)  GoTo 230
    KAPHEL = 1
    If (KPERIH == 2)  GoTo 240
    Call DtoYMDHM (APHEL2, JEXYR,JEXMON,JEXDAT,JEXHR,JEXMIN)
    If (JEXYR == IYEAR)  KAPHEL = 2
    GoTo 240
230 Call DtoYMDHM (APHEL2, JAPYR,JAPMON,JAPDAT,JAPHR,JAPMIN)
    If (JAPYR == IYEAR)  KAPHEL = 1
!  Summer Solstice
240 TAofSS = TAofVE + 0.25d0*TWOPI
    EAofSS = dATan2 (dSin(TAofSS)*BSEMI, dCos(TAofSS)+ECCEN)
    MAofSS = EAofSS - ECCEN*dSin(EAofSS)
    If (MAofSS-MAofVE < 0)  MAofSS = MAofSS + TWOPI
    SUMSOL = VEREQX + (MAofSS-MAofVE)*EDAYzY/TWOPI
    Call DtoYMDHM (SUMSOL, JSSYR,JSSMON,JSSDAT,JSSHR,JSSMIN)
!  Autumnal Equinox
    TAofAE = TAofVE + 0.5d0*TWOPI
    EAofAE = dATan2 (dSin(TAofAE)*BSEMI, dCos(TAofAE)+ECCEN)
    MAofAE = EAofAE - ECCEN*dSin(EAofAE)
    If (MAofAE-MAofVE < 0)  MAofAE = MAofAE + TWOPI
    AUTEQX = VEREQX + (MAofAE-MAofVE)*EDAYzY/TWOPI
    Call DtoYMDHM (AUTEQX, JAEYR,JAEMON,JAEDAT,JAEHR,JAEMIN)
!  Winter Solstice
    TAofWS = TAofVE + 0.75d0*TWOPI
    EAofWS =dATan2 (dSin(TAofWS)*BSEMI, dCos(TAofWS)+ECCEN)
    MAofWS = EAofWS - ECCEN*dSin(EAofWS)
    If (MAofWS-MAofVE < 0)  MAofWS = MAofWS + TWOPI
    WINSOL = VEREQX + (MAofWS-MAofVE)*EDAYzY/TWOPI
    Call DtoYMDHM (WINSOL, JWSYR,JWSMON,JWSDAT,JWSHR,JWSMIN)

! vernal equinox and northern summer solstice days

    veqday = 0.0d0; ssday = 0.0d0
    if(isleap(YearCE)) then
        ndays_in_year = 366
        veqday = veqday + nd_leap(1) + nd_leap(2) + dble(JVEDAT) + dble(JVEHR)/24.0d0 + dble(JVEMIN)/1440.0d0
        ssday = ssday + nd_leap(1) + nd_leap(2) + nd_leap(3) + nd_leap(4) + nd_leap(5) &
            + dble(JSSDAT) + dble(JSSHR)/24.0d0 + dble(JSSMIN)/1440.0d0
        perihelion = dble(accumday_leap(JPRMON)) + dble(JPRDAT) + dble(JPRHR)/24.0d0 + dble(JPRMIN)/1440.0d0
        aphelion   = dble(accumday_leap(JAPMON)) + dble(JAPDAT) + dble(JAPHR)/24.0d0 + dble(JAPMIN)/1440.0d0
    else
        ndays_in_year = 365
        veqday = veqday + nd_noleap(1) + nd_noleap(2) + dble(JVEDAT) + dble(JVEHR)/24.0d0 + dble(JVEMIN)/1440.0d0
        ssday = ssday + nd_noleap(1) + nd_noleap(2) + nd_noleap(3) + nd_noleap(4) + nd_noleap(5) &
            + dble(JSSDAT) + dble(JSSHR)/24.0d0 + dble(JSSMIN)/1440.0d0
        perihelion = dble(accumday_noleap(JPRMON)) + dble(JPRDAT) + dble(JPRHR)/24.0d0 + dble(JPRMIN)/1440.0d0
        aphelion   = dble(accumday_noleap(JAPMON)) + dble(JAPDAT) + dble(JAPHR)/24.0d0 + dble(JAPMIN)/1440.0d0
    end if
    
    contains
    logical(4) function isleap(yearCE)
! is yearCE a leap year? -- no year zero or Gregorian-Julian adjustment
! NOTE:  Year CE/AD = 0 is assumed to exist, and is equivalent to 1950 BP (-1950)

    integer(4), intent(in)  :: yearCE

    isleap = .false.
    if (mod(yearCE, 4) .eq. 0) isleap = .true.
    if (mod(yearCE, 100) .eq. 0) isleap = .false.
    if (mod(yearCE, 400) .eq. 0) isleap = .true.

    end function isleap
    
    real(8) function vernal (iyear, edayzy)
! subroutines based on
! SREVENTS.FOR    Solar EVENTS each year    2012/05/29
! https://data.giss.nasa.gov/ar5/SOLAR/SREVENTS.FOR downloaded 2017-09-12

!  For a given year, vernal calculates an approximate time of vernal
!  equinox in days measured from 2000 January 1, hour 0.
!
!  Vernal assumes that vernal equinoxes from one year to the next
!  are separated by exactly 365.2425 days, a tropical year
!  [Explanatory Supplement to the Astronomical Ephemeris].  If the
!  tropical year is 365.2422 days, as indicated by other references,
!  then the time of the vernal equinox will be off by 2.88 hours in
!  400 years.
!
!  Time of vernal equinox for year 2000 A.D. is March 20, 7:36 GMT
!  [NASA Reference Publication 1349, Oct. 1994].  Vernal assumes
!  that vernal equinox for year 2000 will be on March 20, 7:30, or
!  79.3125 days from 2000 January 1, hour 0.  Vernal equinoxes for
!  other years returned by vernal are also measured in days from
!  2000 January 1, hour 0.  79.3125 = 31 + 29 + 19 + 7.5/24.

    real(8), parameter      :: ve2000=79.3125d0
    real(8), intent(in)     :: edayzy
    integer(4), intent(in)  :: iyear

    vernal = ve2000 + dble((iyear-2000))*edayzy

    end function vernal

end subroutine GISS_srevents

!*****************************************************************************************************************

!*****************************************************************************************************************
SUBROUTINE DtoYMDHM (DAY, IYEAR,IMONTH,IDATE,IHOUR,IMINUT)
!*****************************************************************************************************************
! subroutines based on
! SREVENTS.FOR    Solar EVENTS each year    2012/05/29
! https://data.giss.nasa.gov/ar5/SOLAR/SREVENTS.FOR downloaded 2017-09-12
!
!  DtoYMDHM receives DAY measured since 2000 January 1, hour 0 and
!  returns YEAR, MONTH, DATE, HOUR and MINUTE based on the Gregorian
!  calendar.

    implicit none
    
    external DtoYMD

    real(8), intent(in)     :: day
    integer(4), intent(out) :: IYEAR,IMONTH,IDATE,IHOUR,IMINUT

    real(8)                 :: date

    call DtoYMD (DAY, IYEAR,IMONTH,DATE)
    IDATE  = int(DATE+1.)
    IMINUT = Nint ((DATE-IDATE+1)*24*60)
    IHOUR  = IMINUT / 60
    IMINUT = IMINUT - IHOUR*60
    Return

end subroutine DtoYMDHM
!*****************************************************************************************************************
SUBROUTINE DtoYMD (DAY, IYEAR,IMONTH,DATE)
!*****************************************************************************************************************
! subroutines based on
! SREVENTS.FOR    Solar EVENTS each year    2012/05/29
! https://data.giss.nasa.gov/ar5/SOLAR/SREVENTS.FOR downloaded 2017-09-12
!
!  For a given DAY measured from 2000 January 1, hour 0, determine
!  the IYEAR (A.D.), IMONTH and DATE (between 0. and 31.).

    implicit none

    real(8), intent(in)     :: day
    integer(4), intent(out) :: iyear, imonth
    real(8), intent(out)    :: date

    real(8), PARAMETER  :: JDAY4C = 365*400 + 97, &     !  number of days in 4 centuries
                            JDAY1C = 365*100 + 24, &    !  number of days in 1 century
                            JDAY4Y = 365*  4 +  1, &    !  number of days in 4 years
                            JDAY1Y = 365               !  number of days in 1 year

    integer(4)  :: m
    INTEGER(4)  :: JDSUMN(12),JDSUML(12), n4year, n1year
    INTEGER(8)  :: N4CENT, n1cent
    real(8)     :: day4c, day1c, day4y, day1y

    DATA JDSUMN /0,31,59, 90,120,151, 181,212,243, 273,304,334/
    DATA JDSUML /0,31,60, 91,121,152, 182,213,244, 274,305,335/

    N4CENT = FLOOR(DAY/JDAY4C)
    DAY4C  = DAY - N4CENT*JDAY4C
    N1CENT = int8((DAY4C-1)/JDAY1C)
    If (N1CENT > 0)  GoTo 10
!  First of every fourth century: 16??, 20??, 24??, etc.
    DAY1C  = DAY4C
    N4YEAR = int(DAY1C/JDAY4Y)
    DAY4Y  = DAY1C - N4YEAR*JDAY4Y
    N1YEAR = int((DAY4Y-1)/JDAY1Y)
    If (N1YEAR > 0)  GoTo 200
    GoTo 100
!  Second to fourth of every fourth century: 21??, 22??, 23??, etc.
10 DAY1C  = DAY4C - N1CENT*JDAY1C - 1
    N4YEAR = int((DAY1C+1)/JDAY4Y)
    If (N4YEAR > 0)  GoTo 20
!  First four years of every second to fourth century when there is
!  no leap year: 2100-2103, 2200-2203, 2300-2303, etc.
    DAY4Y  = DAY1C
    N1YEAR = int(DAY4Y/JDAY1Y)
    DAY1Y  = DAY4Y - N1YEAR*JDAY1Y
    GoTo 210
!  Subsequent four years of every second to fourth century when
!  there is a leap year: 2104-2107, 2108-2111 ... 2204-2207, etc.
20 DAY4Y  = DAY1C - N4YEAR*JDAY4Y + 1
    N1YEAR = int((DAY4Y-1)/JDAY1Y)
    If (N1YEAR > 0)  GoTo 200
!
!  Current year is a leap frog year
!
100 DAY1Y = DAY4Y
    Do 120 M=1,11
120 If (DAY1Y < JDSUML(M+1))  GoTo 130
!     M=12
130 IYEAR  = 2000 + int(N4CENT*400) + int(N1CENT*100) + N4YEAR*4 + N1YEAR
    IMONTH = M
    DATE   = DAY1Y - JDSUML(M)
    Return
!
!  Current year is not a leap frog year
!
200 DAY1Y  = DAY4Y - N1YEAR*JDAY1Y - 1
210 Do 220 M=1,11
220 If (DAY1Y < JDSUMN(M+1))  GoTo 230
!     M=12
230 IYEAR  = 2000 + int(N4CENT)*400 + int(N1CENT)*100 + N4YEAR*4 + N1YEAR
    IMONTH = M
    DATE   = DAY1Y - JDSUMN(M)
    Return

end subroutine DtoYMD
!*****************************************************************************************************************
subroutine mon_to_day_ts(nt,imonlen,xm_in,xfill,no_negatives,smooth,restore,ndtot,nw,nsw,xd_out)
! Daily interpolation of a monthly time series
! Interpolation is done one year at a time, and so there can be small discontinuities between years.
! This version makes one pass over the input time series, and optionally smooths and restores the long-term mean.
!*****************************************************************************************************************

    implicit none
    
    external hdaily
    
    integer(4), parameter   :: nm=12, ndl=366
    integer(4), intent(in)  :: nt, ndtot, nw, nsw
    real(8), intent(in)     :: xm_in(nt), xfill
    integer(4), intent(in)  :: imonlen(nt)
    logical, intent(in)     :: no_negatives,smooth,restore
    real(8), intent(out)    :: xd_out(ndtot)
    
    real(8)                 :: xm(nm),xdh(ndl),xd_temp(ndtot)
    real(8)                 :: xm_ltm,xd_ltm,monlentot,ltmdiff
    real(8)                 :: pi,x,wgt(nw),wsum
    integer(4)              :: iml(nm),nd
    
    integer(4)              :: nyrs
    integer(4)              :: i,j,jj,jjj,js,m,mm,n,nn,nzero,nnonfill,nfill
    
    nyrs=nt/nm    
    
    pi=4.0d0*datan(1.0d0)
    
    ! generate smoothing weights
    do j=1,nw
        jj=j-(nw/2)-1
        x=(dble(jj)/((nw-1)/2.0d0))*4.0d0
        wgt(j)=(1.0d0/2.0d0*pi)*(exp(-0.5d0*x**2))
    end do
       
    ! intialize output variables
    xd_out=xfill; xd_temp=xfill
    
    ! main loop 
    n=0; mm=0
    do nn=1,nyrs         
        ! collect nm monthly values to interpolate this year, along with month lengths
        nd=0; nfill=0
        do m=1,nm
            mm=mm+1
            xm(m)=xm_in(mm)
            if (xm(m) .eq. xfill) nfill=nfill+1
            iml(m)=imonlen(mm)
            nd=nd+iml(m)
        end do
        
        ! check for fill value in any month, skip whole year if found
        if (nfill .eq. 0) then      
            
            ! mean-preserving daily interpolation
            call hdaily(nm,nd,xm,iml,no_negatives,xdh)           

            ! save daily values
            do i=1,nd
                n=n+1
                xd_out(n)=xdh(i)
            end do
        else        
            do i=1,nd
                n=n+1
                xd_out(n)=xfill
            end do      
        end if
    end do
    
    if (smooth) then
    
        ! save xd_out
        xd_temp = xd_out
        
        ! smooth across years
        n=imonlen(1)+imonlen(2)+imonlen(3)+imonlen(4)+imonlen(5)+imonlen(6) &
            +imonlen(7)+imonlen(8)+imonlen(9)+imonlen(10)+imonlen(11)+imonlen(12)
        mm=12
        do nn=1,nyrs-1
            jjj=n-(nsw/2)-1       
            do js=1,nsw
                jjj=jjj+1
                wsum=0.0d0; xd_out(jjj)=0.0d0
                jj=jjj-((nw-1)/2)-1
                do j=1,nw
                    jj=jj+1
                    if (xd_temp(jj) .ne. xfill) then
                        xd_out(jjj)=xd_out(jjj)+xd_temp(jj)*wgt(j)
                        wsum=wsum+wgt(j)
                    end if
                end do
                if (wsum.ne.0.0d0) then
                    xd_out(jjj)=xd_out(jjj)/wsum
                else
                    xd_out(jjj)=xfill
                end if

            end do

            n=n+imonlen(mm+1)+imonlen(mm+2)+imonlen(mm+3)+imonlen(mm+4)+imonlen(mm+5)+imonlen(mm+6) &
                +imonlen(mm+7)+imonlen(mm+8)+imonlen(mm+9)+imonlen(mm+10)+imonlen(mm+11)+imonlen(mm+12)
            mm=mm+12
        end do
        
    end if
    
    if (restore) then
    
        ! restore long-term mean
        xm_ltm=0.0d0; xd_ltm=0.0d0; monlentot=0.0d0
        if (no_negatives) then
            do mm=1,nt
                if (xm_in(mm).gt.0.0d0) then
                    if (xm_in(mm).ne.xfill) then
                        xm_ltm=xm_ltm+xm_in(mm)*dble(imonlen(mm))
                        monlentot=monlentot+dble(imonlen(mm))
                    end if
                end if
            end do
            nzero=0
            do n=1,ndtot
                if (xd_out(n).gt.0.0d0) then
                    if (xd_out(n).ne. xfill) then
                        xd_ltm=xd_ltm+xd_out(n)
                        nzero=nzero+1
                    end if
                end if
            end do

            if (monlentot.ne.0.0d0) then
                xm_ltm=xm_ltm/monlentot
            else 
                xm_ltm=xfill
            end if
            if (nzero.ne.0) then
                xd_ltm=xd_ltm/dble(nzero)
            else
                xd_ltm=xfill
            end if
            ltmdiff=xm_ltm-xd_ltm
        else
            do mm=1,nt      
                if (xm_in(mm).ne.xfill) then
                    xm_ltm=xm_ltm+xm_in(mm)*dble(imonlen(mm))
                    monlentot=monlentot+dble(imonlen(mm))
                end if
            end do
            nnonfill=0
            do n=1,ndtot
                if (xd_out(n).ne.xfill) then
                    xd_ltm=xd_ltm+xd_out(n)
                    nnonfill = nnonfill +1
                end if
            end do

            if (monlentot.ne.0.0d0) then
                xm_ltm=xm_ltm/monlentot
            else 
                xm_ltm=xfill
            end if
            if (nnonfill.ne.0) then
                xd_ltm=xd_ltm/dble(nnonfill)
            else
                xd_ltm=xfill
            end if
            ltmdiff=xm_ltm-xd_ltm
        end if
    
        if (no_negatives) then
            do n=1,ndtot
                if (xd_out(n).ne.xfill) then
                    if (xd_out(n).gt.0.0d0) then
                        if (xd_out(n)+ltmdiff.gt.0.0d0) then
                            xd_out(n)=xd_out(n)+ltmdiff
                        end if
                    end if
                end if
            end do
        else
            do n=1,ndtot
                if (xd_out(n).ne.xfill) then
                    xd_out(n)=xd_out(n)+ltmdiff
                end if
            end do
        end if
        
    end if
    
end subroutine mon_to_day_ts
!*****************************************************************************************************************
subroutine day_to_mon_ts(ny,ndays,rmonbeg,rmonend,ndtot,xd,xfill,xm_adj)
! aggregation of pseudo- or actual daily data to months using a paleo calendar
!*****************************************************************************************************************
    implicit none
    
    integer(4), parameter       :: nm=12, nd=366
    integer(4), intent(in)      :: ny, ndtot        ! number of years, total number of days
    integer(4), intent(in)      :: ndays(ny)        ! number of days in each year
    real(8), intent(in)         :: rmonbeg(ny,nm), rmonend(ny,nm)   ! beginning and ending days of each month
    real(8), intent(in)         :: xd(ndtot)        ! daily values
    real(8), intent(in)         :: xfill            ! _FillValue
    real(8), intent(out)        :: xm_adj(ny*nm)    ! (aggregated) average monthly values
    
    ! variables used to calculate monthly means
    integer(4)              :: ibegday, iendday             ! beginning day and ending day of each year
    integer(4)              :: ibeg(nm), iend(nm)           ! beginning and ending (integer) day of each month
    integer(4)              :: ndays_in_month(nm)           ! integer number of days in month
    real(8)                 :: xdx(-29:nd+30)               ! daily data for current year, padded by data from adjacent years
    real(8)                 :: wgt(-29:nd+30), wsum         ! weights (for interpolating over fractional days)
    integer(4)              :: nfill                        ! number of days with fill values
    
    integer(4)              :: n, m, i, nn
    
    ! loop over years, collecting daily data for each year, and getting monthly means
    iendday = 0; nn = 0
    xm_adj=0.0d0
    do n=1,ny
        ibegday = iendday + 1
        iendday = ibegday + ndays(n) - 1
        
        if (ny .eq. 1) then       ! single-year Aclim data  
            ! wrap the input daily data
            xdx(-29:0)=xd(ndays(n)-30+1:ndays(n))
            xdx(1:ndays(n))=xd(1:ndays(n))
            xdx(ndays(n)+1:ndays(n)+30)=xd(1:30)
        else 
            ! copy current year into xdx
            xdx(1:ndays(n)) = xd(ibegday:iendday)
            ! pad beginning and end of xdx
            if (n .eq. 1) then
                xdx(-29:0) = xd(ndays(n)-30+1:ndays(n))
                xdx(ndays(n)+1:ndays(n)+30) = xd(iendday+1:iendday+30)
            elseif (n .eq. ny) then
                xdx(-29:0) = xd(ibegday-30:ibegday-1)
                xdx(ndays(n)+1:ndays(n)+30) = xd(ibegday+1:ibegday+30)
            else
                xdx(-29:0) = xd(ibegday-30:ibegday-1)
                xdx(ndays(n)+1:ndays(n)+30) = xd(iendday+1:iendday+30)
            end if
        end if
    
        ! integer beginning and end of each month, and number of days in each month
        ! ndays_in_month should be equal to the integer month length + 1
        ibeg=ceiling(rmonbeg(n,:)); iend=ceiling(rmonend(n,:)); ndays_in_month=(iend-ibeg+1)
 
        ! monthly means
        do m=1,nm
            nn = nn + 1
            nfill = 0; wgt=1.0d0; wsum=0.0d0
            wgt(ibeg(m))=abs(rmonbeg(n,m)-dble(ibeg(m)))
            wgt(iend(m))=abs(rmonend(n,m)-dble(iend(m)-1))
            do i=ibeg(m),iend(m)
                if (xdx(i) .ne. xfill) then
                    xm_adj(nn)=xm_adj(nn)+xdx(i)*wgt(i)
                    wsum=wsum+wgt(i)
                else
                    nfill = nfill + 1
                end if
            end do
            if (wsum .ne. 0.0d0 .and. nfill .eq. 0) then
                xm_adj(nn)=xm_adj(nn)/wsum
            else
                xm_adj(nn)=xfill
            end if
        end do
    end do 

end subroutine day_to_mon_ts
!*****************************************************************************************************************
subroutine hdaily(nm,nd,xm,monlen,no_negatives,xdh)
!*****************************************************************************************************************
    implicit none
    
    external harmonic_coeffs, xdhat, dzero
    
    integer(4), parameter   :: nh=6
    integer(4), intent(in)  :: nd,nm
    real(8), intent(in)     :: xm(nm)
    integer(4), intent(in)  :: monlen(nm)
    logical, intent(in)     :: no_negatives
    real(8), intent(out)    :: xdh(nd)

    real(8)             :: a(0:nh),b(0:nh)
!    integer(4)          :: m
    
    ! interpolate daily values
    call harmonic_coeffs(nm,xm,a,b)
    call xdhat(nm,nd,monlen,a,b,xdh)
    if (no_negatives) call dzero(nm,nd,monlen,xm,xdh)

end subroutine hdaily
!*****************************************************************************************************************
subroutine harmonic_coeffs(nm,y,a,b)
! Calculates a's and b's of an "adjusted" harmonic fit to monthly values of a variable,
! which preserves the monthly (and annual) mean values by interpolated daily values
! adapted from Epstein (1991, On obtaining daily climatological values from monthly means,
! J. Climate 4:365-368).
!*****************************************************************************************************************
    implicit none
    
    integer(4), parameter   :: nh=6
    integer(4), intent(in)  :: nm
    real(8), intent(in)     :: y(nm)
    real(8), intent(out)    :: a(0:nh),b(0:nh)
    
    real(8)                 :: pi
    real(8)                 :: asum,bsum,c0,c1,c2,c3,c4
    integer(4)              :: j,t
    
    a=0.0d0; b=0.0d0
    pi=4.0d0*datan(1.0d0)
    
    ! a0
    do t=1,nm
        a(0)=a(0)+y(t)
    end do
    a(0)=a(0)/dble(nm) 
    
    ! a's and b's
    do j=1,nh-1
        asum=0.0d0; bsum=0.0d0
        c1=pi*(dble(j)/dble(nm))
        do t=1,nm
            c0=dble(t)/dble(nm)
            c2=(2.0d0*pi*dble(j)*dble(t))/dble(nm) 
            asum=asum+y(t)*dcos(c2)/(dble(nh))
            bsum=bsum+y(t)*dsin(c2)/(dble(nh))
        end do
        a(j)=(c1/dsin(c1))*asum
        b(j)=(c1/dsin(c1))*bsum
    end do
    
    asum=0.0d0
    do t=1,nm
        c3=cos(pi*dble(t))/dble(nm) 
        asum=asum+(y(t)*c3)
    end do
    c4=((pi/2.0d0)/sin(pi/2.0d0))
    a(nh)=c4*asum !((pi/2.0)/sin(pi/2.0))*asum
    b(nh)=0.0d0
    
end subroutine harmonic_coeffs
!*****************************************************************************************************************
subroutine xdhat(nm,nd,monlen,a,b,yhat)
! Calculates/interpolates pseudo-daily values of variable using the a's and b's from harmonic_coeffs()
! adapted from Epstein (1991, On obtaining daily climatological values from monthly means,
! J. Climate 4:365-368).
!*****************************************************************************************************************
    implicit none
    
    integer(4), parameter   :: nh=6
    integer(4), intent(in)  :: nm,nd
    integer(4), intent(in)  :: monlen(nm)
    real(8), intent(in)     :: a(0:nh),b(0:nh)
    real(8), intent(out)    :: yhat(nd)
    
    integer(4)              :: i,j,m,ii
    real(8)                 :: t,pi
    real(8)                 :: c2
    
    pi=4.0d0*datan(1.0d0)
    yhat=0.0
    ii=0
    do i=1,nm
        do m=1,monlen(i)
            ii=ii+1
            t=(dble(i)-0.5d0)+(dble(m)-0.5d0)/dble(monlen(i))
            do j=0,nh
                c2=((2.0d0*pi*dble(j)*t)/dble(nm))
                yhat(ii)=yhat(ii)+a(j)*dcos(c2)+b(j)*dsin(c2)
            end do
        end do
    end do
    
end subroutine xdhat
!*****************************************************************************************************************
subroutine dzero(nm,nd,monlen,xm,xd0)
! enforces 0.0 values of interpolated daily data when the monthly mean is 0.0
!*****************************************************************************************************************
    implicit none
    
    integer(4), intent(in)  :: nm,nd
    integer(4), intent(in)  :: monlen(nm)
    real(8), intent(in)     :: xm(nm)
    real(8), intent(inout)  :: xd0(nd)
    
    real(8)                 :: xdm(nm),diff,totaldiff
    integer(4)              :: i,m,j,nonzero(nm),l,maxiter=30
    integer(4)              :: imonlen(nm)
    
    imonlen = int(monlen)
    
    ! zero all negative daily values
    do i=1,nd
        if (xd0(i).le.0.0) xd0(i)=0.0
    end do
    
    do l=1,maxiter  
      ! zero daily values in months where xm=0.0
      i=0
      do m=1,nm
          do j=1,imonlen(m)
              i=i+1
              if (xm(m).eq.0.0) xd0(i)=0.0          
          end do
      end do
  
      i=0
      xdm=0.0d0; nonzero=0; totaldiff=0.0d0
      do m=1,nm
          do j=1,imonlen(m)
              i=i+1
              xdm(m)=xdm(m)+xd0(i)
              if (xdm(m).gt.0.0d0) nonzero(m)=nonzero(m)+1
          end do
          xdm(m)=xdm(m)/dble(monlen(m))
          diff=dabs((xm(m)-xdm(m)))
          totaldiff=totaldiff+diff
      end do
      if (totaldiff.le.0.0001) exit
    
      i=0
      do m=1,nm
          if (nonzero(m).ne.0) then
              do j=1,imonlen(m)
                  i=i+1
                  if (xd0(i).gt.0d0) then
                      xd0(i)=xd0(i)+(xm(m)-xdm(m)) !/nonzero(m)
                  end if
              end do
          end if
      end do    
    end do    
  
  ! zero daily values in months where xm=0.0
    i=0
    do m=1,nm
        do j=1,imonlen(m)
            i=i+1
            if (xm(m).eq.0.0) xd0(i)=0.0          
        end do
    end do
    
    ! zero all negative daily values
    do i=1,nd
        if (xd0(i).le.0.0) xd0(i)=0.0
    end do
    
end subroutine dzero
