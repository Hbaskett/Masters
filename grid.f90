! grid.f90 with kzz, dzz, vz, thermal and buoyancy diffusion
!> \file
!grid src with vz arrays added
!! Contains module mod_grid, subroutine mod_grid::init_grid() \n
!! Initialization of the pressure, temperature, optical depth, and frequency grid
!<

!> hold the parameters of the grid and the fields
!<
module mod_grid

  integer :: ndepth  !< number of depth for the atmosphere
  real    :: Teff    !< Effective temperature at the top of the atmosphere
  real    :: Tinit   !< Initial temperature for newton raphson iterations
  real    :: grav    !< Surface gravity
  real    :: logg    !< Surface gravity in log10(cm/s**2)
  real    :: Pmin    !< Pressure at the top of the atmosphere
  real    :: Pmax    !< Pressure at the bottom of the atmosphere
  real    :: taumin  !< Optical depth at the top of the atmosphere
  real    :: taumax  !< Optical depth at the bottom of the atmosphere
  integer :: nfreq   !< number of frequencies
  real    :: numin   !< Minimum frequency
  real    :: numax   !< Maximum frequency
  real    :: dnu     !< Frequency resolution
  real    :: nustd   !< Frequency used as reference for the optical depth

  integer :: ndim    !< Switch to the 2d version if =2
  real    :: Rmax    !< maximum extention of the atm if ndim=2
  integer :: ntheta  !< number of angles if ndim=2
  integer :: ndepth_2din
  real    :: u_cst               !< wind constant if ndim=2 and f2din=None
  real    :: alpha_dayside_cst   !< alpha constant if ndim=2 and f2din=None dayside
  real    :: alpha_nightside_cst !< alpha constant if ndim=2 and f2din=None night side
  real    :: unitU_2din          !< unit 2d wind
  real    :: unitP_2din          !< unit pressure profile wind
  real    :: pmin_2din           !< minimal pressure of 2d input profile
  real    :: pmax_2din           !< maximal pressure of 2d input profile
  real    :: theta_min !< minimum angle in 2d
  real    :: theta_max !< maximum angle in 2d
  real,dimension(:),allocatable    :: theta   !< value of the angles at cell centre
  real,dimension(:),allocatable    :: thetaf  !< value of the angles at cell interface
  real,dimension(:,:),allocatable  :: ppf_2d      !< pressure in 2d
  real,dimension(:,:),allocatable  :: ttf_2d      !< temperature in 2d
  real,dimension(:,:),allocatable  :: rhof_2d     !< density in 2d
  real,dimension(:,:),allocatable  :: cpf_2d      !< specific heat capacity in 2d
  real,dimension(:,:),allocatable  :: muf_2d      !< mean molecular weight in 2d
  real,dimension(:,:),allocatable  :: utf_2d      !< zonal wind velocity in 2d
  real,dimension(:,:),allocatable  :: dupf_2d     !< derivative of meridional velocity in 2d
  real,dimension(:,:),allocatable  :: urf_2d      !< vertical velocity in 2d

  real,dimension(:,:),allocatable  :: frad_2d     !< radiative flux in 2d
  real,dimension(:,:),allocatable  :: kapstd_2d   !< standard opacity in 2d
  real,dimension(:,:),allocatable  :: taufstd_2d  !< optical depth in 2d
  real,dimension(:,:,:),allocatable :: kapnu_2d   !< opacities in 2d
  real,dimension(:,:,:),allocatable :: epsfnu_2d  !< epsfnu in 2d
  real,dimension(:,:,:),allocatable :: dsnu_2d    !<  source term in 2d

  real,dimension(:),allocatable :: p_2din               !< pressure profile for 2d input
  real,dimension(:),allocatable :: u_2din               !< input wind profile
  real,dimension(:),allocatable :: alpha_dayside_2din   !< alpha input profile dayside
  real,dimension(:),allocatable :: alpha_nightside_2din !< alpha input profile nightside

  real    :: unitnugrid !< unit of the frequencies used in the namelist grid
  real    :: unitPgrid !< unit of the pressures used in the namelist grid
  real    :: unitKgrid !< unit of the temperatures used in the namelist grid

  logical        :: init_gradad !< init grad_ad if not in the init profile
  logical        :: corr_k   !< switch to use the k-correlated coefficients
  integer        :: nband    !< number of bands for the correlated-k
  integer        :: nkmix    !< number of coeff for the mixture for the correlated-k
  real :: unitnuband         !< units of the frequencies in nubandmin/max
  real   ,dimension(:),allocatable :: nubandmin !< minimum frequency in the band
  real   ,dimension(:),allocatable :: nubandmax !< maximum frequency in the band
  integer,parameter :: nbandmax = 6000 !< total number of bands

  integer :: nk            !< max number of coeff in a band
  integer :: nband_std      !< reference standard band for opacities
  real    :: Rp             !< Planetary radius at the bottom of the atmosphere
  real    :: pp_Rp          !< Pressure at which the radius is defined
  real    :: unitppRp       !< unit for Pressure at which the radius is defined

  real    :: period         !< orbital period
  real    :: omega          !< rotation rate of the planet


  logical :: check_tau      !< Update taustd on the fly
  real    :: error_tau
  integer :: mod_tau 		!< frequency to update taustd

  logical :: um_input    !< True if output from UM is used as input
  real    :: long_spec   !< Longitude at which to calculate 1D spectrum
  real    :: lat_spec    !< Latitude at which to calculate 1D spectrum
  real    :: long_obs    !< Longitude in degrees facing the observer. Used for 3D emission, transmission.
  real    :: lat_obs     !< Latitude in degrees facing the observer. Only used for 3D emission.
  real    :: scst_um     !< Stellar constant of UM

  real   ,dimension(:,:),allocatable :: weights !< weights for the correlated k
  integer,dimension(:,:),allocatable :: nkmax   !< maximum number of coeff

  real   ,dimension(:,:),allocatable :: weightsnu !< weights for the correlated k, allong the profil
  real   ,dimension(:,:),allocatable :: weightsfnu !< weights for the correlated k, allong the profil


  real,dimension(:  ),allocatable   :: rr     !< Radius at layer centre
  real,dimension(:  ),allocatable   :: rrf    !< Radius at layer interface
  real,dimension(:  ),allocatable   :: drf    !< thickness of layer, on cell face
  real,dimension(:  ),allocatable   :: dr     !< thickness of layer, on cell centre

  real,dimension(:  ),allocatable :: pp      !< Pressure profile of the atmosphere
  real,dimension(:  ),allocatable :: ppf     !< Pressure profile at layer interface
  real,dimension(:  ),allocatable :: tt      !< Temperature profile of the atmosphere
  real,dimension(:  ),allocatable :: ttf     !< Temperature profile at layer interface
  real,dimension(:  ),allocatable :: nn      !< Total number density at layer centre
  real,dimension(:  ),allocatable :: nnf     !< Total number density at intefaces
  real,dimension(:  ),allocatable :: gr      !< Gravity at layer centre
  real,dimension(:  ),allocatable :: grf     !< Gravity at interfaces
  real,dimension(:  ),allocatable :: phif     !< Gravity potential at interfaces

  real,dimension(:  ),allocatable :: radtime  !< radiative time profile at layer interface

  real,dimension(:  ),allocatable :: ppf_file     !< Pressure profile at layer interface form input file
  real,dimension(:  ),allocatable :: ttf_file     !< Temperature profile at layer interface from input file
  real,dimension(:  ),allocatable :: kzz_file     !< Eddy diffusion coefficient profile at layer interface from input file
  real,dimension(:  ),allocatable :: vz_file      !< Vertical advection coefficient profile at layer inferface from input file

  real :: time_sample_um    !< Time at which to samle UM output

  integer :: ntime  !< Number of output times in UM output
  integer :: nlat   !< Number of latitude points in UM output
  integer :: nlong  !< Number of longitude points in UM output
  integer :: npseudo
  integer :: nsw_band_um !< Number of SW bands in hires UM file
  integer :: nlw_band_um !< Number of LW bands in hires UM file

  real,dimension(:),allocatable :: time_um      !< Time in UM output
  real,dimension(:),allocatable :: height_um    !< Height in UM output
  real,dimension(:),allocatable :: heightf_um   !< Height at layer interface in UM output
  real,dimension(:),allocatable :: latitude_um  !< Latitude in UM output
  real,dimension(:),allocatable :: longitude_um !< Longitude in UM output
  real,dimension(:),allocatable :: pseudo_um
  real,dimension(:),allocatable :: sw_nubandmin_um
  real,dimension(:),allocatable :: sw_nubandmax_um
  real,dimension(:),allocatable :: lw_nubandmin_um
  real,dimension(:),allocatable :: lw_nubandmax_um
  real,dimension(:),allocatable :: irad_um

  real,dimension(:,:,:),allocatable :: pp_um  !< Pressure profiles from the UM
  real,dimension(:,:,:),allocatable :: tt_um  !< Temperature profiles from the UM
  real,dimension(:,:,:),allocatable :: ppf_um !< Pressure profiles at layer interface from the UM
  real,dimension(:,:,:),allocatable :: ttf_um !< Temperature profiles at layer interface from the UM
  real,dimension(:,:,:),allocatable :: sw_flux_um  !< Shortwave top of the atmosphere flux from the UM
  real,dimension(:,:,:),allocatable :: lw_flux_um  !< Longwave top of the atmosphere flux from the UM

  real,dimension(:  ),allocatable   :: rho     !< Density profile of the atmosphere
  real,dimension(:  ),allocatable   :: rhof    !< Density profile at layer interface
  real,dimension(:,:,:),allocatable :: rhof_um  !< Density profile at layer interface for the UM input
  real,dimension(:  ),allocatable   :: fconv   !< Convective flux profile at interfaces
  real,dimension(:  ),allocatable   :: kconv   !< Convective mixing profile at interfaces
  real,dimension(:,:),allocatable   :: dfconv  !< Derivatives of the convective flux profile at interfaces

  real,dimension(:  ),allocatable :: nuf_glob  !< Global array of frequency at bin interface
  real,dimension(:  ),allocatable :: nuf       !< local-cpu array of frequency at bin interface
  real,dimension(:  ),allocatable :: nu        !< local-cpu array of frequency at bin center
  real,dimension(:  ),allocatable :: taufstd  !< standard optical depth profile at interfaces between layers
  real,dimension(:  ),allocatable :: kapstd   !< standard opacity
  real,dimension(:  ),allocatable :: kapfstd  !< standard opacity at interfaces
  real,dimension(:  ),allocatable :: dkapstd  !< Derivatives of the opacity
  real,dimension(:  ),allocatable :: dkapfstd !< Derivatives of the opacity at interfaces

  real,dimension(:,:),allocatable :: epsnu !< Photon destruction probability =1 for no scattering
  real,dimension(:,:),allocatable :: kscatnu   !< Spectral scattering opacity
  real,dimension(:,:),allocatable :: epsfnu !< Photon destruction probability =1 for no scattering
  real,dimension(:,:),allocatable     :: kscatfnu  !< Spectral scattering opacity
  real,dimension(:,:),allocatable     :: kapnu     !< Spectral opacity
  real,dimension(:,:),allocatable     :: kapfnu    !< Spectral opacity at interfaces
  real,dimension(:,:,:,:),allocatable :: kapfnu_um !< Spectral opacity at interfaces for the UM
  real,dimension(:,:),allocatable     :: dkapnu    !< Derivatives of the spectral opacity
  real,dimension(:,:),allocatable     :: dkapfnu   !< Derivatives of the spectral opacity at interfaces
  real,dimension(:  ),allocatable :: frad      !< Radiative flux profile at interfaces
  real,dimension(:  ),allocatable :: fstar     !< irradiation flux profile at interfaces
  real,dimension(:  ),allocatable :: fuv       !< irradiation uv flux profile at interfaces
  real,dimension(:,:),allocatable :: fnu       !< Spectral flux profile at interfaces
  real,dimension(:,:),allocatable :: fstar_nu  !< Spectral flux profile at interfaces for the star
  real,dimension(:  ),allocatable :: qrad      !< Heating rate profile at interfaces
  real,dimension(:,:),allocatable :: qnu       !< Spectral heating rate profile at interfaces
  real,dimension(:  ),allocatable :: gradad    !< Local adiabatic gradient

  real,dimension(:,:),allocatable :: cf  !< Contribution Function as a function of wavelength and pressure
  real,dimension(:,:),allocatable :: ncf !< normalised contribution function
  real,dimension(:,:),allocatable :: taub !< lbl or correlated k optical depth

  real,dimension(:  ),allocatable :: kzz    !< eddy diffusion coefficient profile
  real                            :: kzzcst !< constant eddy diffusion coefficient
  logical                         :: mixing !< include vertical diffusion

  real,dimension(:  ),allocatable :: vz    !< vertical advection coefficient profile
  real                            :: vzcst !< constant vertical advecetion coefficient
  logical                         :: wind  !< include vertical advection

  logical :: isothermal !< to do an isothermal atmosphere

  real,dimension(:  ),allocatable :: err_hydro_prf !< Error in hydrostatic equilibrium along profile
  real,dimension(:  ),allocatable :: err_energy_prf !< Error in energy balance along profile
  
  real                            :: resolution    !< Resolution lambda/deltalambda
  real                            :: wvl_min       !< Minimum wavelength in microns
  real                            :: wvl_max       !< Maximum wavelength microns   
  logical :: uniform_R  !< logical for lambda/dlambda grid to use in ATMO

contains

  !> initialize the grid and the fields
  !<
  subroutine init_grid

    use mod_param
    use mod_cst
    use mod_util

    implicit none
    include 'netcdf.inc'

    real    :: dlogP,dlogtau,drr,dtheta,dlogwvl,dlogwvl_rounded
    integer,parameter :: nhead_sw = 48 !Number of header lines in sw spec file for UM
    integer,parameter :: nhead_lw = 52 !Number of header lines in lw spec file for UM
    integer :: stat,i,j,id_file,id_var,id_var2,id_var3,nlevel,ncol
    ! NetCDF variable IDs for UM input
    integer :: id_dim(4),id_dim2(4),id_dim3(4),idepth
    ! Name of netCDF dimension or variable
    character(len=32) :: name_dimvar
    character(len=3)  :: fifile
    ! Time index at which to sample UM output
    integer :: itime_um(1),it,iband,num_wvl
    real,dimension(:),allocatable   :: log_wvl, wvl    
    real,dimension(:),allocatable   ::  wvl_band_min,  wvl_band_max

    namelist /grid/ ndepth,Teff,Tinit,logg,Pmin,Pmax,Taumin,taumax,numin,numax,nfreq,nustd,unitnugrid,unitPgrid,unitKgrid,corr_k,&
    nband,nkmix,unitnuband,nband_std,Rp,um_input,long_spec,lat_spec,time_sample_um,long_obs,lat_obs,ndim,Rmax,ntheta,theta_min,&
    theta_max,period,check_tau,mod_tau,pp_Rp,unitppRp,u_cst,alpha_dayside_cst,alpha_nightside_cst,unitU_2din,unitP_2din,&
    ndepth_2din,pmin_2din,pmax_2din, isothermal, kzzcst, mixing, vzcst, wind, resolution, wvl_min , wvl_max, uniform_R, init_gradad


    ndepth 	 = 100
    Teff     = 1000.
    Tinit    = 1000.
    Pmin     = 1.0E-6
    Pmax     = 1.0E3
    Taumin   = 1.0E-5
    Taumax   = 100.
    numin    = 0.
    numax    = 1.0E6
    nustd    = 0.8E6
    nfreq    = 100
    logg     = 2.99
    Rp       = 0.
    period   = 1.
    pp_Rp    = 0.
    unitppRp = 1.0E6


    check_tau = .false.
    mod_tau = 5
    isothermal = .true.

    um_input          = .false.
    time_sample_um    = 0.
    long_spec   = 0.
    lat_spec    = 0.
    long_obs          = 0.
    lat_obs           = 0.

    corr_k     = .false.
    nband      = 32
    nkmix      = 15
    unitnuband = 0.01

    nband_std  = 19
    init_gradad = .true.

    ndim = 1
    ntheta = 5
    Rmax = 2.
    theta_min = 0.
    theta_max = pi

    ndepth_2din = 50
    u_cst = 1.
    alpha_dayside_cst   = 1E10
    alpha_nightside_cst = 1E10
    pmin_2din = 1.0E-9
    pmax_2din = 1.0E6

    unitU_2din = 1.   ! cm/s     in cgs
    unitP_2din = 1.   ! dyne/cm2 in cgs

    kzzcst = 0.
    mixing = .False.
    
    vzcst = 0.
    wind = .False.

    resolution = 1000.  !xshooter 10000
    wvl_min = 2.0E-1 !xshooter 0.2
    wvl_max = 3E+1 !xshooter 2.5
    uniform_R = .false.

    ! read grid namelist in the file fparam
    open(1,file=fparam)
    rewind(1) ; read(1,NML=grid,iostat=stat)
     IF (stat/=0) THEN
      BACKSPACE(1)
      READ(1,fmt='(A)') line
      WRITE(*,'(A)') 'Invalid line in namelist: '//trim(line)
      STOP
     END IF
    close(1)

    !convert units 2din
    u_cst =u_cst*unitU_2din/(unitL/unitT)
    pmin_2din = pmin_2din*unitP_2din/(unitP)
    pmax_2din = pmax_2din*unitP_2din/(unitP)

    !compute frequency resolution
    dnu = (numax-numin)/nfreq

    !convert to code units
    grav = 10**logg
    grav = grav/(unitL/unitT**2)
    Rp   = Rp*Rjup
    Rmax   = Rmax*Rjup

    period = period*day
    omega = 1./period

    !units of the namelist
    unitnugrid = 0.01   !m-1 -> cm-1
    unitPgrid  = 1.0E6  !bar -> dyne/cm2
    unitKgrid  = 1.

    pp_Rp = pp_Rp*unitppRp/(unitP)

    !convert from m-1 to cm-1 to code units
    numin = numin*unitnugrid/(1./unitL)
    numax = numax*unitnugrid/(1./unitL)
    nustd = nustd*unitnugrid/(1./unitL)
    dnu   = dnu  *unitnugrid/(1./unitL)

    !convert from bars to cgs to code units
    Pmin = Pmin*unitPgrid/unitP
    Pmax = Pmax*unitPgrid/unitP
    Teff = Teff*unitKgrid/unitK

    ! if correlated k used, initialize the bands
    if (corr_k) then
     if (.not. uniform_R) then
       allocate(nubandmin(nbandmax),nubandmax(nbandmax))
       nubandmin = 0.; nubandmax = 0.

       if (nband .le. 32) then

          nk = 31 !24 !21 ! changed to 31 for Fe

          nubandmax = (/   21700.,    50000.,    96200.,   155000.,   191600., 227300.,&
               263200.,   304100.,   334600.,   399200.,   460800.,   495000.,&
               562700.,   627700.,   668000.,   751900.,   835400.,   909100.,&
               995000.,  1041700.,  1098900.,  1162800.,  1273900.,  1342300.,&
               1481500.,  1634000.,  1748300.,  2020200.,  2500000.,  2800000.,&
               3176100.,  5000000./)

          nubandmin =  (/ 3.10000000e+03,   2.17000000e+04,   5.00000000e+04, &
               9.62000000e+04,   1.55000000e+05,   1.91600000e+05,&
               2.27300000e+05,   2.63200000e+05,   3.04100000e+05,&
               3.34600000e+05,   3.99200000e+05,   4.60800000e+05,&
               4.95000000e+05,   5.62700000e+05,   6.27700000e+05,&
               6.68000000e+05,   7.51900000e+05,   8.35400000e+05,&
               9.09100000e+05,   9.95000000e+05,   1.04170000e+06,&
               1.09890000e+06,   1.16280000e+06,   1.27390000e+06,&
               1.34230000e+06,   1.48150000e+06,   1.63400000e+06,&
               1.74830000e+06,   2.02020000e+06,   2.50000000e+06,&
               2.80000000e+06,   3.17610000e+06 /)




       else if (nband .le. 500) then

          nk = 21 !19

          nubandmin(1) = 0.!1.0E2
          nubandmax(1) = 1.0E4

          do iband = 2,500

             nubandmin(iband) = nubandmin(iband-1) + 1.0E4
             nubandmax(iband) = nubandmax(iband-1) + 1.0E4

          enddo

           nubandmin(1) = 1.0E2

       else if (nband .le. 5000) then

          nk = 17 

          nubandmin(1) = 0.!1.0E2
          nubandmax(1) = 1.0E3

          do iband = 2,5000

             nubandmin(iband) = nubandmin(iband-1) + 1.0E3
             nubandmax(iband) = nubandmax(iband-1) + 1.0E3

          enddo

          nubandmin(1) = 1.0E1

       endif
     endif
       
     if (uniform_R) then
       
       nk = 17 !10 for xshooter in order to optimize in memory
       
       !! Construct wvl grid same as the k-coeff file and then convert to wavenumber grid 

       dlogwvl = 1/resolution
       num_wvl = NINT((log(wvl_max) - log(wvl_min))/dlogwvl) !! not used ahead but should be equal to nband + 1

       dlogwvl_rounded = (log(wvl_max) - log(wvl_min))/nband
       !print *, "nband =", nband
       allocate(log_wvl(nband+1),wvl(nband+1)) !! num_wvl = nband + 1
       log_wvl(1) = log(wvl_min)
       
       do i = 2, nband+1
        log_wvl(i) = log_wvl(i-1) + dlogwvl_rounded
       enddo
       wvl = exp(log_wvl)
     
       !print *, "wvl =", wvl
       
       allocate(wvl_band_min(nband),wvl_band_max(nband))
       wvl_band_min = 0.; wvl_band_min = 0.
       allocate(nubandmin(nband),nubandmax(nband))
       nubandmin = 0.; nubandmax = 0.
       
       do iband = 1, nband !! 1 less because bands are one less than total number of wavelength points
        wvl_band_min(iband) = wvl(iband)
        wvl_band_max(iband) = wvl(iband+1)
        nubandmax(iband) = ((1.0*1E2)/(wvl_band_min(iband)*1E-4))  !max and min interchange for wavenumber extra 1E2 for ATMO units
        nubandmin(iband) = ((1.0*1E2)/(wvl_band_max(iband)*1E-4))  !max and min interchange in wavenumber extra 1E2 for ATMO units
       enddo
       
       nubandmax = nubandmax((nband):1:-1)
       !print *, "nubandmax = ", size(nubandmax)
       nubandmin = nubandmin((nband):1:-1)
       
        
     endif   

     nubandmin = nubandmin*unitnuband/(1./unitL)
     nubandmax = nubandmax*unitnuband/(1./unitL)
     
    endif

    ! input wind and alpha profile
    if (ndim==2) then

       ! constant profiles if no input profile
       if (f2din=='None') then

          allocate(p_2din(ndepth_2din+1),u_2din(ndepth_2din+1))
          allocate(alpha_dayside_2din(ndepth_2din+1))
          allocate(alpha_nightside_2din(ndepth_2din+1))

          dlogP = (log10(pmax_2din)-log10(pmin_2din))/(ndepth_2din)

          p_2din(1) = pmin_2din

          do idepth=2,ndepth_2din+1
             p_2din(idepth) = p_2din(idepth-1)*10.**dlogP
          enddo

          u_2din              (1:ndepth_2din+1) =               u_cst
          alpha_dayside_2din  (1:ndepth_2din+1) =   alpha_dayside_cst
          alpha_nightside_2din(1:ndepth_2din+1) = alpha_nightside_cst

       else

          call nf(nf_open(f2din,nf_nowrite,id_file))

          call nf(nf_inq_dimid (id_file,'nlevel_2din',id_var))
          call nf(nf_inq_dimlen(id_file,id_var,nlevel))

          ndepth_2din = nlevel-1

          allocate(p_2din(ndepth_2din+1),u_2din(ndepth_2din+1))
          allocate(alpha_dayside_2din(ndepth_2din+1))
          allocate(alpha_nightside_2din(ndepth_2din+1))

          ! Read pressures input profile
          call nf(nf_inq_varid(id_file,'p_2din' ,id_var))
          call nf(nf_get_vara_double(id_file,id_var,(/1/),(/nlevel/),p_2din(1:nlevel)))

          ! Read wind input profile
          call nf(nf_inq_varid(id_file,'u_2din' ,id_var))
          call nf(nf_get_vara_double(id_file,id_var,(/1/),(/nlevel/),u_2din(1:nlevel)))

          ! Read alpha dayside input profile
          call nf(nf_inq_varid(id_file,'alpha_dayside_2din' ,id_var))
          call nf(nf_get_vara_double(id_file,id_var,(/1/),(/nlevel/),alpha_dayside_2din(1:nlevel)))

          ! Read alpha nightside input profile
          call nf(nf_inq_varid(id_file,'alpha_nightside_2din' ,id_var))
          call nf(nf_get_vara_double(id_file,id_var,(/1/),(/nlevel/),alpha_nightside_2din(1:nlevel)))
          p_2din = p_2din*unitP_2din/unitP
          u_2din = u_2din*unitU_2din/(unitL/unitT)

       endif


    endif



    ! Isothermal atmosphere if no input profile is specified
    if (fin=='None') then

       ! allocate arrays
       allocate(tt(ndepth),pp(ndepth),rho(ndepth),nn(ndepth),gr(ndepth))
       allocate(ttf(ndepth+1),ppf(ndepth+1),rhof(ndepth+1),nnf(ndepth+1),grf(ndepth+1),phif(ndepth+1))
       allocate(rr(ndepth),rrf(ndepth+1),dr(ndepth),drf(ndepth+1))
       allocate(frad(ndepth+1),qrad(ndepth+1),fconv(ndepth+1),fuv(ndepth+1),fstar(ndepth+1),kconv(ndepth+1))
       allocate(taufstd(ndepth+1))
       allocate(kapstd(ndepth),kapfstd(ndepth+1))
       allocate(dkapstd(2*ndepth),dkapfstd(2*(ndepth+1)))
       allocate(nuf_glob(nfreq+1))
       allocate(dfconv(ndepth+1,2*(ndepth+1)))
       allocate(gradad(ndepth+1))
       allocate(kzz(ndepth+1))
       allocate(vz(ndepth+1))
       allocate(err_hydro_prf(ndepth+1))
       allocate(err_energy_prf(ndepth+1))

       if (corr_k) then
       allocate(taub(ndepth+1,nband))
       allocate(cf(ndepth,nband))
       allocate(ncf(ndepth,nband))
       else
       allocate(taub(ndepth+1,nfreq))
       allocate(cf(ndepth,nfreq))
       allocate(ncf(ndepth,nfreq))
       endif


       fstar = 0. ; fuv = 0.

       gradad = 0.

       nn = 0. ; nnf = 0.
       rho = 0. ; rhof = 0.
       gr  = grav ; grf = grav
       phif = 0.
       kzzcst = 0.
       vzcst = 0.

       ! build an isothermal atmosphere between Pmin and Pmax
       dlogP = (log10(Pmax)-log10(Pmin))/ndepth
       pp(1)  = Pmin * 10**(0.5*dlogP)

       ppf(1) = Pmin
       ppf(ndepth+1) = Pmax
       do i =2,ndepth
         pp (i) = pp (i-1) * 10**dlogP
         ppf(i) = ppf(i-1) * 10**dlogP
       end do

       ! Set (constant) Kzz
       kzz = kzzcst/(unitL**2/unitT)

       ! Set (constant) vz
       vz = vzcst/(unitL/unitT)

       !if isothermal
       if (isothermal) then
         tt(1)  = teff
         ttf(1) = teff
         ttf(ndepth+1) = teff

         do i = 2,ndepth
           tt(i) = teff
           ttf(i) = teff
         enddo
        !Use initial temperature guess
        else
        tt(1) = Tinit
        ttf(1) = Tinit
        ttf(ndepth+1) = Tinit
        do i = 2,ndepth
          tt(i) = Tinit
          ttf(i) = Tinit
        end do

        end if

    ! Input file is output from the UM
    else if (um_input) then

      ! Open netCDF file with UM output
      call nf(nf_open(fin,nf_nowrite,id_file))

      ! Get variable IDs
      call nf(nf_inq_varid(id_file,'theta',id_var))

      ! Get variable's dimension IDs
      call nf(nf_inq_vardimid(id_file,id_var,id_dim))

      ! Get length of dimensions
      ! Dimensions are ordered as follows:
      ! id_dims(4) = time
      ! id_dims(3) = height
      ! id_dims(2) = latitude
      ! id_dims(1) = longitude
      call nf(nf_inq_dimlen(id_file,id_dim(4),ntime))
      call nf(nf_inq_dimlen(id_file,id_dim(3),ndepth))
      call nf(nf_inq_dimlen(id_file,id_dim(2),nlat))
      call nf(nf_inq_dimlen(id_file,id_dim(1),nlong))
      nlevel = ndepth + 1

      ! Allocate grid arrays
      allocate(time_um(ntime))
      allocate(height_um(ndepth))
      allocate(heightf_um(nlevel))
      allocate(latitude_um(nlat))
      allocate(longitude_um(nlong))

      ! Allocate pressure and temperature arrays
      allocate(pp_um(nlong,nlat,ndepth))
      allocate(tt_um(nlong,nlat,ndepth))
      allocate(ppf_um(nlong,nlat,nlevel))
      allocate(ttf_um(nlong,nlat,nlevel))

      ! Read grid arrays. Make sure to get right correct grid for theta by
      ! using the name of the dimensions of theta.

      call nf(nf_inq_dimname(id_file,id_dim(4),name_dimvar))
      call nf(nf_inq_varid(id_file,trim(name_dimvar),id_var))
      call nf(nf_get_vara_double(id_file,id_var,1,ntime,time_um))

      call nf(nf_inq_dimname(id_file,id_dim(3),name_dimvar))
      call nf(nf_inq_varid(id_file,trim(name_dimvar),id_var))
      call nf(nf_get_vara_double(id_file,id_var,1,ndepth,height_um))

      call nf(nf_inq_dimname(id_file,id_dim(2),name_dimvar))
      call nf(nf_inq_varid(id_file,trim(name_dimvar),id_var))
      call nf(nf_get_vara_double(id_file,id_var,1,nlat,latitude_um))

      call nf(nf_inq_dimname(id_file,id_dim(1),name_dimvar))
      call nf(nf_inq_varid(id_file,trim(name_dimvar),id_var))
      call nf(nf_get_vara_double(id_file,id_var,1,nlong,longitude_um))

      ! Convert longitude and latitude to radians for convenience
      longitude_um = longitude_um*pi/180.
      latitude_um = latitude_um*pi/180.
      long_spec = long_spec*pi/180.
      lat_spec = lat_spec*pi/180.
      long_obs = long_obs*pi/180.
      lat_obs = lat_obs*pi/180.

      ! Convert latitude from being between -pi/2 to pi/2 to 0 and pi to be
      ! consistent with a spherical coordinate system where z is towards the
      ! north pole. At the moment latitude_um(1) is closest to the south pole,
      ! latitude_um(nlat) is closest to the north pole.
      latitude_um = pi/2. - latitude_um
      lat_spec = pi/2. - lat_spec
      lat_obs = pi/2. - lat_obs

      ! Find correct time index at which to sample UM output
      itime_um = minloc(abs(time_um - time_sample_um))
      if (mype == cputerm) then
        write(*,'(a,f8.3,a)') ' Samping UM output at ', time_um(itime_um), ' days'
      endif

      ! Read potential temperature and pressure
      call nf(nf_inq_varid(id_file,um_theta_name,id_var))
      call nf(nf_get_vara_double(id_file,id_var,(/1,1,1,itime_um/), &
          (/nlong,nlat,ndepth,1/),tt_um))
      call nf(nf_inq_varid(id_file,um_p_name,id_var))
      call nf(nf_get_vara_double(id_file,id_var,(/1,1,1,itime_um/), &
          (/nlong,nlat,ndepth,1/),pp_um))
      call nf(nf_close(id_file))

      ! At the moment latitude pi is the first element in latitude_um. Reverse
      ! the order so that latitude is strictly increasing.
      latitude_um = latitude_um(nlat:1:-1)
      tt_um = tt_um(:,nlat:1:-1,:)
      pp_um = pp_um(:,nlat:1:-1,:)

      ! Convert potential temperature to actual temperature
      tt_um = tt_um*(pp_um/pp_ref_surf_um)**(gas_constant_um/specific_heat_um)

      ! Convert temperature and pressure from SI to internal code units
      pp_um = pp_um*newton/metre**2
      tt_um = tt_um*kelvin
      height_um = height_um*metre

      ! Invert depth dimension so that the bottom of the atmosphere is located at ndepth
      pp_um = pp_um(:,:,ndepth:1:-1)
      tt_um = tt_um(:,:,ndepth:1:-1)
      height_um = height_um(ndepth:1:-1)

      ! Calculate height at layer interfaces
      heightf_um(1) = height_um(1) + (height_um(1) - height_um(2))/2.0
      heightf_um(2:nlevel-1) = (height_um(1:ndepth-1) + height_um(2:ndepth))/2.0
      heightf_um(nlevel) = height_um(ndepth) - (height_um(ndepth-1) - &
          height_um(ndepth))/2.0

      ! Interpolate pressure and temperature at levels
      ppf_um(:,:,1) = 10**((LOG10(pp_um(:,:,2)) - LOG10(pp_um(:,:,1)))/ &
          (height_um(2) - height_um(1))*(heightf_um(1) - height_um(1)) + &
          LOG10(pp_um(:,:,1)))
      ttf_um(:,:,1) = (tt_um(:,:,2) - tt_um(:,:,1))/(height_um(2) - &
          height_um(1))*(heightf_um(1) - height_um(1)) + tt_um(:,:,1)
      do i=1,nlong
        do j=1,nlat
          ppf_um(i,j,2:nlevel-1) = 10**( &
              (LOG10(pp_um(i,j,2:ndepth)) - LOG10(pp_um(i,j,1:ndepth-1)))/ &
              (height_um(2:ndepth) - height_um(1:ndepth-1))* &
              (heightf_um(2:nlevel-1) - height_um(1:ndepth-1)) + &
              LOG10(pp_um(i,j,1:ndepth-1)))

          ttf_um(i,j,2:nlevel-1) = (tt_um(i,j,2:ndepth) - &
              tt_um(i,j,1:ndepth-1))/(height_um(2:ndepth) - &
              height_um(1:ndepth-1))*(heightf_um(2:nlevel-1) - &
              height_um(1:ndepth-1)) + tt_um(i,j,1:ndepth-1)
        end do
      end do
      ppf_um(:,:,nlevel) = 10**( &
          (LOG10(pp_um(:,:,ndepth)) - LOG10(pp_um(:,:,ndepth-1)))/ &
          (height_um(ndepth) - height_um(ndepth-1))* &
          (heightf_um(nlevel) - height_um(ndepth-1)) + &
          LOG10(pp_um(:,:,ndepth-1)))
      ttf_um(:,:,nlevel) = (tt_um(:,:,ndepth) - tt_um(:,:,ndepth-1))/ &
          (height_um(ndepth) - height_um(ndepth-1))*(heightf_um(nlevel) - &
          height_um(ndepth-1)) + tt_um(:,:,ndepth-1)

      ! Allocate remaining 3D arrays
      allocate(rhof_um(nlong,nlat,ndepth+1))

      ! Allocate 1D arrays
      allocate(tt(ndepth),pp(ndepth),rho(ndepth),nn(ndepth),gr(ndepth))
      allocate(ttf(ndepth+1),ppf(ndepth+1),rhof(ndepth+1),nnf(ndepth+1),grf(ndepth+1),phif(ndepth+1))
      allocate(rr(ndepth),rrf(ndepth+1),dr(ndepth),drf(ndepth+1))
      allocate(frad(ndepth+1),qrad(ndepth+1),fconv(ndepth+1),kconv(ndepth+1),fstar(ndepth+1),fuv(ndepth+1))
      allocate(taufstd(ndepth+1))
      allocate(kapstd(ndepth),kapfstd(ndepth+1))
      allocate(dkapstd(2*ndepth),dkapfstd(2*(ndepth+1)))
      allocate(nuf_glob(nfreq+1))
      allocate(dfconv((ndepth+1),2*(ndepth+1)))

       if (corr_k) then
       allocate(taub(ndepth+1,nband))
       allocate(cf(ndepth,nband))
       allocate(ncf(ndepth,nband))
       else
       allocate(taub(ndepth+1,nfreq))
       allocate(cf(ndepth,nfreq))
       allocate(ncf(ndepth,nfreq))
       endif

      ! Set pressure and temperature temporarily so it is not undefined
      call copy_um_pt(1, 1)

    else

       call nf(nf_open(fin,nf_nowrite,id_file))

       call nf(nf_inq_dimid (id_file,'nlevel',id_var))
       call nf(nf_inq_dimlen(id_file,id_var,nlevel))

       ! to read old files
       !  call nf(nf_inq_dimid (id_file,'ndepth',id_var))
       !  call nf(nf_inq_dimlen(id_file,id_var,nlevel))

!       ndepth = nlevel -1

       allocate(ttf_file(nlevel),ppf_file(nlevel),kzz_file(nlevel),vz_file(nlevel))
       allocate(tt(ndepth),pp(ndepth),rho(ndepth),nn(ndepth),gr(ndepth))
       allocate(ttf(ndepth+1),ppf(ndepth+1),rhof(ndepth+1),nnf(ndepth+1),grf(ndepth+1),phif(ndepth+1))
       allocate(rr(ndepth),rrf(ndepth+1),dr(ndepth),drf(ndepth+1))
       allocate(frad(ndepth+1),qrad(ndepth+1),fconv(ndepth+1),kconv(ndepth+1),fstar(ndepth+1),fuv(ndepth+1))
       allocate(taufstd(ndepth+1))
       allocate(kapstd(ndepth),kapfstd(ndepth+1))
       allocate(dkapstd(2*ndepth),dkapfstd(2*(ndepth+1)))
       allocate(nuf_glob(nfreq+1))
       allocate(dfconv((ndepth+1),2*(ndepth+1)))
       allocate(gradad(ndepth+1))
       allocate(kzz(ndepth+1))
       allocate(vz(ndepth+1))
       allocate(err_hydro_prf(ndepth+1))
       allocate(err_energy_prf(ndepth+1))

       if (corr_k) then
       allocate(taub(ndepth+1,nband))
       allocate(cf(ndepth,nband))
       allocate(ncf(ndepth,nband))
       else
       allocate(taub(ndepth+1,nfreq))
       allocate(cf(ndepth,nfreq))
       allocate(ncf(ndepth,nfreq))
       endif

       fstar = 0. ; fuv = 0.

       nn = 0. ; nnf = 0.
       rho = 0. ; rhof = 0.
       gr = grav ; grf = grav
       phif = 0.
       kzz = 0.
       vz = 0.

       ! Read gas temperatures
       call nf(nf_inq_varid(id_file,'temperature' ,id_var))
       call nf(nf_get_vara_double(id_file,id_var,1,nlevel,ttf_file))

       ! Read gas pressures
       call nf(nf_inq_varid(id_file,'pressure',id_var))
       call nf(nf_get_vara_double(id_file,id_var,1,nlevel,ppf_file))

       if (.not. init_gradad) then

         call try_gradad(nf_inq_varid(id_file,'gradad',id_var))
         call nf(nf_get_vara_double(id_file,id_var,1,nlevel,gradad))

       endif

       if (mixing) then
         if (kzzcst==0.) then
           ! Read kzz profile from file
           call nf(nf_inq_varid(id_file,'kzz' ,id_var))
           call nf(nf_get_vara_double(id_file,id_var,1,nlevel,kzz_file))
         else
           ! Set to vertically constant value
           kzz_file = kzzcst
         end if
       end if

       if (wind) then
         if (vzcst==0.) then
           ! Read vz profile from file
           call nf(nf_inq_varid(id_file,'vz' ,id_var))
           call nf(nf_get_vara_double(id_file,id_var,1,nlevel,vz_file))
         else
           ! Set to vertically constant value
           vz_file= vzcst
         end if
       end if

       call nf(nf_close(id_file))

       ! Convert to code units
       ttf_file = ttf_file*1./unitK
       ppf_file = ppf_file*1./unitP
       kzz_file = kzz_file/(unitL**2/unitT)
       kzzcst = kzzcst/(unitL**2/unitT)
       vz_file = vz_file/(unitL/unitT)
       vzcst = vzcst/(unitL/unitT)


       ! Number of grid points match
       if (ndepth == nlevel - 1) then
          ppf = ppf_file
          ttf = ttf_file
          kzz = kzz_file
          vz = vz_file
       else
         ! Build new pressure array
         dlogP = (log10(ppf_file(nlevel))-log10(ppf_file(1)))/ndepth

         ppf(1) = ppf_file(1)
         ppf(ndepth+1) = ppf_file(nlevel)

         do i=2,ndepth
            ppf(i) = ppf(i-1) * 10**dlogP
         enddo
         do i=1,ndepth+1
            call interp_1d(log10(ppf_file), ttf_file, log10(ppf(i)), ttf(i))
            call interp_1d(log10(ppf_file), kzz_file, log10(ppf(i)), kzz(i))
            call interp_1d(log10(ppf_file), vz_file, log10(ppf(i)), vz(i))
         enddo

       endif

       do i =1,ndepth

          pp(i) = sqrt(ppf(i)*ppf(i+1))
          tt(i) = sqrt(ttf(i)*ttf(i+1))

       enddo

       deallocate(ttf_file,ppf_file,kzz_file,vz_file)

    endif

    !compute optical depth resolution
    dlogtau = (log10(taumax)-log10(taumin))/ndepth

    taufstd(1) = taumin
    taufstd(ndepth+1) = taumax

    do i =2,ndepth

       taufstd(i) = taufstd(i-1)*10.**dlogtau

    enddo

    if (ndim==2 .and. .not. um_input) then

       allocate(theta(0:ntheta+1))
       allocate(thetaf(0:ntheta+2))
       allocate(ppf_2d(ndepth+1,0:ntheta+2))
       allocate(ttf_2d(ndepth+1,0:ntheta+2))
       allocate(rhof_2d(ndepth+1,0:ntheta+2))
       allocate(utf_2d(ndepth+1,0:ntheta+2))
       allocate(urf_2d(ndepth+1,0:ntheta+2))
       allocate(dupf_2d(ndepth+1,0:ntheta+2))
       allocate(cpf_2d(ndepth+1,0:ntheta+2))
       allocate(muf_2d(ndepth+1,0:ntheta+2))
       allocate(frad_2d(ndepth+1,0:ntheta+2))
       allocate(taufstd_2d(ndepth+1,0:ntheta+2))
       allocate(kapstd_2d(ndepth,0:ntheta+2))

       thetaf(1) = 0.
       thetaf(ntheta+1) = 2.*pi
       dtheta = (thetaf(ntheta+1)-thetaf(1))/ntheta

       do it = 2,ntheta+1

          thetaf(it) = thetaf(it-1)+dtheta
          theta(it-1) = 0.5*(thetaf(it-1)+thetaf(it))

       enddo
       thetaf(0) = thetaf(1)-dtheta
       thetaf(ntheta+2) = thetaf(ntheta+1)+dtheta

       theta(0) = theta(1)-dtheta
       theta(ntheta+1) = theta(ntheta)+dtheta

       rrf(ndepth+1) = Rp
       rrf(1) = Rmax
       drr = (rrf(1)-rrf(ndepth+1))/ndepth

       grf(ndepth+1) = grav
       phif(ndepth+1) = -(grav*rrf(ndepth+1))

       do idepth = ndepth,1,-1

          rrf(idepth) = rrf(idepth+1)+drr
          rr(idepth) = sqrt(rrf(idepth)*rrf(idepth+1))

          grf(idepth) = grf(idepth+1)*rrf(idepth+1)**2/rrf(idepth)**2
          gr(idepth) = grf(idepth+1)*rrf(idepth+1)**2/rr(idepth)**2
          phif(idepth) = phif(idepth+1)*rrf(idepth+1)/rrf(idepth)

       enddo

       if (fin=='None') then

          do it = 0,ntheta+2
             ppf_2d(:,it) = ppf(:)
             ttf_2d(:,it) = Teff
          enddo

          utf_2d = 0. ; rhof_2d = 0. ; cpf_2d = 0.; frad_2d = 0. ; muf_2d = 0. ; dupf_2d = 0. ; urf_2d = 0.


       else

          call nf(nf_open(fin,nf_nowrite,id_file))

          call nf(nf_inq_dimid (id_file,'nlevel',id_var))
          call nf(nf_inq_dimlen(id_file,id_var,nlevel))

          call nf(nf_inq_dimid (id_file,'ncol',id_var))
          call nf(nf_inq_dimlen(id_file,id_var,ncol))

          if(.not. ndepth == nlevel-1) then
             if (mype==cputerm) then
                write(*,'(A)') 'number of layers incompatible with input file'
                write(*,intpmfmt)'ndepth',ndepth
                write(*,intpmfmt)'nlevel',nlevel
                write(*,*)
                write(*,separfmt)
                stop
             endif
          endif

          ! if(.not. ntheta == ncol-1) then
          !    if (mype==cputerm) then
          !       write(*,'(A)') 'number of columns incompatible with input file'
          !       write(*,intpmfmt)'ntheta',ntheta
          !       write(*,intpmfmt)'ncol',ncol
          !       write(*,*)
          !       write(*,separfmt)
          !       stop
          !    endif
          ! endif

          ncol = min(ncol,ntheta+1)

          ! Read gas temperatures
          call nf(nf_inq_varid(id_file,'temperature_2d' ,id_var))
          call nf(nf_get_vara_double(id_file,id_var,(/1,1/),(/nlevel,ncol/),ttf_2d(1:nlevel,1:ncol)))

          ! Read gas pressures
          call nf(nf_inq_varid(id_file,'pressure_2d',id_var))
          call nf(nf_get_vara_double(id_file,id_var,(/1,1/),(/nlevel,ncol/),ppf_2d(1:nlevel,1:ncol)))

          ! Read wind
          call nf(nf_inq_varid(id_file,'wind_2d',id_var))
          call nf(nf_get_vara_double(id_file,id_var,(/1,1/),(/nlevel,ncol/),utf_2d(1:nlevel,1:ncol)))

          ! Read wind
          call nf(nf_inq_varid(id_file,'dwind_2d',id_var))
          call nf(nf_get_vara_double(id_file,id_var,(/1,1/),(/nlevel,ncol/),dupf_2d(1:nlevel,1:ncol)))

          ! ! Read wind
          call nf(nf_inq_varid(id_file,'windz_2d',id_var))
          call nf(nf_get_vara_double(id_file,id_var,(/1,1/),(/nlevel,ncol/),urf_2d(1:nlevel,1:ncol)))

          if (ncol <ntheta+1) then

             do it =ncol+1,ntheta+1

                ttf_2d(:,it) = ttf_2d(:,ncol)
                ppf_2d(:,it) = ppf_2d(:,ncol)
                utf_2d(:,it) = utf_2d(:,ncol)
                urf_2d(:,it) = urf_2d(:,ncol)
                dupf_2d(:,it) = dupf_2d(:,ncol)


             enddo

          endif


          call nf(nf_close(id_file))

          ! Convert to code units
          ttf_2d = ttf_2d*1./unitK
          ppf_2d = ppf_2d*1./unitP
          utf_2d = utf_2d*(unitT/unitL)
          urf_2d = urf_2d*(unitT/unitL)
          dupf_2d = dupf_2d*(unitT/unitL)

          ppf_2d(ndepth+1,:) = pmax
       endif
    endif


    if (mype == cputerm) then
       write(*,headrfmt) 'Init Grid'
       write(*,intpmfmt) 'Ndepth',ndepth
              if (ndim==2) then
                 write(*,intpmfmt) 'Ntheta',ntheta
                 write(*,exppmfmt) 'Rmax',Rmax/Rjup,'radius at the top of the atmosphere'
              else
                 write(*,exppmfmt) 'taumin',taumin,'Standard optical depth at the top of the atmosphere'
                 write(*,exppmfmt) 'taumax',taumax,'Standard optical depth at the bottom of the atmosphere'
              endif
       write(*,exppmfmt) 'Rp',Rp/Rjup,'Planetary radius at the bottom of the atmosphere'


       if (fin=='None') then
          write(*,exppmfmt) 'Pmin',Pmin*unitP/1.0E6,'[bar] Pressure at the top of the isothermal atmosphere'
          write(*,exppmfmt) 'Pmax',Pmax*unitP/1.0E6,'[bar] Pressure at the bottom of the isothermal atmosphere'
          write(*,fltpmfmt) 'Teff',Teff*unitK,'[K] Temperature of the isothermal atmosphere'
       endif
       write(*,*)
       if (corr_k) then
          write(*,logpmfmt) 'corr_k',corr_k,'if T uses the correlated-k coefficients'
          write(*,intpmfmt) 'nband',nband,'total number of bands used'
          write(*,intpmfmt) 'nband_std',nband_std,'Reference band for the optical depth scale'
          write(*,intpmfmt) 'Nfreq',nfreq,'nb of freq. in each band used for the Planck function'
       else
          write(*,intpmfmt) 'Nfreq',nfreq
          write(*,exppmfmt) 'Numin',numin/unitL*100.,'[m-1] Minimum frequency'
          write(*,exppmfmt) 'Numax',numax/unitL*100.,'[m-1] Maximum frequency'
          write(*,exppmfmt) 'Nustd',nustd/unitL*100.,'[m-1] Reference frequency for the optical depth'
          write(*,exppmfmt) 'dnu',dnu/unitL*100.,'[m-1] Frequency resolution'
       endif
       write(*,*)
       write(*,fltpmfmt) 'Teff',Teff*unitK,'[K] Effective temperature'
       write(*,fltpmfmt) 'logg',logg,'Log surface gravity'
       write(*,fltpmfmt) 'grav',grav*unitL/unitT**2/100.,'[m s-2] Surface gravity'
       write(*,logpmfmt) 'init_gradad',init_gradad,'if F gradad is read from fin'
       write(*,*)
       write(*,separfmt)
    endif



    ! build the frequency band

    dnu = (numax-numin)/nfreq

    nuf_glob(1) = numin
    nuf_glob(nfreq+1) = numax

    do i =2,nfreq

       nuf_glob(i) = nuf_glob(i-1) + dnu

    enddo




  end subroutine init_grid


  ! Copy UM output to 1D arrays for use by chemistry module
  subroutine copy_um_pt(ilong, ilat)

    implicit none

    ! Longitude and latitude indices
    integer,intent(in) :: ilong, ilat

    pp(:) = pp_um(ilong,ilat,:)
    ppf(:) = ppf_um(ilong,ilat,:)
    tt(:) = tt_um(ilong,ilat,:)
    ttf(:) = ttf_um(ilong,ilat,:)

  end subroutine copy_um_pt

  ! Copy 1D profiles to 3D arrays for storage
  subroutine save_um_grid_vars(ilong, ilat)

    implicit none

    ! Longitude and latitude indices
    integer,intent(in) :: ilong, ilat

    ! Loop index
    integer :: idepth
    integer :: inu

    do idepth=1,ndepth+1
      kapfnu_um(ilong,ilat,idepth,:) = kapfnu(idepth,:)
      rhof_um(ilong,ilat,idepth) = rhof(idepth)
    end do

  end subroutine save_um_grid_vars


    !> check and print the error messages of ncdf called functions but continue
    !<
    subroutine try_gradad(errcode)

      use mod_param,only:mype
      implicit none
      include 'netcdf.inc'
      integer,intent(in) :: errcode

      if(errcode.eq.NF_NOERR)then
         init_gradad = .false.
      else
         init_gradad = .true.
      endif

    end subroutine try_gradad

    subroutine update_tau

 	implicit none

 	real :: dlogP_ideal, dlogP_actual,dlogtau
 	integer :: idepth

	dlogP_ideal = 0.
	dlogP_actual = 1.

	dlogP_ideal = 0.

	!Calculate ideal dlogP
	dlogP_ideal = (log10(ppf(ndepth+1))-log10(ppf(1)))/ndepth

	!Calculate average dlogP
	dlogP_actual = log10(ppf(2))-log10(ppf(1))

	error_tau = abs(dlogP_actual - dlogP_ideal)
	if (error_tau > 5E-2) then
		if (dlogP_actual > dlogP_ideal) then

			taumin = taumin/1.5
                        write(*,*) 'Decreasing Value of Taumin', taumin

		endif

		if (dlogP_actual < dlogP_ideal) then
			taumin = taumin*1.5
                        write(*,*) 'Increasing Value of Taumin', taumin
		endif

		!write(*,*) 'Updated Value of Taumin: ', taumin
	endif

	!compute optical depth resolution
	dlogtau = (log10(taumax)-log10(taumin))/ndepth

	taufstd(1) = taumin
	taufstd(ndepth+1) = taumax

	do idepth =2,ndepth

		taufstd(idepth) = taufstd(idepth-1)*10.**dlogtau

	enddo


 end subroutine update_tau


end module mod_grid
