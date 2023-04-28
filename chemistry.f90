! chemistry.f90 kzz, dzz, vz, thermal and buoyancy diffusion
!> \file
!! Contains modules mod_chem, subroutine mod_chem::init_chem() \n
!! mod_chem::get_eq_chemistry
!! initialize the physical constants and units
!<

!> Hold subroutines which solve for local chemical equilibrium,  \n
!! currently Gibbs minimisation and analytical formula of Burrows and Sharp 1999
!<
MODULE mod_chem

INTEGER,PARAMETER :: nele    = 23 !< number of elemental species
INTEGER,PARAMETER :: lname   = 10 !< maximum name length
INTEGER,PARAMETER :: nparam  = 21 !< number of NASA coefficient parameters
INTEGER,PARAMETER :: nmolmax = 500 !< max number of molecules allowed

INTEGER :: nmol !< number of molecules
INTEGER :: nmolg !< number of gas phase molecules

INTEGER,ALLOCATABLE :: mol_ele(:,:) !< map of elements in molecules

REAL  :: MdH       !< metallicity
REAL  :: P0        !< standard pressure
REAL  :: chem_conv !< control step size

REAL  :: Aele(nele)   !< element abundances
REAL  :: ele_xfactor(nele) !< multiplying factor for elements
REAL  :: elemass(nele)!< element masses
REAL  :: rel(nele)    !< element radii
REAL  :: elsig(nele)  !< element sigma for molecular diffusion (atomic diffusion volume)
REAL  :: Tcrit(nmolmax) !< species removed from gas phase at critical temperature
REAL  :: Pcrit(nmolmax) !< species removed from gas phase at critical pressure
REAL  :: Acst(nmolmax)  !< manually set abundance of species
REAL  :: X_massfrac     !< hydrogen mass fraction
REAL  :: Y_massfrac     !< helium mass fraction
REAL  :: Z_massfrac     !< metal mass fraction
REAL  :: tfreeze_eq

REAL,ALLOCATABLE :: dmol(:) !< molecule diameter
REAL,ALLOCATABLE :: molmass(:) !< molecule masses
REAL,ALLOCATABLE :: therm_diff(:) !< thermal diffusion parameter
REAL,ALLOCATABLE :: sigma_v(:) !< Diffusion volumes

REAL,ALLOCATABLE :: Amol(:,:)  !< molecule abundances on cell centres
REAL,ALLOCATABLE :: Amolf(:,:) !< molecule abundances on cell faces
REAL,ALLOCATABLE :: Aele_lev(:,:) !< element abundances on each cell face
REAL,ALLOCATABLE :: cp_chem(:) !< total heat capacity
REAL,ALLOCATABLE :: mmw(:)     !< mean molecule weight
REAL,ALLOCATABLE :: molParam(:,:) !< NASA polynomial coefficients

CHARACTER(lname)  :: chem !< method of chemistry; 'eq' <- Gibbs , 'neq' <- kinetics, 'cst' <- constant
CHARACTER(lname)  :: elename(nele) !< element names
CHARACTER(256)     :: fAin !< input abundance file
CHARACTER(256)     :: fAeqout !< output abundance file for equilibrium chemistry
CHARACTER(256)     :: fAneqout     !< name of the output neq abundance file
CHARACTER(256)     :: fAele !< input chemical elements abundances file

CHARACTER(60)     :: fcoeff !< nasa coefficients file
CHARACTER(60)     :: fcoeffnine !< nasa coefficients file (9-parameters)
CHARACTER(2),ALLOCATABLE      :: molphase(:) !< phase of each species (gas or condensed)

CHARACTER(lname),ALLOCATABLE :: molname(:) !< molecule names

LOGICAL :: print_chem
LOGICAL :: cond_NH3 !< flag to include condensation of NH3 (ana)
LOGICAL :: cond_H2O !< flag to include condensation of H2O (ana)
LOGICAL :: rainout  !< flag to include rainout chemistry
LOGICAL :: calc_cp_chem !< flag to calculate heat capacity from chemistry
LOGICAL :: cp_chem_reac !< include reaction component in heat capacity calculation
LOGICAL :: calc_hform !< flag to calculate enthalpy of formation of each species
LOGICAL :: print_thermodynamics ! flag to print thermodynamics
LOGICAL :: init_eq !< flag to initialise neq with eq abundances at run time

LOGICAL :: chem_smooth !< flag to vertically smooth the abundance profiles
INTEGER :: kerchem_smooth !< number of neighbors used for smoothing (1 means smoothing on 3 cells)
LOGICAL :: change_chem_smooth !< flag to turn on chem smooth after a given number of iterations
INTEGER :: nit_chem_smooth !< number of iterations after which to turn on chem smooth
INTEGER :: imol_smooth_start
INTEGER :: imol_smooth_end
INTEGER :: nb_chem_smooth

logical :: freeze_rainout !< logical used to freeze rainout of a given element
integer :: rainout_level !< rainout level used to freeze rainout
real    :: rainout_min_ele !< the minimum element abundance before an element is rained out
character(60) :: elename_freeze !< name of element for which rainout is frozen
integer :: nelename_freeze !< integer used for indexing the element for rainout freeze

REAL :: accuracy_chem !< tolerance of the newton solver for equilibrium chemical abundances
REAL :: nitmax_chem   !< maximum number of iterations for the newton solver for equilibrium chemical abundances

CONTAINS

  !----------------------
  !Initialise chemistry
  !----------------------
  SUBROUTINE init_chem

  USE mod_param
  USE mod_cst
  USE mod_grid

  IMPLICIT NONE

  INTEGER :: stat,i,iend
  REAL    :: unitTcrit,unitPcrit

  NAMELIST /chemistry/ chem, fAin, fAeqout, fAneqout, fcoeff, fcoeffnine, MdH, Tcrit, &
                       Pcrit, tfreeze_eq, Acst, print_chem, cond_NH3, cond_H2O, chem_conv, rainout, &
                       calc_cp_chem, cp_chem_reac, calc_hform, print_thermodynamics, init_eq, chem_smooth, kerchem_smooth, &
                       elename_freeze, change_chem_smooth, nit_chem_smooth, rainout_min_ele, ele_xfactor, accuracy_chem, nitmax_chem,imol_smooth_start,imol_smooth_end,nb_chem_smooth, fAele

  !Initialise namelist variables
  chem    = 'eq'

  fAin    = 'None'
  fAeqout = 'chem_eq.ncdf'
  fAneqout = 'chem_neq.ncdf'
  fcoeff  = '../../chem/coeff_NASA_sc.dat'
  fcoeffnine = '../../chem/coeff9_NASA_sc.dat'
  fAele = 'None'

  MdH     = 0.
  ! Solar element abundances by default
  ele_xfactor = 1.

  Tcrit = 0.
  Pcrit = 1.0E6
  Acst = 0.
  tfreeze_eq = 0.
  chem_conv = 2.0

  print_chem = .True.
  cond_NH3   = .False.
  cond_H2O   = .False.
  rainout    = .False.
  calc_cp_chem    = .False.
  cp_chem_reac = .False.
  calc_hform = .False.
  print_thermodynamics = .False.
  init_eq = .false.

  chem_smooth = .False.
  kerchem_smooth = 1
  change_chem_smooth = .False.
  nit_chem_smooth = 1000
  imol_smooth_start = 1
  imol_smooth_end   = 1E6
  nb_chem_smooth = 1

  freeze_rainout = .False.
  elename_freeze = 'None'
  rainout_min_ele = 1e-15

  accuracy_chem = 1E-6
  nitmax_chem   = 1000

  !Read namelist
  OPEN(1,file=fparam)
  REWIND(1) ; READ(1,NML=chemistry,iostat=stat)
   IF (stat/=0) THEN
    BACKSPACE(1)
    READ(1,fmt='(A)') line
    WRITE(*,'(A)') 'Invalid line in namelist: '//trim(line)
    STOP
   END IF
  CLOSE(1)

  ! override print statement if debug = true
  IF (debug>0) print_chem =.true.

  !Convert Tcrit and Pcrit to code units
  unitTcrit = 1.     !default is Kelvin
  unitPcrit = 1.0E6  !default is bar
  Tcrit = Tcrit*unitTcrit/unitK
  Pcrit = Pcrit*unitPcrit/unitP

  !Set element abundances
  CALL set_elements

  !Read molecule list
  CALL read_molecules

  if(imol_smooth_end>nmol) imol_smooth_end = nmol

  !Read NASA coefficients
  ALLOCATE(molParam(nmol,nparam))
  CALL read_NASA_coeffs(nmol,molname,molParam)

  !calculate mol_ele, molmass, dmol
  CALL get_molmass_dmol

  IF (calc_cp_chem) THEN
      ALLOCATE(cp_chem(ndepth+1))
      cp_chem = 0.
  END IF

  ALLOCATE(mmw(ndepth+1))
  mmw = 0.

  ! Allocate and initialise element on levels array
  ALLOCATE(Aele_lev(ndepth+1,nele))
  DO i = 1, ndepth+1
    Aele_lev(i,:) = Aele(:)
  END DO

  !Write information to screen
  IF (mype == cputerm) THEN
    WRITE(*,headrfmt) 'Chemistry'
    WRITE(*,chrpmfmt) 'chem',chem,'ana, eq, neq, or cst'
    WRITE(*,logpmfmt) 'print_chem',print_chem,'print the iterations of the chemistry'
    WRITE(*,logpmfmt) 'chem_smooth',chem_smooth,'if T abundances are vertically smoothed'
    if (chem_smooth) THEN
      WRITE(*,intpmfmt) 'kerchem_smooth',kerchem_smooth,'number of neighbours used for smoothing'
    endif
    WRITE(*,strpmfmt) 'fAin',   'Input  Abundances: '//fAin
    WRITE(*,strpmfmt) 'fAout',  'Output Abundances: '//fAeqout
    WRITE(*,strpmfmt) 'fAele',    'Input Elements  Abundances: '//fAele
    WRITE(*,strpmfmt) 'fcoeff', 'Thermodynamic data: '//fcoeff
    if (.not. elename_freeze == 'None') then
      write(*,chrpmfmt) 'elename_freeze',elename_freeze,'Element for which rainout is frozen'
    endif
    write(*,*)
    ! Write element abundances
    WRITE(*,intpmfmt) 'nele', nele, 'Number of elements'
    WRITE(*,fltpmfmt) 'MdH',MdH,'Metallicity [M/H] = 0 for the solar'
    WRITE(*,'(8x,A6,6x,A11,2x,A18,2x,A18)') 'Name','Abundance/H', 'Relative Abundance', 'Multiplying Factor'
    DO i =1,nele
      WRITE(*,'(5x,i3,a2,1x,A10,es10.3,6x,es10.3,6x,es10.3)') i, ': ',elename(i), Aele(i), Aele(i)/sum(Aele(:)), ele_xfactor(i)
    END DO
    WRITE(*,fltpmfmt) 'C/O',Aele(iele('C'))/Aele(iele('O')),'Carbon-oxygen ratio for these parameters'
    WRITE(*,fltpmfmt) 'X', X_massfrac, 'Hydrogen mass fraction for these parameters'
    WRITE(*,fltpmfmt) 'Y', Y_massfrac, 'Helium mass fraction for these parameters'
    WRITE(*,fltpmfmt) 'Z', Z_massfrac, 'Metal mass fraction for these parameters'
    WRITE(*,*)
    WRITE(*,intpmfmt) 'Nmol',nmol,'Number of molecules used for chemistry'
    iend = 0
    WRITE(*,'(5x)',advance='no')
    DO i =1,nmol
      iend = iend + 1
      IF (iend>5) THEN
        iend = 1
        WRITE(*,*)
        WRITE(*,'(5x)',advance='no')
      END IF
      WRITE(*,'(i3,a2,A10,1x)',advance='no') i,': ',molname(i)
    ENDDO
    WRITE(*,*)
    DO i =1,nmol
      IF (.not. Acst(i) ==0.) THEN
        WRITE(*,'(i3,a2,A10,A7,es10.3)') i,': ',molname(i),' Acst: ',Acst(i)
      END IF
      IF (Tcrit(i)>0.) THEN
        WRITE (*,'(5x,i3,a2,A10,A12,es10.3)') i,': ',molname(i),' Tcrit [K]: ',Tcrit(i)
      END IF
      IF (Pcrit(i)<1.0E6*unitPcrit/unitP) THEN
        WRITE (*,'(5x,i3,a2,A10,A14,es10.3)') i,': ',molname(i),' Pcrit [bar]: ',Pcrit(i)*unitP/1.0E6
      END IF
    END DO
    WRITE(*,*)
    WRITE(*,separfmt)
  END IF

  END SUBROUTINE init_chem

  !----------------------
  !Get chemical equilibrium abundances using either Gibbs
  !minimisation ('eq') or analytical formula ('ana')
  !----------------------
  SUBROUTINE get_eq_chemistry

  !USE mod_chem_solver
  USE mod_param
  USE mod_grid
  USE mod_cst

  IMPLICIT NONE

  INTEGER :: i,idepth

  !Call chemistry schemes
  !'eq'  - Gibbs minimisation
  !'ana' - analytical formula B&S1999
  !'man' - manually set from namelist
  !'relax' - chemical relaxation

  !Gibbs Minimisation
  IF (chem == 'eq' .or. ((chem == 'neq' .or. chem == 'relax') .and. init_eq)) THEN

    !Standard pressure
    CALL get_gibbs_abundances

    !If Acst(i), override molecule abundance
    DO i=1,nmol
      IF (Acst(i)>0.) THEN
        Amolf(:,i) = Acst(i)
      END IF
    END DO

  !Analytical, Burrows and Sharp 1999
  ELSE IF (chem == 'ana') THEN

    CALL get_analytical_abundances!(ndepth+1,nmol,nele,Aele,Amolf,Acst,ppf,ttf,atm,elename,molname,mype,cputerm)

    !If Acst(i), override molecule abundance
    DO i=1,nmol
      IF (Acst(i)>0.) THEN
        Amolf(:,i) = Acst(i)
      END IF
    END DO

    !Manually set all abundances
  ELSE IF (chem == 'man') THEN

    !Do not call chemistry schemes and manually set abundances
    !Set all abundances to zero initially
    Amolf(:,:) = 0.

    !Set abundances to values from namelist (Acst)
    DO i = 1,nmol
      IF (Acst(i)>0.) THEN
      Amolf(:,i) = Acst(i)
      END IF
    END DO

    IF (mype==cputerm) THEN
      WRITE(*,headrfmt) 'Manually set abundances.'
      WRITE(*,*) 'Note: abundances scaled to ensure sum(Amolf) = 1.'
      DO i = 1,nmol
        IF (Acst(i) > 0.) THEN
          WRITE(*,'(a,a,a,f8.3)') 'Abundance of ', molname(i),' = ', Amolf(1,i)
        END IF
      END DO
      WRITE(*,separfmt)
    END IF

  END IF !End call chemistry schemes

  !Check that sum(Amolf) = 1.0
  DO i = 1,ndepth+1
    Amolf(i,:) = Amolf(i,:)/sum(Amolf(i,:))
  END DO

  !If local temperature < Tcrit or local pressure > Pcrit set abundance
  !to zero
  DO idepth=1,ndepth+1
    DO i=1,nmol

      IF (ttf(idepth)<Tcrit(i) .or. ppf(idepth)>Pcrit(i)) THEN
        Amolf(idepth,i) = 0.
      END IF

    END DO
  END DO

  !interpolate to get Amol at cell center, using geometric mean
  DO idepth = 1,ndepth
    Amol(idepth,:) = sqrt(Amolf(idepth,:)*Amolf(idepth+1,:))
  END DO

  !vertically smooth the abundance profiles if asked for
  IF (chem_smooth) CALL smooth_chemistry

  !write chemistry to file
  CALL write_chemistry

  IF (mype == cputerm .and. print_chem) THEN
    WRITE(*,*)
    WRITE(*,*)
    WRITE(*,separfmt)
  END IF

  END SUBROUTINE get_eq_chemistry

  !----------------------
  !Subroutine to vertically smooth the chemical abundance profiles
  !----------------------
  SUBROUTINE smooth_chemistry

    USE mod_grid

    implicit none

    real,dimension(ndepth+1,nmol)::Amolf_smooth
    real,dimension(ndepth,nmol)::Amol_smooth

    integer :: imin,imax,ifmin,ifmax,idepth,iker,ismooth

	kerchem_smooth = max(kerchem_smooth ,1)

    imin = kerchem_smooth+1
    imax = ndepth-kerchem_smooth

    ifmin = kerchem_smooth+1
    ifmax = ndepth+1-kerchem_smooth

    do ismooth = 1,nb_chem_smooth

    do idepth = ifmin,ifmax

      Amolf_smooth(idepth,imol_smooth_start:imol_smooth_end) = Amolf(idepth,imol_smooth_start:imol_smooth_end)

        do iker = 1,kerchem_smooth

          Amolf_smooth(idepth,imol_smooth_start:imol_smooth_end) = Amolf_smooth(idepth,imol_smooth_start:imol_smooth_end) +  Amolf(idepth-iker,imol_smooth_start:imol_smooth_end)+ Amolf(idepth+iker,imol_smooth_start:imol_smooth_end)

        enddo

        Amolf_smooth(idepth,imol_smooth_start:imol_smooth_end) = Amolf_smooth(idepth,imol_smooth_start:imol_smooth_end)/(2*kerchem_smooth+1)

    enddo

    do idepth = ifmin,ifmax

      Amolf(idepth,imol_smooth_start:imol_smooth_end) = Amolf_smooth(idepth,imol_smooth_start:imol_smooth_end)

    enddo

    do idepth = imin,imax

      Amol_smooth(idepth,imol_smooth_start:imol_smooth_end) = Amol(idepth,imol_smooth_start:imol_smooth_end)

      do iker = 1,kerchem_smooth

        Amol_smooth(idepth,imol_smooth_start:imol_smooth_end) = Amol_smooth(idepth,imol_smooth_start:imol_smooth_end) +  Amol(idepth-iker,imol_smooth_start:imol_smooth_end)+ Amol(idepth+iker,imol_smooth_start:imol_smooth_end)

      enddo

      Amol_smooth(idepth,imol_smooth_start:imol_smooth_end) = Amol_smooth(idepth,imol_smooth_start:imol_smooth_end)/(2*kerchem_smooth+1)

    enddo

    do idepth = imin,imax

      Amol(idepth,imol_smooth_start:imol_smooth_end) = Amol_smooth(idepth,imol_smooth_start:imol_smooth_end)

   enddo

   enddo

  END SUBROUTINE

  !----------------------
  !Subroutine to write chemistry abundances to file
  !----------------------
  SUBROUTINE write_chemistry

  USE mod_grid
  USE mod_util
  USE mod_param
  USE mod_cst

  IMPLICIT NONE
  INCLUDE 'netcdf.inc'

  INTEGER :: id_file,id_var1,id_var2,id_var3,id_var4,id_var5,id_var6,id_var7,id_lev,id_mol,id_name,id_ele
  CHARACTER(256) :: Aout

  IF ((chem=='eq') .OR. (chem=='ana')) THEN
    Aout = fAeqout
  ELSE
    Aout = fAneqout
  END IF

  IF (mype==cputerm) THEN

    CALL nf(nf_create(Aout,nf_clobber,id_file))

    CALL nf(nf_def_dim(id_file,'nlevel',ndepth+1,id_lev))
    CALL nf(nf_def_dim(id_file,'nmol',nmol,id_mol))
    CALL nf(nf_def_dim(id_file,'lname',lname,id_name))
    CALL nf(nf_def_dim(id_file,'nele',nele,id_ele))

    CALL nf(nf_def_var(id_file,'abundances',nf_double,2,(/id_lev,id_mol/),id_var1))
    CALL nf(nf_def_var(id_file,'pressure',nf_double,1,id_lev,id_var2))
    CALL nf(nf_def_var(id_file,'temperature',nf_double,1,id_lev,id_var7))
    IF (calc_cp_chem) THEN
        CALL nf(nf_def_var(id_file,'cp',nf_double,1,id_lev,id_var3))
    END IF
    CALL nf(nf_def_var(id_file,'mean_mol_mass',nf_double,1,id_lev,id_var4))
    CALL nf(nf_def_var(id_file,'molname',nf_char,2,(/id_name,id_mol/),id_var5))
    CALL nf(nf_def_var(id_file,'element_abundances',nf_double,2,(/id_lev,id_ele/),id_var6))
    CALL nf(nf_enddef(id_file))

    CALL nf(nf_put_vara_double(id_file,id_var1,(/1,1/),(/ndepth+1,nmol/),Amolf(:,1:nmol)))
    CALL nf(nf_put_vara_double(id_file,id_var2,1,ndepth+1,ppf*unitP))
    CALL nf(nf_put_vara_double(id_file,id_var7,1,ndepth+1,ttf*unitK))

    IF (calc_cp_chem) THEN
        CALL nf(nf_put_vara_double(id_file,id_var3,1,ndepth+1,cp_chem*unitE/unitM))
    END IF
    CALL nf(nf_put_vara_double(id_file,id_var4,1,ndepth+1,mmw))
    CALL nf(nf_put_vara_text(id_file,id_var5,(/1,1/),(/lname,nmol/),molname(1:nmol)))
    CALL nf(nf_put_vara_double(id_file,id_var6,(/1,1/),(/ndepth+1,nele/),Aele_lev(:,1:nele)))

    CALL nf(nf_close(id_file))
  END IF

  END SUBROUTINE write_chemistry

  !----------------------
  !Function to get the index of a molecule
  !----------------------
  FUNCTION imol(molrequest)

  USE mod_param

  IMPLICIT NONE

  INTEGER :: i,imol
  CHARACTER(*) :: molrequest

  imol = 0

  DO i = 1,nmol
    IF (molrequest == molname(i)) THEN
      imol = i
    END IF
  END DO

  IF (molrequest == 'HV        ') THEN !photochem reaction
    imol = -1
  ELSE IF (molrequest == '          ') THEN
    imol = -1
  END IF

  IF (imol==0) THEN
    IF (mype==cputerm) THEN
      WRITE(*,headrfmt) 'Error imol'
      WRITE(*,'(A)') 'Requested molecule not found'
      WRITE(*,chrpmfmt) 'Requested mol.',molrequest
      WRITE(*,*)
      WRITE(*,separfmt)
      STOP
    END IF
  END IF

  END FUNCTION imol

  !----------------------
  !Function to test if molecule is included
  !----------------------
  FUNCTION ismol(molrequest)

  USE mod_param

  IMPLICIT NONE

  INTEGER :: i
  LOGICAL :: ismol
  CHARACTER(*) :: molrequest

  ismol = .False.

  DO i = 1,nmol
    IF (molrequest == molname(i)) THEN
      ismol = .True.
    END IF
  END DO

  END FUNCTION ismol

  !----------------------
  !Function to get the index of an element
  !----------------------
  FUNCTION iele(elerequest)

  USE mod_param

  IMPLICIT NONE

  INTEGER :: i,iele
  CHARACTER(*) :: elerequest

  iele = 0

  DO i = 1,nele
    IF (elerequest == elename(i)) THEN
      iele = i
    END IF
  END DO

  IF (iele==0) THEN
    IF (mype==cputerm) THEN
      WRITE(*,headrfmt) 'Error iele'
      WRITE(*,'(A)') 'Requested element not found'
      WRITE(*,chrpmfmt) 'Requested ele.',elerequest
      WRITE(*,*)
      WRITE(*,separfmt)
      STOP
    END IF
  END IF

  END FUNCTION iele

END MODULE mod_chem


!----------------------
!Set element abundances
!----------------------
SUBROUTINE set_elements

USE mod_chem

IMPLICIT NONE

INTEGER :: iel

!Local Variables
REAL      :: abundance
CHARACTER(lname) :: elt_name


! Element names; abundances relative to hydrogen; masses; radii
! Abundances from Asplund+09 (photosphere) - Ar, Al, Ca, Rb
! Abundances from Asplund+09 (meteorites) - Li, Cs
! Abundances from Lodders+09 - C, He, N, O, P, S, Si, Ti, V, Cl, K, Na, Mg, F, Ca, Fe
! Sigma values from Poling et al. 2000, "The Properties of Gases and Liquids", McGraw Hill, Table 11.1
elename(1)  = 'H'   ; Aele(1)  = 1.0000E+00 ; elemass(1)  = 1.007947  ; rel(1)  = 31.   ; elsig(1)  = 2.31
elename(2)  = 'He'  ; Aele(2)  = 9.6917E-02 ; elemass(2)  = 4.0026022 ; rel(2)  = 28.   ; elsig(2)  = 2.67
elename(3)  = 'C'   ; Aele(3)  = 2.7759E-04 ; elemass(3)  = 12.01078  ; rel(3)  = 69.   ; elsig(3)  = 15.9
elename(4)  = 'N'   ; Aele(4)  = 8.1846E-05 ; elemass(4)  = 14.00672  ; rel(4)  = 71.   ; elsig(4)  = 4.54
elename(5)  = 'O'   ; Aele(5)  = 6.0618E-04 ; elemass(5)  = 15.99943  ; rel(5)  = 66.   ; elsig(5)  = 6.11
elename(6)  = 'Na'  ; Aele(6)  = 2.2279E-06 ; elemass(6)  = 22.989769 ; rel(6)  = 166.  ; elsig(6)  = 0.
elename(7)  = 'K'   ; Aele(7)  = 1.4518E-07 ; elemass(7)  = 39.0983   ; rel(7)  = 203.  ; elsig(7)  = 0.
elename(8)  = 'Si'  ; Aele(8)  = 3.8610E-05 ; elemass(8)  = 28.08553  ; rel(8)  = 111.  ; elsig(8)  = 0.
elename(9)  = 'Ar'  ; Aele(9)  = 2.5119E-06 ; elemass(9)  = 39.9481   ; rel(9)  = 106.  ; elsig(9)  = 16.2
elename(10) = 'Ti'  ; Aele(10) = 9.5367E-08 ; elemass(10) = 47.8671   ; rel(10) = 160.  ; elsig(10) = 0.
elename(11) = 'V'   ; Aele(11) = 1.1059E-08 ; elemass(11) = 50.94151  ; rel(11) = 153.  ; elsig(11) = 0.
elename(12) = 'S'   ; Aele(12) = 1.3183E-05 ; elemass(12) = 32.065    ; rel(12) = 105.  ; elsig(12) = 22.9
elename(13) = 'Cl'  ; Aele(13) = 1.9962E-07 ; elemass(13) = 35.453    ; rel(13) = 102.  ; elsig(13) = 21.0
elename(14) = 'Mg'  ; Aele(14) = 3.9765E-05 ; elemass(14) = 24.3050   ; rel(14) = 141.  ; elsig(14) = 0.
elename(15) = 'Al'  ; Aele(15) = 2.8184E-06 ; elemass(15) = 26.9815   ; rel(15) = 121.  ; elsig(15) = 0.
elename(16) = 'Ca'  ; Aele(16) = 2.3318E-06 ; elemass(16) = 40.078    ; rel(16) = 176.  ; elsig(16) = 0.
elename(17) = 'Fe'  ; Aele(17) = 3.2742E-05 ; elemass(17) = 56.078    ; rel(17) = 152.  ; elsig(17) = 0.
elename(18) = 'Cr'  ; Aele(18) = 4.3652E-07 ; elemass(18) = 51.9961   ; rel(18) = 139.  ; elsig(18) = 0.
elename(19) = 'Li'  ; Aele(19) = 1.8197E-09 ; elemass(19) = 6.941     ; rel(19) = 128.  ; elsig(19) = 0.
elename(20) = 'Cs'  ; Aele(20) = 1.2023E-11 ; elemass(20) = 132.9055  ; rel(20) = 244.  ; elsig(20) = 0.
elename(21) = 'Rb'  ; Aele(21) = 3.3113E-10 ; elemass(21) = 85.4378   ; rel(21) = 220.  ; elsig(21) = 0.
elename(22) = 'F'   ; Aele(22) = 3.1043E-08 ; elemass(22) = 18.9984   ; rel(22) = 64.   ; elsig(22) = 14.7
elename(23) = 'P'   ; Aele(23) = 3.2048E-07 ; elemass(23) = 30.9738   ; rel(23) = 100.  ; elsig(23) = 0.

! Multiply all elements heavier than H and He by metallicity factor
DO iel = 1, nele
  IF (elename(iel) /= 'H' .and. elename(iel) /= 'He') THEN
    Aele(iel) = Aele(iel)*10.**MdH
  END IF
END DO

! Multiply individual element abundance by individual element factor
! ele_xfactor = 1. by default, set in namelist
DO iel = 1, nele
  Aele(iel) = Aele(iel)*ele_xfactor(iel)

! Calculate hydrogen, helium, and metal mass fractions - X, Y and Z respectively
X_massfrac = Aele(1)*elemass(1)/sum(Aele(:)*elemass(:))
Y_massfrac = Aele(2)*elemass(2)/sum(Aele(:)*elemass(:))
Z_massfrac = sum(Aele(3:)*elemass(3:))/sum(Aele(:)*elemass(:))

END DO

! Calculate a value for the atomic diffusion for the elements where no data is available
! Uses the approximate linear relationship found between atomic radius and diffusion volumes
DO iel = 1, nele
  IF (elsig(iel) == 0.) THEN
       elsig(iel) = 0.3576*rel(iel) - 10.794
  END IF
END DO

END SUBROUTINE set_elements

!----------------------
!Read molecule list
!----------------------
SUBROUTINE read_molecules

USE mod_grid
USE mod_param
USE mod_util
USE mod_chem

IMPLICIT NONE

INTEGER :: i,nlevel,id_var,id_file

INCLUDE 'netcdf.inc'

!If molecule list not supplied, force user to provide one
IF (fAin == 'None' .and. ((chem == 'eq') .or. (chem == 'neq') .or. (chem == 'cst'))) THEN

  IF (mype == cputerm) THEN
    WRITE(*,headrfmt) 'Error init_chemistry'
    WRITE(*,'(A)') 'No molecule list provided'
    WRITE(*,chrpmfmt) 'fAin',fAin
    WRITE(*,'(A)') 'Provide netcdf file with list of molecules to include'
    WRITE(*,*)
    WRITE(*,separfmt)
    STOP
  END IF

  !Set up analytical chemistry
  ELSE IF (chem == 'ana') THEN

    nmol = 12

    ALLOCATE(Amol(ndepth,nmol))
    ALLOCATE(Amolf(ndepth+1,nmol))
    ALLOCATE(molname(nmol))

    molname(1) = 'H2'
    molname(2) = 'He'
    molname(3) = 'H2O'
    molname(4) = 'CO'
    molname(5) = 'CH4'
    molname(6) = 'NH3'
    molname(7) = 'N2'
    molname(8)  = 'Na'
    molname(9)  = 'K'
    molname(10) = 'Si'
    molname(11) = 'TiO'
    molname(12) = 'VO'

    Amolf = 0.0
    Amol  = 0.0

    !Read in chemistry input file
    !For chem = 'eq' this is list of molecules
    !For chem = 'neq' this is list of molecules + eq abundances
  ELSE IF ((chem == 'eq') .or. (chem == 'neq') .or. (chem == 'cst') .or. (chem == 'man') .or. (chem == 'relax')) THEN

    !Read netcdf file
    CALL nf(nf_open(fAin,nf_nowrite,id_file))
    !Get number of molecules
    CALL nf(nf_inq_dimid (id_file,'nmol',id_var))
    CALL nf(nf_inq_dimlen(id_file,id_var,nmol))

    !Allocate arrays
    ALLOCATE(Amol(ndepth,nmol))
    ALLOCATE(Amolf(ndepth+1,nmol))
    ALLOCATE(molname(nmol))

    !Get number of vertical levels
    CALL nf(nf_inq_dimid(id_file,'nlevel',id_var))
    CALL nf(nf_inq_dimlen(id_file,id_var,nlevel))
    !Get molecule names
    CALL nf(nf_inq_varid(id_file,'molname' ,id_var))
    CALL nf(nf_get_vara_text(id_file,id_var,(/1,1/),(/lname,nmol/),molname(1:nmol)))

    !If non-equilibrium or constant abundances, read abundances
    IF ((chem == 'neq' .or. chem == 'cst') .and. .not. init_eq) THEN

      IF (nlevel == ndepth+1) THEN
        CALL nf(nf_inq_varid(id_file,'abundances' ,id_var))
        CALL nf(nf_get_vara_double(id_file,id_var,(/1,1/),(/ndepth+1,nmol/),Amolf(:,1:nmol)))

        !Geometric mean at cell centre
        DO i =1,ndepth
          Amol(i,:) = sqrt(Amolf(i,:)*Amolf(i+1,:))
        END DO

      ELSE IF (ndepth == 1) THEN

        ! Box model
        CALL nf(nf_inq_varid(id_file,'abundances' ,id_var))
        CALL nf(nf_get_vara_double(id_file,id_var,(/1,1/),(/2,nmol/),Amolf(:,1:nmol)))

        Amol(1,:) = Amolf(1,:)

      ELSE

        IF (mype==cputerm) THEN
          WRITE(*,headrfmt) 'Error init_neq_chem'
          WRITE(*,*) nlevel
          WRITE(*,*) ndepth+1
          write(*,*) fAin
          WRITE(*,'(A)') 'The input abundances have not the same number of layers as the input pt profile'
          WRITE(*,chrpmfmt) 'chem',chem
          WRITE(*,'(A)') 'Run equilibrium chemistry before neq to initialize the abundances'
          WRITE(*,*)
          WRITE(*,separfmt)
          STOP
        END IF

      END IF

    ELSE
      Amol  = 0.0
      Amolf = 0.0
    END IF

    !Close file
    CALL nf(nf_close(id_file))
  END IF

END SUBROUTINE read_molecules

!----------------------
!Read NASA coefficients
!----------------------
SUBROUTINE read_NASA_coeffs(nm,molname_loc,params)

USE mod_param
USE mod_cst
USE mod_chem, ONLY : P0, fcoeff, fcoeffnine, nparam, lname

IMPLICIT NONE

INTEGER          :: nm
CHARACTER(lname) :: molname_loc(nm)
REAL             :: params(nm,nparam)

!Local Variables
INTEGER   :: i
REAL      :: Tmin,Tmax,Tmean
LOGICAL   :: mol_inc
CHARACTER(lname) :: name_loc

!Set standard pressure
P0 = 1013250./unitP

params = 0.

DO i =1,nm

  mol_inc = .False.

  OPEN(1,file=fcoeff)
  fileloop: DO

    READ(1,*,END=1) name_loc,Tmin,Tmax,Tmean

    IF (molname_loc(i) == name_loc) THEN

      mol_inc = .True.

      params(i,1) = Tmin
      params(i,2) = Tmax
      params(i,3) = Tmean

      ! First two parameters are zero
      READ(1,*) params(i,6:12)
      READ(1,*) params(i,15:21)

      EXIT fileloop

    END IF

  END DO fileloop

1  CLOSE(1)

  !If not found in 7-coeff file, check 9-coeff file
  IF (.not. mol_inc) THEN
    OPEN(1,file=fcoeffnine)
    fileloop2: DO

      READ(1,*,END=2) name_loc,Tmin,Tmax,Tmean

      IF (molname_loc(i) == name_loc) THEN

        mol_inc = .True.

        params(i,1) = Tmin
        params(i,2) = Tmax
        params(i,3) = Tmean

        READ(1,*) params(i,4:12)
        READ(1,*) params(i,13:21)

        EXIT fileloop2

      END IF

    END DO fileloop2

2  CLOSE(1)
  END IF


  !Check molecule is included in NASA coeff file
  IF (.not. mol_inc) THEN
    WRITE(*,headrfmt) 'Error init_chem'
    WRITE(*,'(A)') 'Warning: Thermochemical data for '
    WRITE(*,chrpmfmt) 'molecule: ', molname_loc(i)
    WRITE(*,'(A)') 'not found.Check the NASA coefficients files.'
    WRITE(*,separfmt)
    STOP
  END IF


END DO

END SUBROUTINE read_NASA_coeffs

!----------------------
!Map elements in the molecules         : mol_ele
!Calculate molecule diameter           : dmol
!Calculate molecule mass               : molmass
!Calculate thermal diffusion parameter : therm_diff
!Caluclate the diffusion volumes       : sigma_v
!----------------------
SUBROUTINE get_molmass_dmol

USE mod_cst
USE mod_param
USE mod_chem

IMPLICIT NONE

INTEGER :: i,i2,ic,iel,e1,e2,elefactor1,elefactor2,lname_loc
REAL    :: molmass_loc,Vmol, sigma_v_loc
CHARACTER(lname)  :: name_loc,name_loc2
LOGICAL      :: skip_loop,cond

ALLOCATE(mol_ele(nmol,0:nele)) !0th index is for electron excess
ALLOCATE(molmass(nmol))
ALLOCATE(dmol(nmol))
ALLOCATE(molphase(nmol))
ALLOCATE(therm_diff(nmol))
ALLOCATE(sigma_v(nmol))

mol_ele = 0
molmass = 0.0
dmol = 0.0
therm_diff = 0.0
sigma_v = 0.0


molphase(:) = ' '

! compute molecules weight
DO i=1,nmol

  ! initialise all phases of molecules to gas phase
  molphase(i) = 'g'

  molmass_loc = 0.
  sigma_v_loc = 0.
  Vmol = 0.

  name_loc = molname(i)

  DO ic = 1,lname

    ! if a '/' character appear in the name of the molecule, this molecule is a condensate
    IF (name_loc(ic:ic) == '(') THEN
      molphase(i) = 'c'
      EXIT
    END IF

    IF (name_loc(ic:ic) == '-') THEN
      EXIT
    ENDIF

    ! Count electron excess (positive ions have negative excess)
    IF (name_loc(ic:ic) == '+') THEN
      mol_ele(i,0) = mol_ele(i,0) - 1
    ELSE IF (name_loc(ic:ic) == '_') THEN
      mol_ele(i,0) = mol_ele(i,0) + 1
    END IF

    skip_loop=.false.

    DO iel = 1,nele

      ! if 2-first characters are an element, count the number of element
      ! and increment the mass
      IF (ic+1.le.lname) THEN
        IF (name_loc(ic:ic+1)==elename(iel)) THEN
          skip_loop = .true.
          READ(name_loc(ic+2:ic+2),*,IOSTAT=e1) elefactor1
          ! if next character is an integer
          IF (e1==0) THEN
            READ(name_loc(ic+2:ic+3),*,IOSTAT=e2) elefactor2
            !if second next character is also an integer
            IF (e2==0) THEN

              molmass_loc = molmass_loc + elefactor2*elemass(iel)
              mol_ele(i,iel) = mol_ele(i,iel)+elefactor2
              Vmol = Vmol + 4./3.*pi*rel(iel)**3*elefactor2
              sigma_v_loc = sigma_v_loc + elsig(iel)*elefactor2

            ELSE

              molmass_loc = molmass_loc + elefactor1*elemass(iel)
              mol_ele(i,iel) = mol_ele(i,iel)+elefactor1
              Vmol = Vmol + 4./3.*pi*rel(iel)**3*elefactor1
              sigma_v_loc = sigma_v_loc + elsig(iel)*elefactor1

            END IF

          ELSE

            molmass_loc = molmass_loc + elemass(iel)
            mol_ele(i,iel) = mol_ele(i,iel)+1
            Vmol = Vmol + 4./3.*pi*rel(iel)**3
            sigma_v_loc = sigma_v_loc + elsig(iel)

          END IF
        END IF
      END IF
    END DO

    IF (.not. skip_loop) THEN

      DO iel = 1,nele
        ! if first character is an element, count the number of element
        ! and increment the mass
        IF (name_loc(ic:ic)==elename(iel)) THEN

          READ(name_loc(ic+1:ic+1),*,IOSTAT=e1) elefactor1

          IF (e1==0) THEN

            READ(name_loc(ic+1:ic+2),*,IOSTAT=e2) elefactor2

            IF (e2==0) THEN

              molmass_loc = molmass_loc + elefactor2*elemass(iel)
              mol_ele(i,iel) = mol_ele(i,iel)+elefactor2
              Vmol = Vmol + 4./3.*pi*rel(iel)**3*elefactor2
              sigma_v_loc = sigma_v_loc + elsig(iel)*elefactor2

            ELSE

              molmass_loc = molmass_loc + elefactor1*elemass(iel)
              mol_ele(i,iel) = mol_ele(i,iel)+elefactor1
              Vmol = Vmol + 4./3.*pi*rel(iel)**3*elefactor1
              sigma_v_loc = sigma_v_loc + elsig(iel)*elefactor1

            END IF

          ELSE

            molmass_loc = molmass_loc + elemass(iel)
            mol_ele(i,iel) = mol_ele(i,iel)+1
            Vmol = Vmol + 4./3.*pi*rel(iel)**3
            sigma_v_loc = sigma_v_loc + elsig(iel)
            

          END IF
        END IF
      END DO
    END IF
  END DO

  !Values for the thermal diffusion taken from Hobbs and Shorttle (2019)
   IF (name_loc == 'H' .or. name_loc == 'H2') THEN
      therm_diff(i) = -0.38
  ELSE
      therm_diff(i) = -0.25
  END IF

  molmass(i) = molmass_loc
  dmol(i) = 2.*(Vmol/pi*3./4.)**(1./3.)
  sigma_v(i) = sigma_v_loc

  !convert from picometer to code unit
  dmol(i) = dmol(i)*1E-10/unitL

  !Specific atmoic volumes taken from Poling et al. 2000, "The Properties of Gases and Liquids", McGraw Hill, Table 11.1
  if (name_loc=='H2') then
      sigma_v(i) = 6.12
  elseif (name_loc=='D2') then
      sigma_v(i) = 6.84
  elseif (name_loc=='N2') then
      sigma_v(i) = 18.5
  elseif (name_loc=='O2') then
      sigma_v(i) = 16.3
  elseif (name_loc=='CO') then
      sigma_v(i) = 18.0
  elseif (name_loc=='CO2') then
      sigma_v(i) = 26.9
  elseif (name_loc=='N2O') then
      sigma_v(i) = 35.9
  elseif (name_loc=='NH3') then
      sigma_v(i) = 20.7
  elseif (name_loc=='H2O') then
      sigma_v(i) = 13.1
  elseif (name_loc=='SF6') then
      sigma_v(i) = 71.3
  elseif (name_loc=='Cl2') then
      sigma_v(i) = 38.4
  elseif (name_loc=='Br2') then
      sigma_v(i) = 69.0
  elseif (name_loc=='SO2') then
      sigma_v(i) = 41.8
  endif
    
END DO

! count number of gas species
nmolg = 0

DO i = 1, nmol
  IF (molphase(i) == 'g') nmolg = nmolg + 1
END DO

! look if a condensed phase and a gas phase are present
IF (nmol .gt. nmolg) THEN
  DO i = 1,nmolg
    name_loc =molname(i)
    lname_loc = 0
    DO ic = 1,lname
      lname_loc = lname_loc + 1
      IF (name_loc(ic:ic)==' ') THEN
        lname_loc = lname_loc -1
      END IF
    END DO
    DO i2 = nmolg+1,nmol
      name_loc2 = molname(i2)
      IF (name_loc(1:lname_loc)==name_loc2(1:lname_loc) .and. name_loc2(lname_loc+1:lname_loc+1)=='(') THEN
        molphase(i2) = 'cg'
      END IF
    END DO
  END DO
END IF

! check if all the condensed species are at the end of the abundance file
cond = .false.
DO i=1,nmol
  IF (molphase(i)(1:1) == 'c') THEN
    cond = .true.
  END IF

  IF (molphase(i) == 'g' .and. cond) THEN
    IF (mype==cputerm) THEN
      WRITE(*,headrfmt) 'Error init_chem'
      WRITE(*,'(A)') 'All the condensed species should be at the end of the list'
      WRITE(*,*)
      WRITE(*,separfmt)
      STOP
    END IF
  END IF

END DO

END SUBROUTINE get_molmass_dmol

!----------------------
!Subroutine: calculate thermodynamic quantities of molecules
!             chemical potential, enthalpy, entropy and heat capacity
!----------------------
SUBROUTINE get_thermodynamics(nmol_loc,tloc,mp,mu,h,s,cp)

USE mod_chem, ONLY: nparam,tfreeze_eq,imol,molname

IMPLICIT NONE

INTEGER,INTENT(IN)  :: nmol_loc

REAL,INTENT(IN)  :: tloc
REAL,INTENT(IN)     :: mp(nmol_loc,nparam)
REAL,INTENT(OUT)    :: mu(nmol_loc)
REAL,INTENT(OUT)    :: h(nmol_loc)
REAL,INTENT(OUT)    :: s(nmol_loc)
REAL,INTENT(OUT)    :: cp(nmol_loc)

!Local variables
INTEGER :: i
INTEGER :: istart
REAL    :: tcalc

! Note :- tloc is the local temperature
!      :- tcalc is the temperature used for the thermodynamic calculation

DO i =1,nmol_loc

    IF (tloc>mp(i,3)) THEN
        istart = 4
    ELSE
        istart = 13
    END IF

    tcalc = tloc
    IF (tloc<tfreeze_eq) tcalc = tfreeze_eq
    !IF (tloc<mp(i,1)) tcalc = mp(i,1)
    !IF (tloc>mp(i,2)) tcalc = mp(i,2)

    !Compute chemical potential, enthalpy and heat capacity
    mu(i) = -(mp(i,istart)/2.)*tcalc**(-2.) + (mp(i,istart+1)/tcalc)*(1+log(tcalc)) &
        + mp(i,istart+2)*(1-log(tcalc)) - (mp(i,istart+3)/2.)*tcalc &
        - (mp(i,istart+4)/6.)*tcalc**2. - (mp(i,istart+5)/12.)*tcalc**3. &
        - (mp(i,istart+6)/20.)*tcalc**4. + (mp(i,istart+7)/tcalc) - mp(i,istart+8)
    h(i)  = -mp(i,istart)*tcalc**(-2.) + (mp(i,istart+1)/tcalc)*log(tcalc) &
        + mp(i,istart+2) + (mp(i,istart+3)/2.)*tcalc + (mp(i,istart+4)/3.)*tcalc**2. &
        + (mp(i,istart+5)/4.)*tcalc**3. + (mp(i,istart+6)/5.)*tcalc**4. + (mp(i,istart+7)/tcalc)
    s(i)  = (-mp(i,istart)/2.)*tcalc**(-2.) - mp(i,istart+1)*tcalc**(-1.) + mp(i,istart+2)*log(tcalc) &
        + mp(i,istart+3)*tcalc + (mp(i,istart+4)/2.)*tcalc**2. + (mp(i,istart+5)/3.)*tcalc**3. &
        + (mp(i,istart+6)/4.)*tcalc**4. + mp(i,istart+8)
    cp(i) = mp(i,istart)*tcalc**(-2.) + mp(i,istart+1)*tcalc**(-1.) + mp(i,istart+2) &
         + mp(i,istart+3)*tcalc + mp(i,istart+4)*tcalc**2. + mp(i,istart+5)*tcalc**3. &
         + mp(i,istart+6)*tcalc**4.

    ! Use the polynomial adjusted to the Saumon-Chabrier EOS for H2 for cp/R
    !if (molname(i)=='H2') then
    !  if (tloc>1000.0) then
    !    cp(i) = 0.*tcalc**(-2.) + 0.*tcalc**(-1.) + 3.10501126E+00 &
    !            + 3.87207292E-04*tcalc + 3.65811209E-07*tcalc**2. + (-2.74517202E-10)*tcalc**3. &
    !            + 6.56030779E-14*tcalc**4.
    !  else
    !    cp(i) = 0.*tcalc**(-2.) + 0.*tcalc**(-1.) + 2.46899055E+00 &
    !            + 6.53957408E-03*tcalc + (-1.47074304E-05)*tcalc**2. + 1.42487390E-08*tcalc**3. &
    !            + (-4.90075762E-12)*tcalc**4.
    !  endif
    !else
    ! cp(i) = mp(i,istart)*tcalc**(-2.) + mp(i,istart+1)*tcalc**(-1.) + mp(i,istart+2) &
    !     + mp(i,istart+3)*tcalc + mp(i,istart+4)*tcalc**2. + mp(i,istart+5)*tcalc**3. &
    !     + mp(i,istart+6)*tcalc**4.
    !
    !end if

END DO


END SUBROUTINE get_thermodynamics

!----------------------
!Subroutine: Calculate the enthalpy of formation of species
!
!----------------------
SUBROUTINE get_enthalpy_formation(idepth,nmol_loc,nmolg_loc,nele_loc,im,ie,hdrt,hform,tloc)

USE mod_chem, ONLY : nmol, nele,nparam,lname, iele, mol_ele

IMPLICIT NONE

INTEGER,INTENT(IN) :: idepth
INTEGER,INTENT(IN) :: nmol_loc
INTEGER,INTENT(IN) :: nmolg_loc
INTEGER,INTENT(IN) :: nele_loc
INTEGER,INTENT(IN) :: im(nmol)
INTEGER,INTENT(IN) :: ie(nele)

REAL              :: tloc
REAL, INTENT(IN)  :: hdrt(nmol)     ! enthalpy of molecules
REAL, INTENT(OUT) :: hform(nmol)

!Local variables
INTEGER :: i
INTEGER :: istart
REAL    :: href(nele)              ! enthalpy of reference elements
REAL    :: muref(nele)             ! chemical potential of reference elements
REAL    :: cpref(nele)             ! heat capacity of reference elements
REAL    :: sref(nele)              ! entropy of reference elements
REAL    :: refparams(nele,nParam)  ! NASA coefficients of reference elements

CHARACTER(lname) :: eleref_name(nele)

!Standard states of the elements
eleref_name  = (/'H2        ', &
                 'He        ', &
                 'C(gr)     ', &
                 'N2        ', &
                 'O2        ', &
                 'Na(cr)    ', &
                 'K(cr)     ', &
                 'Si(cr)    ', &
                 'Ar        ', &
                 'Ti(a)     ', &
                 'V(cr)     ', &
                 'S(cr1)    ', &
                 'Cl2       ', &
                 'Mg(cr)    ', &
                 'Al(cr)    ', &
                 'Ca(a)     ', &
                 'Fe(a)     ', &
                 'Cr(a)     ', &
                 'Li(cr)    ', &
                 'Cs(cr)    ', &
                 'Rb(cr)    ', &
                 'F2        ', &
                 'P(cr)     '/)

!Read in NASA coeffs of reference elements
CALL read_NASA_coeffs(nele,eleref_name,refParams)

!Get enthalpy of reference elements
CALL get_thermodynamics(nele,tloc,refParams,muref,href,sref,cpref)

DO i = 1, nmol_loc

    hform(i) = hdrt(i) - mol_ele(im(i),ie(iele('H')))*href(1)/2.  &
                       - mol_ele(im(i),ie(iele('He')))*href(2)    &
                       - mol_ele(im(i),ie(iele('C')))*href(3)     &
                       - mol_ele(im(i),ie(iele('N')))*href(4)/2.  &
                       - mol_ele(im(i),ie(iele('O')))*href(5)/2.  &
                       - mol_ele(im(i),ie(iele('Na')))*href(6)    &
                       - mol_ele(im(i),ie(iele('K')))*href(7)     &
                       - mol_ele(im(i),ie(iele('Si')))*href(8)    &
                       - mol_ele(im(i),ie(iele('Ar')))*href(9)    &
                       - mol_ele(im(i),ie(iele('Ti')))*href(10)   &
                       - mol_ele(im(i),ie(iele('V')))*href(11)    &
                       - mol_ele(im(i),ie(iele('S')))*href(12)    &
                       - mol_ele(im(i),ie(iele('Cl')))*href(13)/2.&
                       - mol_ele(im(i),ie(iele('Mg')))*href(14)   &
                       - mol_ele(im(i),ie(iele('Al')))*href(15)   &
                       - mol_ele(im(i),ie(iele('Ca')))*href(16)   &
                       - mol_ele(im(i),ie(iele('Fe')))*href(17)   &
                       - mol_ele(im(i),ie(iele('Cr')))*href(18)   &
                       - mol_ele(im(i),ie(iele('Li')))*href(19)   &
                       - mol_ele(im(i),ie(iele('Cs')))*href(20)   &
                       - mol_ele(im(i),ie(iele('Rb')))*href(21)   &
                       - mol_ele(im(i),ie(iele('F')))*href(22)/2. &
                       - mol_ele(im(i),ie(iele('P')))*href(1)
END DO



END SUBROUTINE get_enthalpy_formation


!----------------------
!Subroutine: calculate the specific heat capacity of the mixture
!            using the methodology of Gordon & McBride 1994
!----------------------
SUBROUTINE get_cp_mixture(idepth,nmol_loc,nmolg_loc,nele_loc,imol_loc,iele_loc,tloc,hdrt,cp0,lnlam)

USE mod_chem, ONLY : nmol, nele, mol_ele, cp_chem, cp_chem_reac
USE mod_cst, ONLY: unitE,unitM

IMPLICIT NONE

INTEGER,INTENT(IN) :: idepth, nmol_loc,nmolg_loc, nele_loc
INTEGER,INTENT(IN) :: imol_loc(nmol),iele_loc(nele)
REAL,INTENT(IN)    :: hdrt(nmol),cp0(nmol)
REAL,INTENT(IN)    :: lnlam(0:nele+1+nmol)
REAL,INTENT(IN)    :: tloc

!Local Variables
INTEGER :: nit,iel,i,rank,solverstat
INTEGER :: pivot(nmol+nele+1)

REAL    :: err, step, cpr
REAL    :: jacobi(nmol+1+nele,nmol+1+nele)
REAL    :: rhs(nmol+1+nele)
REAL    :: dlnlam(nmol+1+nele)
REAL    :: lnlamdt(nele+1+nmol)

REAL,PARAMETER    :: rgaskg = 8314.51 ! gas constant, J/(kg-mole)/K, Gordon and McBride 1994

! First compute frozen component
DO i=1,nmol_loc
  IF (i .le. nmolg_loc) THEN
    cpr = cpr + exp(lnlam(nele_loc+1+i))*cp0(imol_loc(i))
  ELSE
    cpr = cpr +lnlam(nele_loc+1+i)*cp0(imol_loc(i))
  END IF
END DO

! Compute reaction component
IF (cp_chem_reac) THEN

  ! Initialise some variables
  err = 1.
  nit = 0

  lnlamdt(1:nele_loc+1+nmol_loc) = 0.

  !Calculate derivitive dln(n)/dln(t) using NR iterations
  !Similar methodology to Gibbs minimisation
  DO WHILE (err>1E-6)

    ! Initialise
    err = 0.
    nit = nit + 1

    dlnlam = 0.
    jacobi = 0.
    rhs    = 0.

    !Fill RHS
    IF (nmol_loc==nmolg_loc) THEN
      DO iel =1,nele_loc
        rhs(iel) = - sum(exp(lnlam(nele_loc+2:nele_loc+1+nmolg_loc))* &
        mol_ele(imol_loc(1:nmolg_loc),iele_loc(iel))*lnlamdt(nele_loc+2:nele_loc+1+nmolg_loc))
        err = max(err,abs(rhs(iel))/(lnlam(nele_loc+1)*lnlamdt(nele_loc+1)))
      END DO
    ELSE
      DO iel =1,nele_loc
        rhs(iel) = - sum(exp(lnlam(nele_loc+2:nele_loc+1+nmolg_loc))*mol_ele(imol_loc(1:nmolg_loc),iele_loc(iel))*lnlamdt(nele_loc+2:nele_loc+1+nmolg_loc))  &
        - sum(lnlamdt(nele_loc+2+nmolg_loc:nele_loc+1+nmol_loc)*mol_ele(imol_loc(nmolg_loc+1:nmol_loc),iele_loc(iel)))
        err = max(err,abs(rhs(iel))/exp(lnlam(nele_loc+1))*lnlamdt(nele_loc+1))
      END DO
    END IF

    rhs(nele_loc+1) = exp(lnlam(nele_loc+1))*lnlamdt(nele_loc+1)-sum(exp(lnlam(nele_loc+2:nele_loc+1+nmolg_loc))*lnlamdt(nele_loc+2:nele_loc+1+nmolg_loc))
    err = max(err,abs(rhs(nele_loc+1))/exp(lnlam(nele_loc+1))*lnlamdt(nele_loc+1))

    DO i =1,nmol_loc
      !Gas species
      IF (i<=nmolg_loc) THEN
        rhs(nele_loc+1+i) = hdrt(imol_loc(i)) -(lnlamdt(nele_loc+1+i)+sum(mol_ele(imol_loc(i),iele_loc(1:nele_loc))*lnlamdt(1:nele_loc))-lnlamdt(nele_loc+1))
      !Condensed species
      ELSE
        rhs(nele_loc+1+i) =  hdrt(imol_loc(i)) -(sum(mol_ele(imol_loc(i),iele_loc(1:nele_loc))*lnlamdt(1:nele_loc)))
      END IF
      IF (abs(hdrt(imol_loc(i)))>0.) THEN
        err = max(err,abs(rhs(nele_loc+1+i))/hdrt(imol_loc(i)))
      END IF
    END DO

    !Fill Jacobian
    DO iel = 1,nele_loc
      DO i =1,nmol_loc
        IF (i<=nmolg_loc) THEN
          jacobi(iel,nele_loc+1+i) = mol_ele(imol_loc(i),iele_loc(iel))*exp(lnlam(nele_loc+1+i))
        ELSE
          jacobi(iel,nele_loc+1+i) = mol_ele(imol_loc(i),iele_loc(iel))
        END IF
      END DO
    END DO

    jacobi(nele_loc+1,nele_loc+1) = -exp(lnlam(nele_loc+1))

    DO i = nele_loc+2,nele_loc+1+nmolg_loc
      jacobi(nele_loc+1,i) = exp(lnlam(i))
    END DO

    DO i = 1,nmol_loc
      IF (i<=nmolg_loc) THEN
        jacobi(nele_loc+1+i,nele_loc+1+i) = 1.
        DO iel = 1,nele_loc
          jacobi(nele_loc+1+i,iel) = mol_ele(imol_loc(i),iele_loc(iel))
        END DO
        jacobi(nele_loc+1+i,nele_loc+1) = -1.
      ELSE
        DO iel = 1,nele_loc
          jacobi(nele_loc+1+i,iel) = mol_ele(imol_loc(i),iele_loc(iel))
        END DO
      END IF
    END DO

    rank = nele_loc + nmol_loc + 1

    !Call DGESV to find correction for dln(n)/dln(T)
    CALL DGESV(rank,1,jacobi(1:rank,1:rank),rank,pivot(1:rank),rhs(1:rank),rank,solverstat)

    !Copy solution
    dlnlam(1:rank) = rhs(1:rank)

    ! Check if the solution was successful
    IF (solverstat<0) THEN
      WRITE(*,'(A,1x,I1,1x,A)')'linsolve: argument',-solverstat,'for dgesv has an illegal value.'
      STOP
    ELSE IF (solverstat>0) THEN
      WRITE(*,'(A,1x,I4)')'linsolve: upper triangular matrix is singular! Diagonal element:',solverstat
      STOP
    END IF

    !Apply correction
    step = 1.0
    lnlamdt = lnlamdt + step*dlnlam

  END DO

  ! Add reaction term to frozen term
  DO i=1,nmol_loc
    IF (i .le. nmolg_loc) THEN
      cpr = cpr + exp(lnlam(nele_loc+1+i))*hdrt(imol_loc(i))*lnlamdt(nele_loc+1+i)
    ELSE
      cpr = cpr + hdrt(imol_loc(i))*lnlamdt(nele_loc+1+i)
    END IF
  END DO

END IF

!Specific heat capacity
!1. convert from unitless into J/(kg-mole)/K
cpr = cpr*rgaskg
!2. convert into code units
cp_chem(idepth) = cpr*1.0E4/unitE*unitM

END SUBROUTINE get_cp_mixture

!----------------------
!Subroutine: contains solver for
!minimization of gibbs free energy
!following the method of Gordon & McBride 1994
!----------------------
SUBROUTINE get_gibbs_abundances

USE mod_chem
USE mod_grid
USE mod_param

IMPLICIT NONE

INTEGER :: nmol_loc
INTEGER :: nele_loc
INTEGER :: idepth
INTEGER :: iel
INTEGER :: i, j
INTEGER :: nit
INTEGER :: icond_min
INTEGER :: iprint
INTEGER :: istart, iend
INTEGER :: imol_loc(nmol)
INTEGER :: iele_loc(nele)

INTEGER :: solverstat
INTEGER :: rank
INTEGER :: pivot(0:nmol+nele+1)

REAL :: step
REAL :: step1
REAL :: step2
REAL :: err
REAL :: tloc
REAL :: ploc
REAL :: mtot
REAL :: dGdn
REAL :: dGdn_min
REAL :: fh

!0th index is for charge conservation terms
REAL :: jacobi(0:nele+1+nmol,0:nele+1+nmol)
REAL :: rhs(0:nele+1+nmol)
REAL :: dlnlam(0:nele+1+nmol)
REAL :: lnlam(0:nele+1+nmol)
REAL :: lnlam_old(0:nele+1+nmol)
REAL :: mupot(nmol)
REAL :: cp0(nmol)
REAL :: hdrt(nmol)
REAL :: sdr(nmol)
REAL :: hform(nmol)
CHARACTER :: CR = CHAR(13)

REAL :: Aele_loc(nele)
LOGICAL :: arainout(nele)
INTEGER :: nmolg_loc,nmolc
LOGICAL :: mol_inc, ion_inc

! Charge conservation terms included for ions with
! mole fractions greater than fmin_ion
REAL :: fmin_ion = 0.0

EXTERNAL :: DGESV

mtot = 0.0
nele_loc = 0
iele_loc = 0

arainout = .False.

!Count and map number of elements
DO iel = 1,nele
  IF (sum(mol_ele(:,iel)) >0) THEN
    nele_loc = nele_loc + 1
    iele_loc(nele_loc) = iel
  END IF
END DO

!Calculate total elemental mass
DO iel =1,nele_loc
  mtot = mtot + Aele(iele_loc(iel))*elemass(iele_loc(iel))
END DO

Aele_loc = Aele

!Loop over idepth, calculate for each level in isolation
!Loop from 'bottom' to 'top'
DO idepth = ndepth+1,1,-1

    nmol_loc = 0
    imol_loc = 0

  !Initially only gas species
  j = 1
  DO i = 1,nmolg
      mol_inc = .True.
      DO iel = 1,nele
          IF ((mol_ele(i,iel) > 0) .and. (arainout(iel))) THEN
              mol_inc = .False.
          END IF
      END DO
      IF (mol_inc) THEN
          nmol_loc = nmol_loc + 1
          imol_loc(j) = i
          j = j + 1
      END IF
  END DO

  nmolg_loc = nmol_loc

  !Map which condensed species are available-ignore species which have rained out
  ! nmolc = 0
  ! DO i = nmolg+1,nmol
  !     mol_inc = .True.
  !     DO iel = 1, nele
  !          ! If contains element which has 'rained out', then ignore
  !          IF ((mol_ele(i,iel) > 0 .and. (arainout(iel)))) THEN
  !              mol_inc = .False.
  !          END IF
  !     END DO
  !     !Add to total and map
  !     IF (mol_inc) THEN
  !          nmolc = nmolc + 1
  !          imol_loc(j) = i
  !          j = j + 1
  !     END IF
  ! END DO

  tloc = ttf(idepth)
  ploc = ppf(idepth)

    !Get thermodynamic quantities of molecules
    CALL get_thermodynamics(nmol,tloc,molParam,mupot,hdrt,sdr,cp0)

    ! loop on the condensates
    icond_min  = 0
    nit        = 0

    !Initial estimates for number density, gas species, lagrange multipliers
    lnlam(0:nele_loc) = 0.
    lnlam(nele_loc+1) = log(0.1)
    lnlam(nele_loc+2:nele_loc+1+nmolg_loc) = log(0.1/nmol_loc)

    !Initial estimates for number density, condensed species
    IF (nmol_loc > nmolg_loc) THEN
      lnlam(nele_loc+2+nmolg_loc:nele_loc+1+nmol_loc) = 0.1/nmol_loc
    END IF

    !Loop over condensates, until no more condensates are added
    DO WHILE(icond_min >= 0)

      err = 1.0

      IF (icond_min>0) THEN
        lnlam(nele_loc+1+nmol_loc) = 0.1/nmol_loc
      END IF

      !Loop over chemistry, until converged is reached within accuracy
      DO WHILE (err>accuracy_chem)

        ! Initialise
        nit    = nit + 1
        err    = 0.0
        jacobi = 0.0
        rhs    = 0.0
        ion_inc = .False.

        !Test for non-convergence
        IF (nit>nitmax_chem) THEN
          WRITE(*,*) 'ERROR: Gibbs minimisation reached maximum iterations: ',nit
          STOP
        END IF

        ! Fill righthand side array
        ! 1.0 conservation of charge
        DO j = 1, nmolg_loc
          ! Current spec14s is an ion
          IF (mol_ele(imol_loc(j),0)/=0 .and. &
              exp(lnlam(nele_loc+1+j))/sum(exp(lnlam(nele_loc+2:nele_loc+1+nmolg_loc)))>fmin_ion) THEN
            rhs(0) = rhs(0) - exp(lnlam(nele_loc+1+j))*mol_ele(imol_loc(j),0)
            ion_inc = .True.
          END IF
        END DO
        err = max(err,abs(rhs(0)))

        ! 1.1 fill right hand side for elemental conservation; elemental mass-balance
        IF (nmol_loc==nmolg_loc) THEN
          DO iel =1,nele_loc
            rhs(iel) = - sum(exp(lnlam(nele_loc+2:nele_loc+1+nmolg_loc))*mol_ele(imol_loc(1:nmolg_loc),iele_loc(iel))) &
              + Aele_loc(iele_loc(iel))/mtot
            err = max(err,abs(rhs(iel))/(Aele_loc(iele_loc(iel))/mtot))
          END DO
        ELSE
          ! conservation of elements log(n) for gases just n for condensates
          DO iel =1,nele_loc
            rhs(iel) = - sum(exp(lnlam(nele_loc+2:nele_loc+1+nmolg_loc))*mol_ele(imol_loc(1:nmolg_loc),iele_loc(iel)))  &
              - sum(lnlam(nele_loc+2+nmolg_loc:nele_loc+1+nmol_loc)*mol_ele(imol_loc(nmolg_loc+1:nmol_loc),iele_loc(iel))) &
              + Aele_loc(iele_loc(iel))/mtot
            err = max(err,abs(rhs(iel))/(Aele_loc(iele_loc(iel))/mtot))
          END DO
        END IF

        ! 1.2  conservation of total gas density
        rhs(nele_loc+1) = exp(lnlam(nele_loc+1))-sum(exp(lnlam(nele_loc+2:nele_loc+1+nmolg_loc)))
        err = max(err,abs(rhs(nele_loc+1))/exp(lnlam(nele_loc+1)))

        ! 1.3 molecule equations
        ! Equation 2.10, Gordon & McBride 1994
        DO i =1,nmol_loc
          !Gas species
          IF (i<=nmolg_loc) THEN
            rhs(nele_loc+1+i) = -mupot(imol_loc(i)) - lnlam(nele_loc+1+i)- log(ploc/P0) +lnlam(nele_loc+1) &
              - sum(lnlam(1:nele_loc)*mol_ele(imol_loc(i),iele_loc(1:nele_loc)))
            IF (mol_ele(imol_loc(i),0)/=0 .and. &
              exp(lnlam(nele_loc+1+i))/sum(exp(lnlam(nele_loc+2:nele_loc+1+nmolg_loc)))>fmin_ion) THEN
                rhs(nele_loc+1+i) = rhs(nele_loc+1+i) - lnlam(0)*mol_ele(imol_loc(i),0)
            END IF
            err = max(err,abs(rhs(nele_loc+1+i))/abs(mupot(imol_loc(i))+log(ploc/P0)))
           !Condensed species
          ELSE
            rhs(nele_loc+1+i) = -mupot(imol_loc(i))  &
              -sum(lnlam(1:nele_loc)*mol_ele(imol_loc(i),iele_loc(1:nele_loc)))
            err = max(err,abs(rhs(nele_loc+1+i))/abs(mupot(imol_loc(i))))
          END IF
        END DO

        ! 2. fill the jacobi matrix
        ! 2.0 Charge conservation
        DO i = 1,nmolg_loc
          IF (mol_ele(imol_loc(i),0)/=0 .and. &
              exp(lnlam(nele_loc+1+i))/sum(exp(lnlam(nele_loc+2:nele_loc+1+nmolg_loc)))>fmin_ion) THEN
          jacobi(0,nele_loc+1+i) = mol_ele(imol_loc(i),0)*exp(lnlam(nele_loc+1+i))
          END IF
        END DO

        ! 2.1 terms of the elemental conservation
        DO iel = 1,nele_loc
          DO i =1,nmol_loc
            IF (i<=nmolg_loc) THEN
              jacobi(iel,nele_loc+1+i) = mol_ele(imol_loc(i),iele_loc(iel))*exp(lnlam(nele_loc+1+i))
            ELSE
              jacobi(iel,nele_loc+1+i) = mol_ele(imol_loc(i),iele_loc(iel))
            END IF
          END DO
        END DO

        !!! 2.2 terms of the gas density conservation
        jacobi(nele_loc+1,nele_loc+1) = -exp(lnlam(nele_loc+1))

        DO i = nele_loc+2,nele_loc+1+nmolg_loc
          jacobi(nele_loc+1,i) = exp(lnlam(i))
        END DO

        !!! 2.3 terms of the molecules equations
        DO i = 1,nmol_loc
          IF (i<=nmolg_loc) THEN
            jacobi(nele_loc+1+i,nele_loc+1+i) = 1.
            DO iel = 1,nele_loc
              jacobi(nele_loc+1+i,iel) = mol_ele(imol_loc(i),iele_loc(iel))
            END DO
            jacobi(nele_loc+1+i,nele_loc+1) = -1.
          ELSE
            DO iel = 1,nele_loc
              jacobi(nele_loc+1+i,iel) = mol_ele(imol_loc(i),iele_loc(iel))
            END DO
          END IF
          ! If species is an ion
          IF (mol_ele(imol_loc(i),0)/=0 .and. &
              exp(lnlam(nele_loc+1+i))/sum(exp(lnlam(nele_loc+2:nele_loc+1+nmolg_loc)))>fmin_ion) THEN
            jacobi(nele_loc+1+i,0) = mol_ele(imol_loc(i),0)
          END IF
        END DO

        ! Calculate rank of arrays
        ! Depends on presence of ions with mole fractions > fmin_ion
        IF (ion_inc) THEN
          rank = nele_loc + nmol_loc + 2
          istart = 0
          iend   = rank-1
        ELSE
          rank = nele_loc + nmol_loc + 1
          istart = 1
          iend = rank
        END IF

        !Call linsolve to calculate change in number density of each species, with constraint
        !of elemental mass-balance and total density
        CALL DGESV(rank,1,jacobi(istart:iend,istart:iend),rank,pivot(istart:iend),rhs(istart:iend),rank,solverstat)

        !Copy solution
        dlnlam(istart:iend) = rhs(istart:iend)

        ! Check if the solution was successful
        IF (solverstat<0) THEN
          WRITE(*,'(A,1x,I1,1x,A)')'linsolve: argument',-solverstat,'for dgesv has an illegal value.'
          STOP
        ELSE IF (solverstat>0) THEN
          WRITE(*,'(A,1x,I4)')'linsolve: upper triangular matrix is singular! Diagonal element:',solverstat
          STOP
        END IF

        !Set the step size to aid convergence
        !See Gordon & McBride 1994, section 3.3
        step1 = 1.0
        step2 = 1.0

        step1 = 5.*abs(dlnlam(nele_loc+1))

        DO i=1,nmolg_loc
          IF ((lnlam(nele_loc+1+i)-lnlam(nele_loc+1))>-46) THEN
            step1 = max(step1,abs(dlnlam(nele_loc+1+i)))
          ELSE IF (dlnlam(nele_loc+1+i)>0. .and.(lnlam(nele_loc+1+i)-lnlam(nele_loc+1))<=-46) then
            step2 = min(step2,abs((-lnlam(nele_loc+1+i)+lnlam(nele_loc+1)-23)/(dlnlam(nele_loc+1+i)-dlnlam(nele_loc+1))))
          END IF
        END DO

        step1 = min(chem_conv/step1,1.)
        step  = min(step1,step2)

        !Update number densities using step size and dlnlam from linsolve
        lnlam(nele_loc+2:nele_loc+1+nmol_loc) = lnlam(nele_loc+2:nele_loc+1+nmol_loc) &
          +step*dlnlam(nele_loc+2:nele_loc+1+nmol_loc)
        lnlam(1:nele_loc) = lnlam(1:nele_loc)+step*dlnlam(1:nele_loc)
        lnlam(nele_loc+1) = lnlam(nele_loc+1)+step*dlnlam(nele_loc+1)
        IF (ion_inc) lnlam(0) = lnlam(0) + step*dlnlam(0)

        !Save lnlam to initialise next cell
        lnlam_old = lnlam

        !Print convergence info screen
        IF (mype == cputerm .and. print_chem) THEN
          WRITE(*,'(A)',advance='no') CR
          WRITE(*,1000,advance='no') idepth,nit,maxval(abs(dlnlam(1:nmol_loc)/(lnlam(1:nmol_loc)+1E-40))),err
        ENDIF
1000    FORMAT ((5x,'idepth: ',i5,5x,'Iter. Chem.: ',i5,5x,' dln/ln: ',es10.1,5x,' err.: ',es10.1))

      END DO

      !Update mixing ratios
      Amolf(idepth,:) = 0.

      IF (nmol_loc == nmolg_loc) THEN
        Amolf(idepth,imol_loc(1:nmolg_loc)) = exp(lnlam(nele_loc+2:nele_loc+1+nmolg_loc)) &
         /sum(exp(lnlam(nele_loc+2:nele_loc+1+nmolg_loc)))
      ELSE
        Amolf(idepth,imol_loc(1:nmolg_loc)) = exp(lnlam(nele_loc+2:nele_loc+1+nmolg_loc)) &
         /(sum(exp(lnlam(nele_loc+2:nele_loc+1+nmolg_loc))))
        Amolf(idepth,imol_loc(nmolg_loc+1:nmol_loc)) = lnlam(nele_loc+2+nmolg_loc:nele_loc+1+nmol_loc)&
          /(sum(exp(lnlam(nele_loc+2:nele_loc+1+nmolg_loc))))
      END IF

      !Test for addition/removal of condensates
      !See Gordon & McBride 1994, section 3.4
      dGdn  = 0.0
      dGdn_min  = -1.0
      icond_min = -1

      IF (nmol>nmolg) THEN
        DO i = nmolg+1,nmol
          dGdn = mupot(i) + sum(mol_ele(i,iele_loc(1:nele_loc))*lnlam(1:nele_loc))
          IF (dGdn<dGdn_min .and. abs(dGdn/mupot(i))>1.0E-13 .and. (tloc<=molParam(i,2)) .and. (tloc>=molParam(i,1)) ) THEN

            mol_inc = .true.

            do j = 1,nele
              if((.not.(mol_ele(i,j) ==0)) .and. arainout(j)) then
                  mol_inc = .false.
              endif

            enddo

            if (mol_inc) then
              dGdn_min = dGdn
              icond_min = i
            endif

          END IF
        END DO

        !If Gibbs energy can be reduced, add condensed species
        IF (icond_min>=0) THEN
          nmol_loc = nmol_loc + 1
          imol_loc(nmol_loc) = icond_min
        END IF

        !If abundance of condensate is negative, removed from mixture
        DO i = nmol_loc, nmolg_loc+1,-1
          IF (amolf(idepth,imol_loc(i))<0.) THEN
            amolf(idepth,imol_loc(i)) = 0.0
            imol_loc(i) = imol_loc(nmol_loc)
            nmol_loc = nmol_loc - 1
            icond_min = 1
          END IF
        END DO

      END IF

    END DO !End chemistry loop

    IF (mype == cputerm .and. print_chem) THEN
      WRITE(*,*)
    END IF

    !Reduce elemental abundance due to condensation
    IF (rainout) THEN

      ! Calculate the number fraction of hydrogen atoms (number density of H atoms/number density of gas)
      fh = sum(mol_ele(imol_loc(1:nmol_loc),iele_loc(iele('H')))*Amolf(idepth,imol_loc(1:nmol_loc)))

      DO i = 1,nele_loc
        DO j = nmolg_loc+1,nmol_loc
          ! Reduce element abundance (element abundances are N_i/N_H
          Aele_loc(iele_loc(i)) = Aele_loc(iele_loc(i)) - mol_ele(imol_loc(j),iele_loc(i))*Amolf(idepth,imol_loc(j))/fh
        END DO
        ! If negative abundance, flag element as rained out
        IF (Aele_loc(iele_loc(i)) < rainout_min_ele) THEN
          arainout(iele_loc(i)) = .True.
          Aele_loc(iele_loc(i)) = 0.
          IF ((mype == 0) .AND. print_chem) THEN
            WRITE(*,chrpmfmt) 'Rainout: ',elename(iele_loc(i))
          END IF
          ! Record the rainout level of a given species
          IF (.NOT. freeze_rainout) THEN
            IF (elename(iele_loc(i)) == elename_freeze) THEN
              rainout_level = idepth
              nelename_freeze = iele_loc(i)
            END IF
          END IF
        END IF
      END DO

      IF (freeze_rainout) THEN
        IF (idepth == rainout_level) THEN
          WRITE(*,'("Manually setting rainout level of ",A2,1x,A,1x,i3)') elename_freeze, 'at model level', idepth
          WRITE(*,*)
          arainout(nelename_freeze) = .True.
        END IF
      END IF

    END IF

    !do i = nmolg_loc+1,nmol_loc
    !  if (molname(imol_loc(i)) == 'H2O(s)' .or. molname(imol_loc(i)) == 'H2O(l)' .or. &
    !    & molname(imol_loc(i)) == 'NH3(s)' .or. molname(imol_loc(i)) == 'NH3(l)') then
    !
    !    if (amolf(idepth,imol_loc(i)) /= 0.0) then
    !
    !      write (*,"(i3,3x,a6,3x,'abundance = ',es10.1)") idepth, molname(imol_loc(i)), &
    !        & amolf(idepth,imol_loc(i))

    !    endif
    !  endif

    !enddo

    iele_loc = 0
    nele_loc = 0
    !Recount and remap number of elements
    DO iel = 1,nele
        IF ((sum(mol_ele(:,iel)) > 0) .and. (.not. arainout(iel))) THEN
            nele_loc = nele_loc + 1
            iele_loc(nele_loc) = iel
        END IF
    END DO

    !Renormalise elemental abundances
    !Aele_loc(iele_loc(1:nele_loc)) = Aele_loc(iele_loc(1:nele_loc))/sum(Aele_loc(iele_loc(1:nele_loc)))

    mtot = 0.
    DO iel =1,nele_loc
        mtot = mtot + Aele_loc(iele_loc(iel))*elemass(iele_loc(iel))
    END DO

    !Calculate Cp of layer
    IF (calc_cp_chem) CALL get_cp_mixture(idepth,nmol_loc,nmolg_loc,nele_loc,imol_loc,iele_loc,tloc,hdrt,cp0,lnlam)

    !Calculate enthalpy of formation for the species
    IF (calc_hform) THEN
        CALL get_enthalpy_formation(idepth,nmol_loc,nmolg_loc,nele_loc,imol_loc,iele_loc,hdrt,hform,tloc)
    END IF

    !Calculate mean molecule weight
    mmw(idepth) = sum(Amolf(idepth,:)*molmass(:))

    !Copy final element abundances
    Aele_lev(idepth,:) = Aele_loc(:)

  END DO !End depth loop

  !If condensates present, print to screen
  IF (nmol.gt. nmolg_loc) THEN
    IF (mype == cputerm .and. print_chem) THEN
      WRITE(*,'(a)',advance='no') 'Cond:     '
    END IF

    iprint = 0
    DO i=nmolg+1,nmol
      IF (sum(amolf(:,i))>0. ) THEN
        IF (mype == cputerm .and. print_chem) THEN
          WRITE(*,'(a,1x)',advance='no') molname(i)
          iprint = iprint + 1
          IF (iprint ==7) THEN
            WRITE(*,*)
            WRITE(*,'(a)',advance='no') '          '
            iprint = 0
          END IF
        END IF
      END IF
    END DO
  END IF

END SUBROUTINE get_gibbs_abundances

!----------------------
! Analytical formula of Burrows and Sharp 1999
!----------------------
SUBROUTINE get_analytical_abundances

USE mod_grid
USE mod_chem
USE mod_param
USE mod_cst

IMPLICIT NONE

!local variables
INTEGER :: i,idepth
REAL :: Ael_H,Ael_He,Ael_C,Ael_O,Ael_Si,Ael_N,Ael_V,Ael_Ti,Ael_Na,Ael_K,coefft1,coefft2

!partial pressure profile of H2
REAL :: ppH2(ndepth+1)

! Number of equilibrium constants and coefficients in equilibrium constant
INTEGER,PARAMETER :: neqconst=2,neqcoeff=5

! equilibrium constant
REAL :: K(neqconst,ndepth+1)

! Gas constant in cal/(mol*K) used for equilibrium constants
REAL,PARAMETER :: R_cal=1.9858775

! Average number of oxygen atoms removed per silicon atom
REAL,PARAMETER :: xSi=3.28

! Temperature below which oxygen is removed by silicon atoms
REAL, PARAMETER :: T_remove_O=1.!1500
! Characteristic transition scale in T over which O is removed
REAL,PARAMETER :: T_trans_scale_remove_O=0.0
REAL :: rmO(ndepth+1)

! Temperature below which H2O condenses
REAL :: tcondH2O(ndepth+1)

! Temperature below which NH3 condenses
REAL :: tcondNH3(ndepth+1)

! Assign value to equilibrium coefficients:
REAL,DIMENSION(neqcoeff),PARAMETER :: &
eqcoeff_1 = (/ 1.106131d+6, -5.6895d+4, 62.565d+0, -5.81396d-4, &
2.346515d-8 /), &
eqcoeff_2 = (/ 8.16413d+5, -2.9109d+4, 58.5878d+0, -7.8284d-4, &
4.729048d-8 /)

! Assign value to condensation line parameters
REAL,DIMENSION(4),PARAMETER :: &
param_H2O = (/ 38.84, -3.83, -3.93, -0.20 /), &
param_NH3 = (/ 68.02, -6.31, -6.19, 0.0 /)

REAL :: &
T_crit_TiO = 1800, &!1, &!
P_crit_TiO = 1.0E+07, &!1.0E+10, &!
T_crit_VO  = 1600, &!1, &!
P_crit_VO  = 1.0E+07, &!1.0E+10, &!
T_crit_Na  = 1000, &!1, &!
T_crit_K   = 1100!1!

! Characteristic transition scale in P,T over which the abundance of TiO
! and VO change
REAL :: &
T_trans_scale_TiO = 10.0, &
P_trans_scale_TiO = 0.0, &
T_trans_scale_VO  = 10.0, &
P_trans_scale_VO  = 0.0, &
T_trans_scale_Na  = 10.0, &
T_trans_scale_K   = 10.0

coefft1 = 30.
coefft2 = 30.

! Calculate equilibrium constants
K(1,:) = exp((eqcoeff_1(1)/ttf + eqcoeff_1(2) + eqcoeff_1(3)*ttf + &
& eqcoeff_1(4)*ttf**2 + eqcoeff_1(5)*ttf**3)/(R_cal*ttf))
K(2,:) = exp((eqcoeff_2(1)/ttf + eqcoeff_2(2) + eqcoeff_2(3)*ttf + &
& eqcoeff_2(4)*ttf**2 + eqcoeff_2(5)*ttf**3)/(R_cal*ttf))

! Calculate condensation temperatures
IF (cond_NH3)  THEN
  tcondNH3 = 1.0E4/(param_NH3(1)+param_NH3(2)*log10(ppf*unitP/1.0E6)+ &		     !pressure in bars
    param_NH3(3)*MdH+param_NH3(4)*MdH*log10(ppf*unitP/1.0E6))         !pressure in bars
ELSE
  tcondNH3 = 0.
END IF

IF (cond_H2O) THEN
  tcondH2O = 1.0E4/(param_H2O(1)+param_H2O(2)*log10(ppf*unitP/1.0E6)+ &		     !pressure in bars
    param_H2O(3)*MdH+param_H2O(4)*MdH*log10(ppf*unitP/1.0E6))         !pressure in bars
ELSE
  tcondH2O = 0.
END IF

! compute partial pressures and store them in mmr before conversion to mmr

!elemental abundances relative to H
DO i = 1,nele
  IF (elename(i) == 'H') THEN
    Ael_H = Aele(i)
  END IF
END DO
DO i = 1,nele
  IF (elename(i) == 'He') THEN
    Ael_He = Aele(i)/Ael_H
  ELSE IF (elename(i) == 'C') THEN
    Ael_C = Aele(i)/Ael_H
  ELSE IF (elename(i) == 'N') THEN
    Ael_N = Aele(i)/Ael_H
  ELSE IF (elename(i) == 'O') THEN
    Ael_O = Aele(i)/Ael_H
  ELSE IF (elename(i) == 'Na') THEN
    Ael_Na = Aele(i)/Ael_H
  ELSE IF (elename(i) == 'K') THEN
    Ael_K = Aele(i)/Ael_H
  ELSE IF (elename(i) == 'Si') THEN
    Ael_Si = Aele(i)/Ael_H
  ELSE IF (elename(i) == 'Ti') THEN
    Ael_Ti = Aele(i)/Ael_H
  ELSE IF (elename(i) == 'V') THEN
    Ael_V = Aele(i)/Ael_H
  END IF
END DO

Ael_H = 1.

!compute H2 partial pressure assuming all H is in H2
ppH2 = ppf*(Ael_H/2.0)/(Ael_H/2.0+Ael_He)

!looks if oxygen is removed by Silicon atoms
rmO = 1.0/ &
(exp((ttf - T_remove_O)/T_trans_scale_remove_O) + 1.0)

DO i = 1,nmol
  IF (Acst(i) .ne. 0.) THEN
    Amolf(:,i) = Acst(i)
  ELSE
    IF (molname(i)=='H2') THEN
      Amolf(:,i) = (Ael_H/2.0)/(Ael_H/2.0+Ael_He)
    ELSE IF (molname(i)=='He') THEN
      Amolf(:,i) = Ael_He/(Ael_H/2.0+Ael_He)
    ELSE IF (molname(i)=='CO') THEN
      ! partial pressure relative to H2
      Amolf(:,i) = Ael_C + (Ael_O - xSi*Ael_Si*rmO) + &
        (ppH2/atm)**2/(2*K(1,:)) - sqrt((Ael_C + (Ael_O - &
        xSi*Ael_Si*rmO) + (ppH2/atm)**2/(2*K(1,:)))**2 - &
        4.*Ael_C*(Ael_O - xSi*Ael_Si*rmO))

      ! convert to partial pressures
      Amolf(:,i) = Amolf(:,i)*ppH2(:)/ppf(:)

      ! CO disappears if H2O is removed by condensation
      IF (cond_H2O) THEN
        Amolf(:,i) = Amolf(:,i)/(1+exp(coefft1*(tcondH2O-ttf)/tcondH2O))
      END IF

    ELSE IF (molname(i)=='CH4') THEN

    ! partial pressure relative to H2
    Amolf(:,i) = Ael_C + (Ael_O - xSi*Ael_Si*rmO) + &
      (ppH2/atm)**2/(2*K(1,:)) - sqrt((Ael_C + (Ael_O - &
      xSi*Ael_Si*rmO) + (ppH2/atm)**2/(2*K(1,:)))**2 - &
      4.*Ael_C*(Ael_O - xSi*Ael_Si*rmO))
    Amolf(:,i) = 2*Ael_C - Amolf(:,i)

    ! convert to partial pressures
    Amolf(:,i) = Amolf(:,i)*ppH2(:)/ppf(:)

    ! CO disappears if H2O is removed by condensation
    IF (cond_H2O) THEN
      WHERE (ttf <= tcondH2O)
        Amolf(:,i) = Ael_C/(Ael_H/2.0+Ael_He) + &
          (Amolf(:,i) - Ael_C/(Ael_H/2.0+Ael_He))/(1+exp(coefft1*(tcondH2O-ttf)/tcondH2O))
      END WHERE
    END IF

  ELSE IF (molname(i).eq.'H2O') THEN

    ! partial pressure relative to H2
    Amolf(:,i) = Ael_C + (Ael_O - xSi*Ael_Si*rmO) + &
      (ppH2/atm)**2/(2*K(1,:)) - sqrt((Ael_C + (Ael_O - &
      xSi*Ael_Si*rmO) + (ppH2/atm)**2/(2*K(1,:)))**2 - &
      4.*Ael_C*(Ael_O - xSi*Ael_Si*rmO))
    Amolf(:,i) = 2*(Ael_O - xSi*Ael_Si*rmO) - Amolf(:,i)

    ! convert to partial pressures
    Amolf(:,i) = Amolf(:,i)*ppH2(:)/ppf(:)

    ! H2O condensation
    IF (cond_H2O) THEN
      Amolf(:,i) = Amolf(:,i)/(1+exp(coefft1*(tcondH2O-ttf)/tcondH2O))
    END IF

  ELSE IF (molname(i).eq.'N2') THEN

    ! partial pressure relative to H2
    WHERE (Ael_N/((ppH2/atm)**2/(8.0*K(2,:)))< 1.0E-10)
      Amolf(:,i) = 0.!Ael_N - sqrt(2*Ael_N*(ppH2/atm)**2/(8.0*K(2,:)) )
    ELSE WHERE
      Amolf(:,i) = Ael_N + (ppH2/atm)**2/(8.0*K(2,:)) - &
        sqrt((Ael_N + (ppH2/atm)**2/(8.0*K(2,:)))**2 - Ael_N**2)
    END WHERE

    ! convert to partial pressures
    Amolf(:,i) = Amolf(:,i)*ppH2(:)/ppf(:)

    ! N2 disappears if NH3 is removed by condensation
    IF (cond_NH3) THEN
      Amolf(:,i) = Amolf(:,i)/(1+exp(coefft2*(tcondNH3-ttf)/tcondNH3))
    END IF

  ELSE IF (molname(i).eq.'NH3') THEN

    WHERE (Ael_N/((ppH2/atm)**2/(8.0*K(2,:))) < 1.0E-10)
      Amolf(:,i) = 2*Ael_N
    ELSE WHERE
      ! partial pressure relative to H2
      Amolf(:,i) = 2*(sqrt((ppH2/atm)**2/(8.0*K(2,:))*(2*Ael_N + &
        (ppH2/atm)**2/(8.0*K(2,:)))) - (ppH2/atm)**2/(8.0*K(2,:)))
    END WHERE

    ! convert to abundances
    Amolf(:,i) = Amolf(:,i)*ppH2(:)/ppf(:)

    ! NH3 condensation
    IF (cond_NH3) THEN
      Amolf(:,i) = Amolf(:,i)/(1+exp(coefft2*(tcondNH3-ttf)/tcondNH3))
    END IF

  ELSE IF (molname(i).eq.'Na') THEN

    DO idepth=1,ndepth+1
      Amolf(idepth,i) = Ael_Na/(Ael_H/2.0+Ael_He)*1.0!/ &
        !(EXP(-(ttf(idepth) - T_crit_Na)/T_trans_scale_Na) + 1.0)
    END DO

  ELSE IF (molname(i).eq.'K') THEN

    DO idepth=1,ndepth+1
      Amolf(idepth,i) = Ael_K/(Ael_H/2.0+Ael_He)*1.0!/ &
        !(EXP(-(ttf(idepth) - T_crit_K)/T_trans_scale_K) + 1.0)
    END DO

  ELSE IF (molname(i).eq.'Si') THEN

    DO idepth=1,ndepth+1
      Amolf(idepth,i) = Ael_Si/(Ael_H/2.0+Ael_He)
    END DO

  ELSE IF (molname(i).eq.'TiO') THEN

    DO idepth=1,ndepth+1
      Amolf(idepth,i) = Ael_Ti/(Ael_H/2.0+Ael_He)!*1.0/ &
        !(EXP(-(ttf(idepth) - T_crit_TiO)/T_trans_scale_TiO) + 1.0)
    END DO

  ELSE IF (molname(i).eq.'VO') THEN

    DO idepth=1,ndepth+1
      Amolf(idepth,i) = Ael_V/(Ael_H/2.0+Ael_He)!*1.0/ &
        !(EXP(-(ttf(idepth) - T_crit_VO)/T_trans_scale_VO) + 1.0)
    END DO

  ELSE IF (molname(i) .eq. 'H2O(SL)') THEN

    Amolf(:,i) = 0.

  ELSE IF (molname(i) .eq. 'NH3(SL)') THEN

    Amolf(:,i) = 0.

  ELSE

    IF (mype==cputerm) THEN
      WRITE(*,headrfmt) 'Warning get_Amol'
      WRITE(*,'(A)') 'Requested molecule does not have an analytic abundance inplemented and the specified Abundance Acst is 0.0'
      WRITE(*,'(A)') 'The molecule will be ignored by the code.'
      WRITE(*,intpmfmt) 'imol',i
      WRITE(*,chrpmfmt) 'molname',molname(i)
      WRITE(*,fltpmfmt) 'Acst',Acst(i)
      WRITE(*,*)
      WRITE(*,separfmt)
    END IF

  END IF

END IF
END DO

END SUBROUTINE get_analytical_abundances
