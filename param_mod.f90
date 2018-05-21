module param_mod
  use grid_mod
  implicit none
  !
  ! Time stepping and output parameters
  !
  character(50)   :: runname         ! runname
  real            :: dt              ! time step
  integer         :: tfinal          ! final real time in seconds
  integer         :: outcount = 0           ! number to start data files
!  logical         :: isBinary = .True.      ! flag for binary output
  logical         :: isBinary = .False.      ! flag for binary output

  !
  ! dimensional fluid parameters
  !
  real :: mu       ! poise 0.04 is vis of plasma, 0.01 is vis of H2O
  real :: rho      ! g/cm^3
  real :: tmax     ! characteristic time
  real :: pltdiam  ! characteristic length
  real :: dheight  ! domain height
  real :: umax     !
  real :: ilengthnd ! nondimensional injury length
  real :: rzL       ! location of left side of reaction zone 
  real :: rzR       ! location of right side of reaction zone 
  integer :: rzLi   ! index at left side of reaction zone 
  integer :: rzRi   ! index at right side of reaction zone 

  ! parameters to define the shape of alpha^2 function
  !
  real :: xmid     ! midpoint of x domain
  real :: xsteep   ! steepness of bump in x
  real :: ysteep   ! steepness of bump in y
  real :: xhw      ! half width of bump in x dir

  real :: xleft    ! left edge
  real :: xright   ! right edge
  real :: ytop     ! top edge
  real :: asqmax   ! max of alpha^2

  ! weight for sor
  !
  real :: omegaf 
  real :: omega


  !
  ! Scales of variables in cgs units
  !
  real :: xchar
  real :: tchar
  real :: uchar
  real :: pchar 
  real :: fchar 

  ! Reynolds number
  !
  real :: Re

!real:: asq(1:nxb+1,1:nyb+1,2)
  ! Biochem
  !
!  integer,parameter :: num_pl = 5  ! number of platelet species
!  integer,parameter :: num_se = 9  ! total number of SE-bound species
!  integer,parameter :: num_fp = 33 ! total number of fluid phase species
!  integer,parameter :: num_pb = 39 ! total number of platelet-bound species
  

!  real :: plinit(-1:nxb+2,-1:nyb+2)! initial condition for platelets
!  real :: seinit(1:nxb)            ! initial condition for SE-bound
!  real :: fpinit(-1:nxb+2,-1:nyb+2)! initial configuration for fluid phase
!  real :: pbinit(1:nxb,1:nyb      )! initial condition for platelet-bound

  ! Non-dimensional values for output

  !real :: ndpl(1:num_pl)
  !real :: ndse(1:num_se)
  !real :: ndfp(1:num_fp)
 ! real :: ndpb(1:num_pb)
!
  ! Create types for each 
  !
!    type platelet
!     character(len=20):: name                    ! species
!     real             :: value(-1:nxb+2,-1:nyb+2)! concentration with 2 ghost cells
!     real             :: bv                      ! left boundary value
!     real             :: k1(1:nxb,1:nyb)         ! for RK
!     real             :: k2(1:nxb,1:nyb)         ! for RK
!     real             :: k3(1:nxb,1:nyb)         ! for RK
!     real             :: k4(1:nxb,1:nyb)         ! for RK
!     real             :: temp(1:nxb,1:nyb)       ! for RK
!     real             :: diff(0:nxb+1,0:nyb+1)   ! diffusion coefficient
!     integer          :: eqn_num                 ! equation number
!  end type platelet

!  type sebound
!     character(len=20):: name        ! species
!     real             :: value(1:nxb)! concentration (cell centered)
!     real             :: k1(1:nxb)   ! for RK
!     real             :: k2(1:nxb)   ! for RK
!     real             :: k3(1:nxb)   ! for RK
!     real             :: k4(1:nxb)   ! for RK
!     real             :: temp(1:nxb) ! for RK
!     integer          :: eqn_num     ! equation number
!  end type sebound

!  type fluidphase
!     character(len=20):: name                    ! species
!     real             :: value(-1:nxb+2,-1:nyb+2)! concentration with 2 ghost cells
!     real             :: bv                      ! left boundary value
!     real             :: k1(1:nxb,1:nyb)         ! for RK
!     real             :: k2(1:nxb,1:nyb)         ! for RK
!     real             :: k3(1:nxb,1:nyb)         ! for RK
!     real             :: k4(1:nxb,1:nyb)         ! for RK
!     real             :: temp(1:nxb,1:nyb)       ! for RK
!     real             :: diff(0:nxb+1,0:nyb+1)   ! diffusion coefficient
!     integer          :: eqn_num                 ! equation number
!  end type fluidphase

!  type pltbound
!     character(len=20):: name                    ! species
!     real             :: value(1:nxb,1:nyb)      ! concentration
!     real             :: k1(1:nxb,1:nyb)         ! for RK
!     real             :: k2(1:nxb,1:nyb)         ! for RK
!     real             :: k3(1:nxb,1:nyb)         ! for RK
!     real             :: k4(1:nxb,1:nyb)         ! for RK
!     real             :: temp(1:nxb,1:nyb)       ! for RK
!     integer          :: eqn_num                 ! equation number
!  end type pltbound



  ! make array of all chemical species
  ! so can call by number and loop through all 
  !
 ! type(platelet)      :: pl(1:num_pl)
 ! type(sebound)       :: se(1:num_se)
 ! type(fluidphase)    :: fp(1:num_fp) 
 ! type(pltbound)      :: pb(1:num_pb)
  

! PLATELET PARAMS
!
!  real :: P0                 ! normal plasma concentration
!  real :: Pmaxb              ! max conc of platelets that can fit in fluid
!  real :: Pmaxse             ! max conc of platelets that can fit on se
!  real :: kcoh               ! cohesion rate 1/s
!  real :: krel               ! release rate 1/s
!  real :: kadh_max           ! constant for maximum
!  real :: kadh(1:nxb,1:nyb)  ! adhesion to se 1/s
!  real :: pbdryL(1:nyb)      ! upstream conc. of mobile unactivated platelet 

!  real :: DCP                         ! diff coeff
!  real :: DiffCoeff(0:nxb+1,0:nyb+1)  ! diff coeff for fluid phase species
!  real :: DiffCoeffP(0:nxb+1,0:nyb+1) ! diff coeff for platelets
  
 ! real :: rzheight

!
!   Non-Dimensional Parameters
!
 ! real :: NDadvect    ! non-dim coeff in front of advection term
 ! real :: NDdiffuse   ! non-dim coeff in front of diffusion term


    ! KINETIC RATE CONSTANTS
 !   real :: k7on,k7off
 !   real :: k2cat,k2min,k2plu
 !   real :: k3cat,k3min,k3plu
 !   real :: k8cat,k8min,k8plu,k8starm,k8starp
 !   real :: k9cat,k9min,k9plu

 !   real :: k1cat,k1min,k1plu
 !   real :: k12cat,k12min,k12plu
 !   real :: k14cat,k14min,k14plu
 !   real :: k18cat,k18min,k18plu

 !   real :: k9on,k9off
 !   real :: k10on,k10off
    ! FX on/off
 !   real :: k10zon,k10zoff
 !   real :: k5on,k5off
 !   real :: k8on,k8off
 !   real :: k2on,k2off

 !   real :: ke2on,ke2off !new thrombin binding sites
 !   real :: ke2gp1bon,ke2gp1boff !new thrombin binding sites gp1b


 !   real :: k5cat,k5min,k5plu
 !   real :: k13cat,k13min,k13plu
 !   real :: k6cat,k6min,k6plu
 !   real :: k15cat,k15min,k15plu
 !   real :: ktenmin,ktenplu
 !   real :: kpromin,kproplu
 !   real :: k4cat,k4min,k4plu
 !   real :: k7cat,k7min,k7plu

 !   real :: k17cat,k17min,k17plu
 !   real :: k16cat,k16min,k16plu
 !   real :: k10plu,k10min
 !   real :: k11plu,k11min
 !   real :: k9in,k10in,k2in

 !   real :: kplamin,kplaplu,kadpin
 !   real :: kplaact,ke2actAe2,kact,kactadp

 !   real :: kz1e2plu,kz1e2min,kz1e2cat


    ! FXI NEW

 !     double precision kz11on,kz11off,ke11on,ke11off
 !     double precision kz9me11mp,kz9me11mmi,kz9me11mcat
 !     double precision kz11me11mp,kz11me11mmi,kz11me11mcat
 !     double precision kz11me2mp,kz11me2mmi,kz11me2mcat
 !     double precision kz11e11p,kz11e11mi,kz11e11cat
 !     double precision kz9e11p,kz9e11mi,kz9e11cat
 !     double precision kz11e2p,kz11e2mi,kz11e2cat
 !     double precision k11in
 !     double precision z11_init,e11_init,z11m_init,e11m_init
 !     double precision z11mBe2m_init,z9mBe11m_init
 !     double precision z11Be2_init,z9Be11_init
 !     double precision z11Be11_init,z11mBe11m_init
 !        
 !     double precision kz9e11hp,kz9e11hmi,kz9e11hcat
 !     double precision kz9me11hmp,kz9me11hmmi,kz9me11hmcat
 !     double precision kz9me11msp,kz9me11msmi,kz9me11mscat

 !     double precision kz11e11hp,kz11e11hmi,kz11e11hcat
 !     double precision kz11me11hmp,kz11me11hmmi,kz11me11hmcat
 !     double precision kz11me11msp,kz11me11msmi,kz11me11mscat

 !     double precision ke11he11p,ke11he11mi,ke11he11cat
 !     double precision ke11he11hp,ke11he11hmi,ke11he11hcat
 !     double precision ke11he2p,ke11he2mi,ke11he2cat

 !     double precision ke11hmse11hmp,ke11hmse11hmmi,ke11hmse11hmcat
 !     double precision ke11hmse11msp,ke11hmse11msmi,ke11hmse11mscat
 !     double precision ke11hmse2mp,ke11hmse2mmi,ke11hmse2mcat

 !     double precision ke11ons,ke11offs
 !     double precision ke11hons,ke11hoffs
 !     double precision ke11hon,ke11hoff
    


    ! NORMAL PLASMA CONCENTRATIONS

 !   real :: va,v2,v5,v7,v8,v9,v10
 !   real :: ve7,vf,v11
 !   real :: v11a,v1

    ! SURFACE DENSITIES
    
 !   real :: s7d,s5,s10,s9,s91,s8,s2
 !   real :: s11,s111

  !  real :: s2s,s2g !thrombin on platelet



    ! SURFACE BINDING SITE NUMBERS

  !  real :: N2b,N5b,N8b,N9b,N9starb,N10b
  !  real :: N2se,N5se,N8se,N9se,N9starse,N10se
  !  real :: N5,N11b,N11se,N111b,N111se

  !  real :: N2starb,N2starse ! new thrombin binding sites
  !  real :: N2gp1bb,N2gp1bse ! new thrombin binding sites




    ! PARAMETERS FOR PLATELET AVERAGING
    
  !  real    :: hbreal             ! physical length of one grid cell
  !  integer :: avgipm             ! number of indices out plus/minus to average 
                                  !  1 for 9 cells, 2 for 25 etc
    ! PARAMETER FOR G(ETA)
  !  real    :: gammaeta,betageta !:)
  !  real    :: etat,etastar
  !  real    :: Deta

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !            FOR ADP RELEASE
    !     length of Rtemp  -    5/rlen = 0.25
    ! this means each time to release is .25 seconds
  !  integer,parameter :: rlen = 20
    real    :: tsigma,tcount
  !  real    :: rtemp(1:rlen)
  !  real    :: sigmatemp(1:rlen)
    ! sigmar holds all future information
 !   real    :: sigmar(1:nxb,1:nyb,1:rlen)
 !   ! sigmaz is the CURRENT release value
 !   real    :: sigmaz(1:nxb,1:nyb) 
 !   real    :: ir !integral of rtemp
 !   real    :: Pbold(1:nxb,1:nyb)
 !   real    :: Pbnew(1:nxb,1:nyb)
 !   real    :: ADPrel
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



    real    :: omegaimp

    ! FOR HINDERED TRANSPORT
    real    :: A_h(-1:nxb+2,-1:nyb+2) ! coefficient for hindered advection
    real    :: D_h(-1:nxb+2,-1:nyb+2) ! coefficient for hindered diffusion
    real    :: rh         ! ratio for hinderance
    integer :: adpf                   ! adp flag -- if 1, hinder adp, if 0 don't
    real    :: kh                     ! amount to scale phi_b to get to theta_p
    real    :: gamma_cz
CONTAINS

                                   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! SUBROUTINE: SETPARAMS -- set parameter values and param values that are 
!                          combinations of other parameters
!
subroutine setparams
implicit none
integer :: i,j   ! for looping

  ! FLUID PARAMS
  !
  ! Scales of variables in cgs units
  !

  xchar = pltdiam
  tchar = tmax            ! try one millisecond for advection time scale
  dt    = 1e-8 !1e-6      !dt =  make sure this is less than hb/umax/uchar!!! 
  uchar = xchar/tchar     ! set U=L/T changes are now in the fluidsolver
  !pchar = mu*uchar/xchar  ! High Re scaling
  pchar = rho*uchar**2    ! Low Re scaling
  fchar = pchar/xchar
  write(*,*) 'tmax = ', tmax, 'hb =', hb,'umax =', umax,'uchar = ', uchar
!stop
  ! wall shear rate = 2*umax/radius = 1/tchar
  gamma_cz =10.0
  umax = gamma_cz*dheight/4d0!gamma_cz*ymax/4d0 !1d-1*5d1*dheight/(4d0)!1d-1*0.5*dheight/(2.0*tchar)
  write(*,*) dheight, 'ymax=', ymax
!stop
  !write(*,*) tchar
  !write(*,*) 0.5*dheight/(2.0*tchar)

  write(*,*)'umax,uchar',umax,uchar
  
!stop
  !
  ! Dimensionless parameters
  !
  !   Re  -- Reynolds number
  !

  Re  = (xchar*uchar*rho)/mu   
  print *, Re
write(*,*) 'Re=', Re

  ! INJURY SITE PARAMS
  !
  ! rzL is location of left side of reaction zone
  ! rzR is location of left side of reaction zonde
  ! rzLi,rzRi are the indices at the edges
  !
write(*,*) 'xmax=', xmax
  ilengthnd = 0.25*xmax
  rzL       = 0.5*xmax - 0.5*ilengthnd
  rzR       = 0.5*xmax + 0.5*ilengthnd
  rzLi      = 40 !int(nxb*rzL/xmax)
  rzRi      = 88 !int(nxb*rzR/xmax)  
  write(*,*)'rzLi,rzRi',rzLi,rzRi
  write(*,*)'rzLi,rzRi,hb',rzL,rzR,hb

  
  ! BIOCHEMICAL PARAMS
  !
  !    Platelets
  
 ! P0        = 4.151e-13 ! M 
 ! DCP       = 5.0e-7  ! cm^2/s
 ! Pmaxb     = 1.1071e-10 ! psub in M
 ! Pmaxse    = 1.1071e-10 ! psub in M
 ! kcoh      = 1.0e4   ! cohesion rate 1/s
 ! krel      = 0.0     ! release rate 1/s
 ! kadh_max  = 2.0e10  ! max for kadh
 ! kadh=0.0
 ! kadh(rzLi:rzRi,1:2)=kadh_max
 ! DiffCoeffP=DCP


  !
  ! KINETIC RATES
  !  
!    k7on  = 5.0e7
!    k7off = 5.0e-3
!    k2cat = 5.0
!    k2min = 1.0
!    k2plu = 5.0e6
!    k3cat = 6.1e-2
!    k3min = 1.0
!    k3plu = 3.92e5
!    k8cat = 1.15
!    k8min = 1.0
!    k8plu = 8.95e6
!    k8starm = k8min
!    k8starp = k9plu
!    k9cat = 1.15
!    k9min = 1.0
!    k9plu = 8.95e6

!    k1cat = 5.0
!    k1min = 1.0
!    k1plu = 5.0e6
!    k18cat= 6.1e-2
!    k18min= 1.0
!    k18plu= 3.92e5
!    k12cat= 0.23
!    k12min= 1.0
!    k12plu= 1.73e7
!    k14cat= 0.9
!    k14min= 1.0
!    k14plu= 2.64e7

!    k9off = 2.5e-2
!    k9on  = 1.0e7
!    k10off= 2.5e-2
!    k10on = 1.0e7
!    k10zon = k10on
!    k10zoff = 0.5
!    k5off = 0.17
!    k5on  = 5.7e7
!    k8off = 0.17
!    k8on  = 5.0e7
!    k2off = 5.9
!    k2on  = 1.0e7

    ! new thrombin binding sites
!    ke2on = 1.0e7
!    ke2off = 0.2

!    ke2gp1bon = 1.0e8
!    ke2gp1bon = 0.0
!    ke2gp1boff = 1.0

!    k5cat = 4.6e-2
!    k5min = 1.0
!    k5plu = 1.0e8
!    k13cat= 0.23
!    k13min= 1.0
!    k13plu= 1.73e7
!    k6cat = 2.3e-2
!    k6min = 1.0
!    k6plu = 5.1e7
!    k15cat= 0.9
!    k15min= 1.0
!    k15plu= 2.64e7
!    ktenmin = 0.01
!    ktenplu = 1.0e8
!    kpromin = 0.01
!    kproplu = 1.0e8
!    k4cat = 20.0
!    k4min = 1.0
!    k4plu = 1.31e8
!    k7cat = 30.0
!    k7min = 1.0
!    k7plu = 1.03e8

!    k17cat = 0.5
!    k17min = 1.0
!    k17plu = 1.2e8
!    k16cat = 0.5
!    k16min = 1.0
!    k16plu = 1.2e8
!    k10min = 3.3e-4
!    k10plu = 1.6e7
!    k11min = 1.1e-3
!    k11plu = 1.0e7
!    k9in   = 0.00118
!    k10in  = 0.0067
!    k2in   = 0.017

!    kadpin = 0.0 ! WHAT SHOULD THIS BE?

!    kplamin = 0.0
!    kplaplu = 2.0e10
!    kplaact = 3.0e8
!    ke2actAe2 = 0.5
!    kact      = 0.5
!    kactadp   = 0.34


    ! FXI params

    ! fxi platelet binding constants
    ! on/off standard
!    kz11on=1.0d7
!    kz11off=0.1d0
!    ke11on=1.0d7
!    ke11off=0.017d0

    ! on/off high
    !       kz11on=1.0d8
    !      kz11off=1.0d0
    !       ke11on=1.0d8
    !       ke11off=0.17d0
    ! on/off low
    !       kz11on=1.0d6
    !       kz11off=0.01d0
    !       ke11on=1.0d6
    !       ke11off=0.0017d0
        
    ! fxi kinetic parameters

    ! activation of z9m by e11m

    ! standard      
   
 !   kz9me11mp=0.6d7
 !   kz9me11mmi=1.d0
 !   kz9me11mcat=0.21d0

    !A1 set
    !        kz9me11mp=0.155d7
    !        kz9me11mmi=0.1d0
    !        kz9me11mcat=0.21d0
    ! C1 set
    !        kz9me11mp=1.55d7
    !c        kz9me11mmi=1.0d0
    !c        kz9me11mcat=2.1d0
    !c E1 set
    !c        kz9me11mp=2.6d6
    !c        kz9me11mmi=0.1d0
    !c        kz9me11mcat=0.42d0
    !c
    
    ! activation of z11m by e11m (autoactivation)
    !       
!    kz11me11mp=0.d0
!    kz11me11mmi=0.d0
!    kz11me11mcat=0.d0
    !
    ! activation of z11m by e2m 
    !         
    ! standard      
        
!    kz11me2mp=2.0d7
!    kz11me2mmi=1.0d0
!    kz11me2mcat=1.3d-4

    ! B1 set
    !        kz11me2mp=2.0d6
    !        kz11me2mmi=0.1d0
    !        kz11me2mcat=1.3d-4
    ! D1 set
    !        kz11me2mp=2.0d5
    !        kz11me2mmi=0.01d0
    !        kz11me2mcat=1.3d-3
    ! F1 set
    !        kz11me2mp=2.05d5
    !        kz11me2mmi=0.01d0
    !        kz11me2mcat=2.6d-3
    

    ! activation of z11 by e11 (autoactivation)
    !               
!    kz11e11p=0.d0
!    kz11e11mi=0.d0
!    kz11e11cat=0.d0
    
    ! activation of z9 by e11
       
    ! standard
!    kz9e11p=0.6d7
!    kz9e11mi=1.0d0
!    kz9e11cat=0.21d0

    ! A1 set
    !        kz9e11p=0.155d7
    !        kz9e11mi=0.1d0
    !c        kz9e11cat=0.21d0
    !c C1 set
    !c        kz9e11p=1.55d7
    !c        kz9e11mi=1.0d0
    !c        kz9e11cat=2.1d0
    !c E1 set
    !c        kz9e11p=2.6d6
    !c        kz9e11mi=0.1d0
    !c        kz9e11cat=0.42d0
    !c
    
    ! activation of z11 by e2 
                       
    ! standard
!    kz11e2p=2.0d7
!    kz11e2mi=1.0d0
!    kz11e2cat=1.3d-4

    !c B1 set
    !c        kz11e2p=2.0d6
    !c        kz11e2mi=0.1d0
    !c        kz11e2cat=1.3d-4
    !c D1 set
    !c        kz11e2p=2.0d5
    !c        kz11e2mi=0.01d0
    !c        kz11e2cat=1.3d-3
    !c F1 set
    !c        kz11e2p=2.05d5
    !c        kz11e2mi=0.01d0
    !c        kz11e2cat=2.6d-4
    !c
    
    ! first order inhibition of fxi 

!    k11in = 0.01
    
    ! noxiinh
    !        k11in=0.0d0


!    ke11ons=         ke11on
!    ke11offs=        ke11off
!    ke11hons=       ke11on
!    ke11hoffs=       ke11off
!    ke11hon=         kz11on
!    ke11hoff=        kz11off
    
!    kz9e11hp=         kz9e11p
!    kz9e11hmi=       kz9e11mi
!    kz9e11hcat=      kz9e11cat
    
!    kz9me11hmp=       kz9me11mp
!    kz9me11hmmi=      kz9me11mmi
!    kz9me11hmcat=     kz9me11mcat
    
!    kz9me11msp=       kz9me11mp
!    kz9me11msmi=      kz9me11mmi
!    kz9me11mscat=     kz9me11mcat
    
!    kz11e11hp=       0.d0
!    kz11e11hmi=      0.d0
!    kz11e11hcat=     0.d0
    
!    kz11me11hmp=     0.d0
!    kz11me11hmmi=    0.d0
!    kz11me11hmcat=   0.d0
    
!    kz11me11msp=     0.d0
!    kz11me11msmi=    0.d0
!    kz11me11mscat=   0.d0
    
!    ke11he11p=       0.d0
!    ke11he11mi=      0.d0
!    ke11he11cat=     0.d0
    
!    ke11he11hp=      0.d0
!    ke11he11hmi=     0.d0
!    ke11he11hcat=    0.d0
    
!    ke11he2p=         kz11e2p
!    ke11he2mi=       kz11e2mi
!    ke11he2cat=      kz11e2cat
    
!    ke11hmse11hmp=   0.d0
!    ke11hmse11hmmi=  0.d0
!    ke11hmse11hmcat= 0.d0
    
!    ke11hmse11msp=   0.d0
!    ke11hmse11msmi=  0.d0
!    ke11hmse11mscat= 0.d0
    
 !   ke11hmse2mp=     kz11me2mp
!    ke11hmse2mmi=    kz11me2mmi
!    ke11hmse2mcat=   kz11me2mcat


    ! FIBRINOGEN
!    kz1e2plu = 1.18e7
!    kz1e2min = 1.0
!    kz1e2cat = 84.0


!    v11a=0.0


    ! PLATELET BINDING SITES

    ! set to original value with no new thrombin bs
!    N2b = 1000.0
!    N2starb = 1000.0
!    N2gp1bb = 4000.0
    
!    N5b = 3000.0
!    N8b = 450.0
!    N9b = 250.0
!    N9starb = 250.0
!    N10b = 2700.0


!    N2se = 1000.0
!    N2starse = 1000.0
!    N2gp1bse = 4000.0

!    N5se = 3000.0
!    N8se = 450.0
!    N9se = 250.0
!    N9starse = 250.0
!    N10se = 2700.0

!    N5 = 3000.0

!    N11b   =  1500.0 
!    N11se  =  1500.0
!    N111b  =   250.0
!    N111se =   250.0

! NORMAL VOLUME CONCENTRATIONS - ALL IN M

!    va  = 5.0e-6  ! 5.0  microM
!    v2  = 1.4e-6  ! 1.4  microM
!    v5  = 0.01e-6 ! 0.01 microM 
!    v7  = 0.01e-6 ! 0.01 microM

!    v8  = 1.0e-9  ! 1.0  nM
    !Hem A
    !v8=0.01*1.0e-9 

!    v9  = 0.09e-6 ! 0.09 microM
    !Hem B
    !v9=0.01*0.09e-6 

!    v10 = 0.17e-6 ! 0.17 microM
!    ve7 = 0.1e-9  ! 0.1  nM
!    vf  = 2.5e-9  ! 2.5  nM

!    v11 = 0.03e-6 ! 0.03 microM

    !FIBRINOGEN
!    v1 = 8.8e-6   ! 8.8 microM


   ! NORMAL SURFACE CONCENTRATIONS 
   ! SE IN FMOLES/CM^2
   ! PB IN moles of binding sites/liter

    !s7d = 15.0e-12 ! see notes
    
    ! for testing
    !s7d = 1.0e-12 ! see notes
    !s7d = 30.0e-12 ! see notes


!    s5  = N5b*Pmaxb + N5se*Pmaxse
!    s10 = N10b*Pmaxb + N10se*Pmaxse
!    s9  = N9b*Pmaxb + N9se*Pmaxse
!    s91 = N9starb*Pmaxb + N9starse*Pmaxse
!    s8  = N8b*Pmaxb + N8se*Pmaxse
!    s2  = N2b*Pmaxb + N2se*Pmaxse
    ! new thrombin binding sites
!    s2s  = N2starb*Pmaxb + N2starse*Pmaxse
!    s2g  = N2gp1bb*Pmaxb + N2gp1bse*Pmaxse
!    s11 = N11b*Pmaxb + N11se*Pmaxse
!    s111= N111b*Pmaxb + N111se*Pmaxse



! **********************************************


  ! 
  ! NON-DIMENSIONAL PARAMS
  !
 
  ! pl
 ! ndpl(1) = P0*6.022e17
 ! ndpl(2) = P0*6.022e17
 ! ndpl(3) = Pmaxse*6.022e17
 ! ndpl(4) = Pmaxb*6.022e17
 ! ndpl(5) = 1.0
  ! se
 ! ndse(1) = s7d*1e12 ! put in fmol/cm^2
 ! ndse(2) = s7d*1e12 ! put in fmol/cm^2
 ! ndse(3) = s7d*1e12 ! put in fmol/cm^2
 ! ndse(4) = s7d*1e12 ! put in fmol/cm^2
 ! ndse(5) = s7d*1e12 ! put in fmol/cm^2
 ! ndse(6) = s7d*1e12 ! put in fmol/cm^2
 ! ndse(7) = s7d*1e12 ! put in fmol/cm^2
 ! ndse(8) = s7d*1e12 ! put in fmol/cm^2
 ! ndse(9) = 1.0      ! fraction of the se covered
  ! fp
 ! ndfp(1) = v2*1.0e9 ! put in nM
 ! ndfp(2) = v2*1.0e9 ! put in nM
 ! ndfp(3) = v5*1.0e9 ! put in nM
 ! ndfp(4) = v5*1.0e9 ! put in nM
 ! ndfp(5) = v7*1.0e9 ! put in nM
 ! ndfp(6) = v7*1.0e9 ! put in nM
 ! ndfp(7) = v8*1.0e9 ! put in nM
 ! ndfp(8) = v8*1.0e9 ! put in nM
 ! ndfp(9) = v9*1.0e9 ! put in nM
 ! ndfp(10)= v9*1.0e9 ! put in nM
 ! ndfp(11)= v10*1.0e9! put in nM
 ! ndfp(12)= v10*1.0e9! put in nM
 ! ndfp(13)= v5*1.0e9 ! put in nM
 ! ndfp(14)= v7*1.0e9 ! put in nM
 ! ndfp(15)= v7*1.0e9 ! put in nM
 ! ndfp(16)= v8*1.0e9 ! put in nM
 ! ndfp(17)= v2*1.0e9 ! put in nM
 ! ndfp(18)= vf*1.0e9 ! put in nM
 ! ndfp(19)= vf*1.0e9 ! put in nM
 ! ndfp(20)= va*1.0e6 ! put in microM
 ! ndfp(21)= v11*1.0e9 ! put in nM
 ! ndfp(22)= v11*1.0e9 ! put in nM
 ! ndfp(23)= v11*1.0e9 ! put in nM
 ! ndfp(23)= v11*1.0e9 ! put in nM
 ! ndfp(24)= v9*1.0e9 ! put in nM
 ! ndfp(25)= v9*1.0e9 ! put in nM
 ! ndfp(26)= v11*1.0e9 ! put in nM
 ! ndfp(27)= v11*1.0e9 ! put in nM
 ! ndfp(28)= v11*1.0e9 ! put in nM
 ! ndfp(29)= v11*1.0e9 ! put in nM
 ! ndfp(30)= v11*1.0e9 ! put in nM
 ! ndfp(31)= v11*1.0e9 ! put in nM
 ! ndfp(32)= v1*1.0e9 ! put in nM
 ! ndfp(33)= v1*1.0e9 ! put in nM


  ! pb
 ! ndpb(1) = s2*1.0e9 ! put in nM
!  ndpb(2) = s2s*1.0e9 ! put in nM nd by e2m
!  ndpb(3) = s5*1.0e9 ! put in nM
!  ndpb(4) = s5*1.0e9 ! put in nM
!  ndpb(5) = s8*1.0e9 ! put in nM
!  ndpb(6) = s8*1.0e9 ! put in nM
!  ndpb(7) = s9*1.0e9 ! put in nM
!  ndpb(8) = s9*1.0e9 ! put in nM
!  ndpb(9) = s10*1.0e9! put in nM
!  ndpb(10)= s10*1.0e9! put in nM
!  ndpb(11)= s8*1.0e9 ! put in nM
!  ndpb(12)= s5*1.0e9 ! put in nM
!  ndpb(13)= s5*1.0e9 ! put in nM
!  ndpb(14)= s2s*1.0e9 ! put in nM nd by e2m
!  ndpb(15)= s10*1.0e9! put in nM
!  ndpb(16)= s8*1.0e9 ! put in nM
!  ndpb(17)= s8*1.0e9 ! put in nM
!  ndpb(18)= s8*1.0e9 ! put in nM
!  ndpb(19)= s5*1.0e9 ! put in nM
!  ndpb(20)= s8*1.0e9 ! put in nM
!  ndpb(21)= s91*1.0e9! put in nM
!  ndpb(22)= s8*1.0e9 ! put in nM
!  ndpb(23)= s8*1.0e9 ! put in nM
!  ndpb(24)=P0*6.022e17 ! tracking # plts activated
!  ndpb(25)=P0*6.022e17
!  ndpb(26)=P0*6.022e17
!  ndpb(27)= s11*1.0e9 ! put in nM
!  ndpb(28)= s11*1.0e9 ! put in nM
!  ndpb(29)= s111*1.0e9 ! put in nM
!  ndpb(30)= s111*1.0e9 ! put in nM

!  ndpb(31)= s9*1.0e9 ! put in nM
!  ndpb(32)= s9*1.0e9 ! put in nM
!  ndpb(33)= s11*1.0e9 ! put in nM
!  ndpb(34)= s11*1.0e9 ! put in nM
!  ndpb(35)= s11*1.0e9 ! put in nM
!  ndpb(36)= s111*1.0e9 ! put in nM
!  ndpb(37)= s111*1.0e9 ! put in nM
!  ndpb(38)= s111*1.0e9 ! put in nM
!  ndpb(39)= s2g*1.0e9 ! put in nM




  !NDadvect  = uchar*tchar/xchar 
!  NDdiffuse = tchar*DCP/(xchar**2)

!  write(*,*)'coeffs for adv,diff: ',NDadvect,NDdiffuse
write(*,*) 'dt=', dt, 'hb/(umax/uchar)=', hb/(umax/uchar)  
if (dt.gt. hb/(umax/uchar))then
    write(*,*)'CFL condtion violated=', hb/(umax/uchar)
    write(*,*)'stopping!', 'hb=', hb, 'umax=', umax, 'uchar=', uchar
    stop
  end if

! AVERAGING OF PLATELETS

!  hbreal = hb*xchar
!  avgipm = 1!floor(pltdiam/hbreal)

!  write(*,*)'avgipm,hbreal',avgipm,hbreal

! ADP RELEASE

!  rtemp(1)=0.0;
!  ir=0
!  do i=2,rlen
!     tsigma=(i-1)*5.0/(rlen-1)
!     rtemp(i)=exp(-(tsigma-3.0)**2)
!     ir=ir+rtemp(i)*5.0/(rlen-1)
!  end do
!  rtemp=rtemp/ir
  ! rtemp was tested and is correct 11/7/08

!  Pbold = 0.0
!  Pbnew = 0.0

!  ADPrel = 2.0e-17

  ! G(ETA)
  ! gamma = sqrt(10*D/decaydistance^2)
!  Deta=10.0*DCP
!  gammaeta = Deta/((3.0e-4)**2)
!  write(*,*)'Deta',Deta
!  etat = 0.1
!  etastar = 0.5-etat

!  betageta = 1.0/((1.0-etat)**3/(etastar**3+(1.0-etat)**3))
!  write(*,*)'beta',betageta


  ! for implicit RK
!  omegaimp = 1.0

  ! initialize for hindered transport
!  A_h(:,:)=1.0
!  D_h(:,:)=1.0


end subroutine setparams
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! SUBROUTINE BIOCHEMINIT -- sets up initial values for biochemical variables
!
!subroutine biocheminit
!  integer :: i  ! for looping
!  pbdryL=0.0  
!  plinit(:,:)=0.0  
!  call pltbdry(pbdryL)

 ! to simulate platelets chasing saline comment these lines:

!  do i=1,nxb
!   plinit(i,1:nyb)=pbdryL
!  end do

!  seinit=0.0
!  seinit(rzLi:rzRi)=1.0

!  fpinit=0.0
!  pbinit=0.0
  !fpinit(rzLi:rzRi,1)=1.0-P0/Pmaxse
  !fpinit(rzLi:rzRi,2)=1.0-P0/Pmaxse

  ! INITIALIZE ALL SPECIES - NONDIMENSIONAL
  
  !PLATELETS (5)
  ! Pm,u
  !pl(1)=platelet('PLPmu',0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,1)
!  pl(1)=platelet('PLPmu',plinit,1.0,0.0,0.0,0.0,0.0,0.0,0.5,1)
  ! Pm,a
!  pl(2)=platelet('PLPma',0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.5,2) 
  ! Pse,a
!  pl(3)=platelet('PLPsea',fpinit,0.0,0.0,0.0,0.0,0.0,0.0,1.0,3)
  ! Pb,a
!  pl(4)=platelet('PLPba',0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,4)
  ! eta NOTE EQN NUM
!  pl(5)=platelet('PLeta',0.0,0.0,0.0,0.0,0.0,0.0,0.0,10.0,60)


  ! SUBENDOTHELIUM-BOUND (9)
  ! z7se
!  se(1)=sebound('SEz7se',0.0,0.0,0.0,0.0,0.0,0.0,6)
!  ! e7se
!  se(2)=sebound('SEe7se',0.0,0.0,0.0,0.0,0.0,0.0,7)
  ! z7se:e2
!  se(3)=sebound('SEz7se_e2',0.0,0.0,0.0,0.0,0.0,0.0,8)
  ! z7se:e10
!  se(4)=sebound('SEz7se_e10',0.0,0.0,0.0,0.0,0.0,0.0,9)
  ! z9:e7se
!  se(5)=sebound('SEz9_e7se',0.0,0.0,0.0,0.0,0.0,0.0,10)  
  ! z10:e7se
!  se(6)=sebound('SEz10_e7se',0.0,0.0,0.0,0.0,0.0,0.0,11)  
  ! tfpi:e10:e7se
!  se(7)=sebound('SEtfpi_e10_e7se',0.0,0.0,0.0,0.0,0.0,0.0,12)
  ! e10:e7se
!  se(8)=sebound('SEe10_e7se',0.0,0.0,0.0,0.0,0.0,0.0,13)   
  ! tf
  !se(9)=sebound('SEtf',0.0,0.0,0.0,0.0,0.0,0.0,14) 
!  se(9)=sebound('SEtf',seinit,0.0,0.0,0.0,0.0,0.0,14) 


  ! FLUID PHASE (31)

  ! z2
  !fp(1)=fluidphase('FPz2',0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,15)
 ! fp(1)=fluidphase('FPz2',1.0,1.0,0.0,0.0,0.0,0.0,0.0,1.0,15)

  ! e2
!  fp(2)=fluidphase('FPe2',0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,16)

  ! z5
  !fp(3)=fluidphase('FPz5',0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,17) 
!  fp(3)=fluidphase('FPz5',1.0,1.0,0.0,0.0,0.0,0.0,0.0,1.0,17) 

  ! e5
!  fp(4)=fluidphase('FPe5',0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,18)

  ! z7
!  fp(5)=fluidphase('FPz7',1.0,1.0,0.0,0.0,0.0,0.0,0.0,1.0,19)  
 
  ! e7
  !fp(6)=fluidphase('FPe7',0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,20)
!  fp(6)=fluidphase('FPe7',0.01,0.01,0.0,0.0,0.0,0.0,0.0,1.0,20)

  ! z8
  !fp(7)=fluidphase('FPz8',0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,21) 
!  fp(7)=fluidphase('FPz8',1.0,1.0,0.0,0.0,0.0,0.0,0.0,1.0,21)

  ! e8
!  fp(8)=fluidphase('FPe8',0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,22)

  ! z9
  !fp(9)=fluidphase('FPz9',0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,23) 
!  fp(9)=fluidphase('FPz9',1.0,1.0,0.0,0.0,0.0,0.0,0.0,1.0,23) 

  ! e9
!  fp(10)=fluidphase('FPe9',0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,24)

  ! z10
  !fp(11)=fluidphase('FPz10',0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,25)   
!  fp(11)=fluidphase('FPz10',1.0,1.0,0.0,0.0,0.0,0.0,0.0,1.0,25)

  ! e10
!  fp(12)=fluidphase('FPe10',0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,26)  

  ! z5:e2
!  fp(13)=fluidphase('FPz5_e2',0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,27)

  ! z7:e2
!  fp(14)=fluidphase('FPz7_e2',0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,28) 
  
  ! z7:e10
!  fp(15)=fluidphase('FPz5_e10',0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,29)

  ! z8:e2
!  fp(16)=fluidphase('FPz8_e2',0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,30)

  ! apc
!  fp(17)=fluidphase('FPapc',0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,31)

  ! tfpi
  !fp(18)=fluidphase('FPtfpi',0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,32) 
!  fp(18)=fluidphase('FPtfpi',1.0,1.0,0.0,0.0,0.0,0.0,0.0,1.0,32)

  ! tfpi:e10
!  fp(19)=fluidphase('FPtfpi_e10',0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,33) 
 
  ! adp
!  fp(20)=fluidphase('FPadp',0.0,0.0,0.0,0.0,0.0,0.0,0.0,9.6,34) 
 

  ! FXI STARTING NEW EQN NUMBERS HERE FROM 62

  ! z11
!  fp(21)=fluidphase('FPz11',0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,62)  

  ! e11h
!  fp(22)=fluidphase('FPe11h',0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,63)

  ! e11
!  fp(23)=fluidphase('FPe11',v11a,v11a,0.0,0.0,0.0,0.0,0.0,1.0,64)

  ! z9:ellh
!  fp(24)=fluidphase('FPz9_e11h',0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,65)

  ! z9:ell
!  fp(25)=fluidphase('FPz9_e11',0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,66)

  ! z11:ellh
!  fp(26)=fluidphase('FPz11_e11h',0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,67)

  ! z11:ell
!  fp(27)=fluidphase('FPz11_e11',0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,68)

  ! z11:e2
!  fp(28)=fluidphase('FPz11_e2',0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,69)

  ! e11h:ellh
!  fp(29)=fluidphase('FPe11h_e11h',0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,70)

  ! e11h:ell
!  fp(30)=fluidphase('FPe11h_e11',0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,71)

  ! e11h:e2
!  fp(31)=fluidphase('FPe11h_e2',0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,72)

  ! z1
  !fp(32)=fluidphase('FPz1',1.0,1.0,0.0,0.0,0.0,0.0,0.0,1.0,85)
!  fp(32)=fluidphase('FPz1',0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,85)

  ! z1:e2
!  fp(33)=fluidphase('FPz1_e2',0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,86)

  !********************************************************

  ! PLATELET-BOUND (38)

  ! z2ba
!  pb(1)=pltbound('PBz2ba',0.0,0.0,0.0,0.0,0.0,0.0,35)

  ! e2ba
!  pb(2)=pltbound('PBe2ba',0.0,0.0,0.0,0.0,0.0,0.0,36)  

  ! z5ba
!  pb(3)=pltbound('PBz5ba',0.0,0.0,0.0,0.0,0.0,0.0,37)

  ! e5ba
!  pb(4)=pltbound('PBe5ba',0.0,0.0,0.0,0.0,0.0,0.0,38)

  ! z8ba
!  pb(5)=pltbound('PBz8ba',0.0,0.0,0.0,0.0,0.0,0.0,39)

  ! e8ba
!  pb(6)=pltbound('PBe8ba',0.0,0.0,0.0,0.0,0.0,0.0,40)

  ! z9ba
!  pb(7)=pltbound('PBz9ba',0.0,0.0,0.0,0.0,0.0,0.0,41)

  ! e9ba
!  pb(8)=pltbound('PBe9ba',0.0,0.0,0.0,0.0,0.0,0.0,42) 
 
  ! z10ba
!  pb(9)=pltbound('PBz10ba',0.0,0.0,0.0,0.0,0.0,0.0,43)

  ! e10ba
!  pb(10)=pltbound('PBe10ba',0.0,0.0,0.0,0.0,0.0,0.0,44)

  ! ten		      
!  pb(11)=pltbound('PBten',0.0,0.0,0.0,0.0,0.0,0.0,45)

  ! pro	      
!  pb(12)=pltbound('PBpro',0.0,0.0,0.0,0.0,0.0,0.0,46)  

  ! z2ba:pro	      
!  pb(13)=pltbound('PBz2ba_pro',0.0,0.0,0.0,0.0,0.0,0.0,47)

  ! z5ba:e2ba	      
!  pb(14)=pltbound('PBz5ba_e2ba',0.0,0.0,0.0,0.0,0.0,0.0,48) 

  ! z5ba:e10ba	      
!  pb(15)=pltbound('PBz5ba_e10ba',0.0,0.0,0.0,0.0,0.0,0.0,49) 

  ! z8ba:e2ba	      
!  pb(16)=pltbound('PBz8ba_e2ba',0.0,0.0,0.0,0.0,0.0,0.0,50)  
   
  ! z8ba:e10ba	      
!  pb(17)=pltbound('PBz8ba_e10ba',0.0,0.0,0.0,0.0,0.0,0.0,51)    

  ! z10ba:ten	      
!  pb(18)=pltbound('PBz10ba_ten',0.0,0.0,0.0,0.0,0.0,0.0,52)

  ! apc:e5ba	      
!  pb(19)=pltbound('PBapc_e5ba',0.0,0.0,0.0,0.0,0.0,0.0,53) 

  ! apc:e8ba	      
!  pb(20)=pltbound('PBapc_e8ba',0.0,0.0,0.0,0.0,0.0,0.0,54)

  ! e9starba	      
!  pb(21)=pltbound('PBe9starba',0.0,0.0,0.0,0.0,0.0,0.0,55)

  ! tenstar
!  pb(22)=pltbound('PBtenstar',0.0,0.0,0.0,0.0,0.0,0.0,56) 

  ! z10ba:tenstar	      
!  pb(23)=pltbound('PBz10ba_tenstar',0.0,0.0,0.0,0.0,0.0,0.0,57)   
  
  ! Plts act by adp	      
!  pb(24)=pltbound('PBPltsADP',0.0,0.0,0.0,0.0,0.0,0.0,58)     

  ! Plts act by thrombin	      
!  pb(25)=pltbound('PBPltsE2',0.0,0.0,0.0,0.0,0.0,0.0,59)    
 
  ! Plts binding - NOTE EQN # 61 since eta is # 60
!  pb(26)=pltbound('PBPltsBind',0.0,0.0,0.0,0.0,0.0,0.0,61)   


  ! FXI START FROM 73

  ! z11m
!  pb(27)=pltbound('PBz11m',0.0,0.0,0.0,0.0,0.0,0.0,73) 

  ! e11hm
!  pb(28)=pltbound('PBe11hm',0.0,0.0,0.0,0.0,0.0,0.0,74) 

  ! e11hmstar
!  pb(29)=pltbound('PBe11hms',0.0,0.0,0.0,0.0,0.0,0.0,75) 

  ! e11mstar
!  pb(30)=pltbound('PBe11ms',0.0,0.0,0.0,0.0,0.0,0.0,76)

  ! z9m:ellhm
!  pb(31)=pltbound('PBz9m_e11hm',0.0,0.0,0.0,0.0,0.0,0.0,77)

  ! z9m:ellms
!  pb(32)=pltbound('PBz9m_e11ms',0.0,0.0,0.0,0.0,0.0,0.0,78)

  ! z11m:ellhm
!  pb(33)=pltbound('PBz11m_e11hm',0.0,0.0,0.0,0.0,0.0,0.0,79)

  ! z11m:ellms
!  pb(34)=pltbound('PBz11m_e11ms',0.0,0.0,0.0,0.0,0.0,0.0,80)

  ! z11m:e2m
!  pb(35)=pltbound('PBz11m_e2m',0.0,0.0,0.0,0.0,0.0,0.0,81)

  ! e11hms:e11hm
!  pb(36)=pltbound('PBe11hms_e11hm',0.0,0.0,0.0,0.0,0.0,0.0,82)

  ! e11hms:e11ms
!  pb(37)=pltbound('PBe11hms_e11ms',0.0,0.0,0.0,0.0,0.0,0.0,83)

  ! e11hms:e2m
!  pb(38)=pltbound('PBe11hms_e2m',0.0,0.0,0.0,0.0,0.0,0.0,84)

  ! e2ba_gp1b
!  pb(39)=pltbound('PBe2ba_gp1b',0.0,0.0,0.0,0.0,0.0,0.0,87)

!end subroutine biocheminit
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! SUBROUTINE: OUTPUTSTEPINFO -- output information about time steps and 
!                               data file numbers
!
subroutine outputstepinfo( un )
  integer :: un
  integer :: outshift
  integer :: restartflag

  character(len=20) :: date
  character(len=20) :: time

  ! for restarts, the initial data is NOT output
  !
  if( restartflag.eq.0 )then
     outshift = 0
  else
     outshift = 1
  endif
 
  call date_and_time(date,time)
  write(un,'(a6,a2,a1,a2,a1,a4,a4,a2,a1,a2)')             &
       'Run:  ',date(5:6),'/',date(7:8),'/',date(1:4),' at ', &
       time(1:2),':',time(3:4)
  write(un,'(a,a)') 'runname=',runname
  write(un,'(a,i5)')          'run until time =           ', tfinal
  if( isBinary ) then
     write(un,'(a)')          'binary output'
  else
     write(un,'(a)')          'ascii output'
  end if
  write(un,*)
  write(un,*)

end subroutine outputstepinfo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!
! SUBROUTINE: OUTPUTPARAMS -- write out all parameters
!
subroutine outputparams( un )
  integer :: un
 
  character(len=20) :: date
  character(len=20) :: time
 
 
  write(un,'(a)') 'PARAMETER VALUES'
  write(un,'(a)') '----------------'
  write(un,'(a,es12.5)') 'time step                (dt)   ', dt
  write(un,*)
  write(un,'(a,es12.5)') 'viscosity                (mu)   ', mu
  write(un,'(a,es12.5)') 'density                  (rho)  ', rho
  write(un,*) 

  write(un,*)
  write(un,*)
  write(un,'(a)') 'DIMENSIONLESS PARAMETERS'
  write(un,'(a)') '------------------------'
  write(un,'(a,es12.5)') 'reynolds number                   ', Re
  write(un,'(a,es12.5)') 'max alpha^2                       ', asqmax
  write(un,'(a,es12.5)') 'half width of bump                ', xhw
  write(un,*) 

  write(un,'(a)') 'HINDERED TRANSPORT'
  write(un,'(a)') '------------------------'   
  write(un,'(a,es12.5)') 'adp flag                           ', adpf
  write(un,'(a,es12.5)') 'scale for phi_b                    ', kh

  write(un,*)
  write(un,*)
  write(un,'(a)') 'VARIABLE SCALES'
  write(un,'(a)') '---------------'
  write(un,'(a,es12.5)') 'length   (xchar) ', xchar
  write(un,'(a,es12.5)') 'velocity (uchar) ', uchar
  write(un,'(a,es12.5)') 'time     (tchar) ', tchar
 

 ! write(un,*)
 ! write(un,*)
 ! write(un,'(a)') 'BIOCHEMICAL PARAMETERS'
 ! write(un,'(a)') '---------------'
 ! write(un,'(a,es12.5)') 'test concentration '
 


end subroutine outputparams
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!
! SUBROUTINE: LOADPARAMS -- read parameter values from input files displayed 
!                           in the same form as below
!              This is for parameters that will change with simulation
!
subroutine loadparams

!  read(*,*) runname
runname = 'test'
!  read(*,*) tfinal
 tfinal = 1 
mu = 0.04
rho = 1
tmax = 1.0e3
pltdiam = 1.0e-2
dheight = .01
asqmax = 500
xhw =0.5
ytop =0
xsteep =9
ysteep =15
omegaf =1.8
omega =1.001
rh  =1000.0
adpf =0
kh = 0.75
!  read(*,*) mu
!  read(*,*) rho
!  read(*,*) tmax
!  read(*,*) pltdiam 
!  read(*,*) dheight
!  read(*,*) asqmax
!  read(*,*) xhw
!  read(*,*) ytop
!  read(*,*) xsteep
!  read(*,*) ysteep

!  read(*,*) omegaf
!  read(*,*) omega

!  read(*,*) rh
!  read(*,*) adpf
!  read(*,*) kh
end subroutine loadparams
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!subroutine pltbdry(ppp)
!real    :: ppp(1:nyb),y ! platelet boundary values 
!integer :: i          ! for looping
!real    :: intp

!do i=1,nyb
!   y=(i-0.5)*hb
        
   ! Constant
   !ppp(i)=1.0

   !Eckstein
!   ppp(i)=(1 + 202.0*(abs(y-10)/10)**18*(1-(abs(y-10)/10)))  
   !write(*,*)'ppp(i) before',i,ppp(i)
!end do

!intp = 0.5*hb*(ppp(1)+ppp(nyb)) + hb*sum(ppp(2:nyb-1))

!ppp=20.0*ppp/intp
!write(*,*)'intp',intp
!do i=1,nyb
!   write(*,*)'ppp(i)',i,ppp(i)
!end do
!end subroutine pltbdry
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module param_mod
