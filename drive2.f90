module gel_mod
use param_mod
use grid_mod
use bobconst_mod
use gc_mod
use diffusion_mod
use fluidsolve_mod
use clot_mod
implicit none
!!  use param_mod

!write(*,*) 'WTF'
!call init
!  Nxb = 66
!  Nyb = 33
!  Nz =  100
!  call driver(Nxb,Nyb,Nz)


!-----------------------------------------------------------------------------
!Program gelation
!integer, parameter:: Nxb = 66, Nyb = 33, Nz =  100 
integer,parameter:: Nz = 100
double precision, dimension(1:Nz):: hz
double precision:: Nxbd, Nybd, Nzd
double precision:: Lx, Ly,Lz, hx, hy             !domain parameters
double precision:: to, tend, t, dt1, dto, timerecord, mu_alpha,nu
double precision:: kx, ky
double precision:: DD, dalpha, dbeta, dtc
double precision:: Kappaon, Kappaoff, kappaoc
integer::i, j, N1, N2, ix, initialdata, iy, indexrecord, recordidx, dogelation
integer:: dosource,n, jz, flowon,jk
double precision:: Ls, cfl, q
double precision::alphao,alpha2, beta, gamma, L2, L3, kon, kcat, koff,Ro
double precision, dimension(0:Nxb+1):: xo,Rad,x,r1, betat
double precision, dimension(0:Nxb+1):: alpha, alphat
double precision, dimension(0:Nyb+1):: yo, y
double precision, dimension(0:Nxb+1,0:Nyb+1)::D
double precision, dimension(0:Nxb+1,0:Nyb+1):: dis
double precision, dimension(-1:Nxb+2,-1:Nyb+2):: Ulr, Vlr, Utb, Vtb
double precision, dimension(-1:Nyb+2):: yv
!----------------------------------------Z variables---------
double precision:: Zup,ZupD
double precision, dimension(Nxb, Nyb) :: Z2, ADVz2
double precision, dimension(1:Nyb):: BCxRHSz
double precision, dimension(1:Nxb):: BCyRHSz
double precision, dimension(-1:Nxb+2,-1:Nyb+2):: Z2pv
double precision, dimension(0:Nxb+1,0:Nyb+1):: Z2g, Z2p, Z2gh, Z2ph
double precision, dimension(0:Nyb+1)::BCxz
double precision, dimension(0:Nxb+1)::Az, Bz, Dz
!----------------------------------------E variables---------
double precision:: Eup,kat
double precision, dimension(Nxb,Nyb):: E2, ADVe
double precision, dimension(0:Nxb+1):: Ae, Be, De
double precision, dimension(1:Nyb)::   BCxRHSe
double precision, dimension(1:Nxb)::   BCyRHSe 
double precision, dimension(0:Nyb+1):: BCxe
double precision, dimension(0:Nxb+1,0:Nyb+1):: E2g, E2p, E2gh, E2ph
double precision, dimension(-1:Nxb+2,-1:Nyb+2):: E2pv
!----------------------------------------F variables---------
double precision:: Fup,kf, kfs, Fupd, Kappaf,Kappafs
double precision, dimension(Nxb,Nyb):: F, ADVf, Rkf, Rke, Rkz
double precision, dimension(0:Nxb+1):: Af, Bf, Df
double precision, dimension(1:Nyb)::   BCxRHSf
double precision, dimension(1:Nxb)::   BCyRHSf 
double precision, dimension(0:Nyb+1):: BCxf
double precision, dimension(0:Nxb+1,0:Nyb+1):: Fg, Fp1, Fgh, Fph
double precision, dimension(-1:Nxb+2,-1:Nyb+2):: Fpv
!--------------W parameters------------------------------
integer:: m, nt
double precision, dimension(Nz):: z,phi
double precision, dimension(1:Nxb,1:Nyb,1:Nz):: W, ADVw, Funp, ADVzW, Rxn_zW, LW_z
double precision, dimension(0:Nxb+1,0:Nyb+1,1:Nz):: Wg, Wp, Wold, Wnew, Wnewg
double precision, dimension(-1:Nxb+2,-1:Nyb+2,1:Nz):: Wpv 
double precision, dimension(0:Nxb+1,1:Nz):: Aw, Dw, Bw
double precision, dimension(1:Nyb):: BCxRHSw
double precision, dimension(1:Nxb):: BCyRHSw
double precision, dimension(0:Nxb+1,0:Nyb+1):: Wold1
double precision, dimension(0:Nyb+1)::BCxW
double precision, dimension(1:Nxb, 1:Nyb):: ADVtemp
double precision:: kl, kb, mo, Wup, kbd, kld 
double precision:: chi
!-------------R parameters-----------------------------
double precision, dimension(0:Nxb+1, 0:Nyb+1):: R,Rs, Rg, Rsold, Rold, Rgold, rxnr
double precision, dimension(0:Nxb+1, 0:Nyb+1):: Rsavg, Ravg, Rgavg
double precision, dimension(Nxb,Nyb):: ADVr
double precision:: Rup
double precision, dimension(1:Nyb):: BCxRHSR
double precision, dimension(1:Nxb):: BCyRHSR
double precision, dimension(-1:Nxb+2,-1:Nyb+2):: Rsv
double precision, dimension(0:Nxb+1,0:Nyb+1):: Rgh, Rph, Rsh, Rgph
double precision, dimension(0:Nyb+1)::BCxR
double precision, dimension(0:Nxb+1)::AR, BR, DR 
!------------V parameters------------------------
double precision:: Vup
double precision, dimension(1:Nxb,1:Nyb,1:Nz):: V, ADVV, ADVzV, Wz, RXNv, Vold
double precision, dimension(0:Nxb+1,0:Nyb+1,1:Nz):: Vg, Vp, Vnew, Vnewg
double precision, dimension(-1:Nxb+2,-1:Nyb+2,1:Nz):: Vpv 
double precision, dimension(0:Nxb+1,1:Nz):: AV, DV, BV
double precision, dimension(1:Nyb):: BCxRHSV
double precision, dimension(1:Nxb):: BCyRHSV
double precision, dimension(0:Nyb+1)::BCxV
!-------------Theta -------------------------
double precision, dimension(0:Nxb+1, 0:Nyb+1):: theta,thetas, thetag, thetaold
double precision, dimension(0:Nxb+1, 0:Nyb+1):: thetasold, Tsg, Tp, Tsp
double precision, dimension(Nxb,Nyb):: ADVt, RXNt
double precision:: Tup
double precision, dimension(-1:Nxb+2,-1:Nyb+2):: Tsv
double precision, dimension(0:Nyb+1)::BCxT
double precision, dimension(0:Nxb+1)::At, Bt, Dtheta 
!------------Branch---------------------------
double precision, dimension(0:Nxb+1, 0:Nyb+1):: branch,branchs, branchg, branchold
double precision, dimension(0:Nxb+1, 0:Nyb+1):: Bsg, Bp, Bsp, Branchsold
double precision, dimension(Nxb,Nyb):: ADVb, RXNb
double precision:: Bup
double precision, dimension(-1:Nxb+2,-1:Nyb+2):: Bsv
double precision, dimension(0:Nyb+1)::BCxB
double precision, dimension(0:Nxb+1)::Ab, Bb, Db
double precision, dimension(1:Nxb,1:Nyb):: Wint, Vint
double precision, dimension(1:Nz+1):: htrap
double precision, dimension(1:Nz+1):: Vex,Wex
!---------------------------------------------
double precision, dimension(0:Nxb+1,0:Nyb+1):: S10c, S10old

!------------------------Karin Parameters that are needed----------------------------
!real :: fbg(0:nxb+1,0:nyb+1,2)
!real :: asq(1:nxb+1, 1:nyb+1, 2)
!real :: mgtol   = 1e-10
!real :: difftol = 1e-7
!----Saving files----------------
character(100):: Nyb_char, Dtc_char
character(6)::fileend
character(50)::Zfile, Efile, Ffile, Rfile, Rgfile,Tfile, Sfile, Vfile
character(50):: Wfile,Bfile, Thfile, Bgfile, Thgfile, Bsfile, Tsfile, Parsfile
integer:: Wunit=47, Eunit=48, Vunit = 49, Runit=50, Rgunit=51, Tunit=52, Sunit=53
integer:: Thunit=54, Bunit=55, Thgunit=56, Bgunit=57, Bsunit = 58, Tsunit=59
integer:: Z2unit= 61, Funit = 60, Parsunit = 61
External:: DGTSV

Contains
!----------------------------------------------------------------------------
subroutine driver
double precision:: maxu


!--------------------------------------------------------------------------
write(*,*) 'begin next file'
fileend  ='0609_2'
NxbD      = dble(Nxb)
NybD      = dble(Nyb)
NzD      = dble(Nz)

Ls       = 1.00d-2
ZupD     = 1d-9 ! The dimensional upstream prothrombic concentration
DD       = 1d-7 !Dimensional Diffusion 
alphao   = 1d2! 1d2!1d+00 !1d+03 !BC term
dalpha   = alphao*DD*ZupD/(Ls)  ! Dimensional alpha
write(*,*) 'dalpha=', dalpha
alpha    = 0d0 !alphao
beta     = 5d-4!(alpha(1)*(1d0-1d-1*2d0/pi)-1d-1*(1d0-2d0*1d-1/pi))/
!BC term
dbeta    = beta*ZupD
kcat     = 1d-7
Ro       = dalpha/kcat
write(*,*) 'Ro=', Ro, 'kcat=', kcat, 'koffish=', dbeta*kon-kcat
!koff     = 1d-1 !dbeta*kon - kcat 
kon      = 1d6!(kcat+koff)/dbeta
koff     = dbeta*kon-kcat
gamma    = 1d1
Ly       = 2.5d-1 !2.00d0
Lx       = 6d0!d-
D        = 1d0 !DD/(gamma*LS**2)!1.00d-4/(gamma*LS**2d0)!1.00d-07/(gamma*LS**2)
mu_alpha = 6d0 !Sharpness of Alpha transition between 0 and alpha_o
Kappaon  = ZupD*kon*Ls**2d0/DD !Kappaon2 = 0d0 !kon*Ls*Ro/DD
Kappaoc  = (koff+kcat)*Ls**2d0/DD
Kappaoff = -koff*Ro*Ls/(ZupD*DD) !1d8*(koff*Ro*Ls/DD)
kld      = 4.8d1  !I have no idea what numbers these should be!!!
kbd      = 4.8d3  ! I have no idea what numbers these should be!!
Lz       = 1d0
Zup      = 1d0
Eup      = 0d0
Fup      = 4d0
fupd     = 1d-8
kf       = 1d-1*1d2
kfs      = FupD
Kappaf   = kf*Zupd*Ls**2/(DD*Fupd)
Kappafs  = Kfs/FupD
kat      = 1d-1*Ls**2/DD
nu       = 3d0!4d0
kl       = 5d-1!1d0 !Ls**2.*Fupd*kld/DD!1d0!4.8d3*Ls**2/DD
kb       = 1d0!1d1!Ls**2.*Fupd**2.*kbd/DD!5d-1 ! 4.8d1*Ls**2/DD
flowon   = 0
write(*,*) 'kb=', kb, 'kappaf=', kappaf

!Initial Values
W   = 0d0
!R = 0.0
chi = 1d0
!-----------FK Parameters------------------------------------
Ly = 2d-1
Lx = 2d-1
D  = 1d-2
JK = 1  !-1=cheryl's bcs and diffusion, 1=FK bcs and CZ Diff, 2=FK BC and Diff
!The boundary conditions also have to be altered in both Ghost_cells

!-------space------------
hx = Lx/Nxb
hy = Ly/Nyb
write(*,*) 'hx=', hx, 'hy=', hy
call fun_hz
y =  (/ (hy*(iy) - hy/2d0, iy = 0,Nyb+1)/)
xo =  (/ (hx*(ix) - hx/2d0, ix = 0,Nxb+1)/)

z(1) = hz(1)
do i= 2,Nz
   z(i) = z(i-1) + hz(i)
end do
!write(*,*) 'z=', z
x  =  (/ (hx*(dble(ix)) - hx/2d0, ix = 0,Nxb+1)/)

to = 0d0
dt = minval(hz)/8d0 !min(hy/1d3,minval(hz)/8d0,hx/1d3) !minval(hz)/4d0 
tend = 3d-1 !0.8!1d-1!2.5d-1!2d0*dt!8d-1 !1d0!5d0*dt !1d2* dt !5d-2
!tend = tend+dt/
timerecord= 1d-2
!Aaron parameters-------------------------------------------------

!------------Diffusion Initalization (with BC)----------------------


!write(*,*) 'hy = ', hy
!Bz =1d-1 !0d0 ! -hy*1d-1

z2ph = 1d0

do ix = 0,Nxb+1
   Z2ph(ix,0:Nyb+1)   = zup !cos(pi/2d0*y + pi/2d0)*(-Bz(1)*2d0)/(pi)+1-(Bz(1)*2d0)/(pi)
end do
!write(*,*) Z2ph(40,0:Nyb+1)
Z2ph(1:Nxb,Nyb+1)   = Z2ph(1:Nxb,Nyb)
Z2ph(0,0:Nyb+1)    = Z2ph(1,0:Nyb+1)
Z2ph(Nxb+1,0:Nyb+1) = Z2ph(Nxb,0:Nyb+1)
z2                = z2ph(1:Nxb,1:Nyb)
z2p               = z2ph !zup

!Z2g            = Z2p

betat = beta + 3d0/4d0*Z2ph(0:Nxb+1,1) -1d0/8d0*Z2ph(0:Nxb+1,2) +3d0/8d0*Z2ph(0:Nxb+1,0)
Az     = (8d0*betat - 6d0*alpha*hy)/(8d0*betat+3d0*alpha*hy)
Dz     = alpha*hy/(8d0*betat + 3d0*alpha*hy)
Bz     = 0d0

!Time stepping

kx = dt/(2d0*hx**2)    
ky = dt/(2d0*hy**2)

BCxz  = -1d0

L2 =  2.5d0
L3 =  3.5d0

alpha = alphao*(5d-1*tanh(mu_alpha*(xo-L2))-5d-1*tanh(mu_alpha*(xo-L3)))
!write(*,*) 'alphao=', alphao, 'alpha=', alpha

Z2g = Z2ph
!----------------------------------------------------------------
E2   =  eup
E2g  =  eup
E2ph =  eup
E2gh =  eup
E2p  =  eup
Ae   =  1d0
De   =  0d0
BCxe = -1d0
!--------------------------------------------------------
F   = Fup
Fg  = Fup
Fph = Fup
Fgh = Fup
Fp1  = Fup
!------F BC's----
BCxf    = -1d0
Af      = 1d0
Bf      = 0d0
Df      = 0d0
BCyRHSf = 0d0
BCxRHSf = 2d0*D(1,1:Nyb)*kx*Fup
!-----------------W Boudnary Conditions----------------------
Wup     =  0d0
BCxW    =  -1d0 !-1d0 for a Direclet bc on left boundary
Aw      =  1d0
Bw      =  0d0
Dw      =  0d0
BCxRHSW =  0d0 !RHS addition due to BC in x
BCyRHSW =  0d0 ! RHS addition due to BC in y
Wg      =  0d0
Wp      =  0d0
Wold1   =  0d0
Wold    =  0d0
Wnew    =  0d0
!----------------R Initial and Boundary Conditions---------
Rup     =  0d0
BCxR    =  -1d0 !-1d0 for a Direclet bc on left boundary
Ar      =  1d0 ! Bottom Boundary conditions
Dr      =  0d0 ! For Bottom BC
Br      =  0d0 ! For  Bottom BC
BCxRHSr =  0d0
BCyRHSr =  0d0
R       =  0d0
Rg      =  0d0
Rgh     =  0d0
Rsh     =  0d0
Rsv     =  0d0
S10c     =  0d0
Rgold   =  0d0
Rsold   =  0d0
Rold    =  0d0
!-----------V initial and BCs-----------------------------
Vup     =  0d0
BCxV    =  -1d0 !-1d0 for a Direclet bc on left boundary 
AV      =  1d0
BV      =  0d0
DV      =  0d0
BCxRHSV =  0d0 !RHS addition due to BC in x  
BCyRHSV =  0d0 ! RHS addition due to BC in y  
Vg      =  0d0
Vp      =  0d0
Vold    =  0d0
Vnew    =  0d0

!-------------theta initial and BC---------------------
Tup      = 0d0
theta    = 0d0
thetas   = 0d0
thetag   = 0d0
At       = 1d0
Dtheta   = 0d0
Bt       = 0d0
Tsv      = 0d0
Tsg      = 0d0
thetaold = 0d0

!--------Branch initial and BCs----------------

Bup       = 0d0
Branch    = 0d0
Branchs   = 0d0
Branchg   = 0d0
Ab        = 1d0
Db        = 0d0
Bb        = 0d0
Bsv       = 0d0
Bsg       = 0d0
Branchold = 0d0
S10old    = 0d0
call ghost_cells(Z2p, Az, Bz,Dz, Z2g, Nxb, Nyb,1, zup,jk)
!call ghost_cells(Fp, Af, Bf, Df, Fg, Nxb, Nyb, 1 ,fup)

write(Dtc_char,*) dtc
!--------------------------Velocity Information --------------------
Ulr = 0d0
Utb = 0d0
write(*,*) 'flowon=', flowon
write(*,*) 'gamma=', gamma
if (flowon==1) then
   do i = -1, Nxb+2
      Ulr(i,1:Nyb) = Ls**2/DD*gamma*y(1:Nyb)*(1d0-y(1:Nyb)/Ly)!y*(4d0-8d0*y) !1-(1-yv(0:Nyb+1))**2 
   end do
   do i = -1, Nxb+2
      Utb(i,1:Nyb) = Ls**2/DD*gamma*(y(1:Nyb)-hy/2)*(1d0-(y(1:Nyb)-hy/2)/Ly)!1-(1-(yv(0:Nyb+1)-hy/2))**2 
   end do
end if
!write(*,*) Utb(1,1:Nyb+1)
!stop
!Ulr(Nxb+1,1:Nyb)=y(1:Nyb)*(1d0-4d-1*y(1:Nyb))!1-(1-yv(1:Nyb))**2      
Vtb = 0d0
Vlr = 0d0
if (maxval(Utb)*dt/hx >1d0) then
   write(*,*) 'courant number is bigger than 1', maxval(Utb)*dt/hx
   stop
end if
write(*,*) 'courant number=',maxval(Utb)*dt/hx  !this is the value for parabolic flow  
!stop


!----FK parameters part 2-------------------------------
BCxV = 1d0
BCxW = 1d0
BCxR = 1d0
!---------------------------------------------------------

t = to
dto = dt
write(*,*) TRIM(fileend), '.dat'

write(Parsfile, '(25a)') 'Parm_',TRIM(fileend),'.dat'
open(Parsunit, file=Parsfile)

write(Parsunit,20) "kb ", kb
20 format(a, f5.2)
write(Parsunit,20) "kl ", kl
write(Parsunit,20) "MaxU", maxval(Ulr)
write(Parsunit,30) "(kappaf, kappafs)", kappaf,',', kappafs
30 format(a, E7.2,a, E7.2)
write(Parsunit,22) "kat", kat
22 format(a, E7.2)
write(Parsunit,40) "(Nxb, Nyb, Nz)", Nxb,',', Nyb,',', Nz
40 format(a, I3, a, I3,a, I3)
write(Parsunit,41) "(Lx, Ly, Lz)", Lx,',', Ly,',', Lz
41 format(a, E7.2,a, E7.2,a, E7.2)
write(Parsunit,31) "Injury zone ", L2,',', L3
31 format(a, E7.2,a, E7.2)
write(Parsunit,31) "(alpha, beta)", alphao,',', beta
write(Parsunit,31) "(Zup, Fup)", Zup,',', Fup
write(Parsunit,41) "(to, dt, tend)", to,',', dt, ',',tend
write(Parsunit,23) "Is flow on(1=yes, other=no)?", flowon
23 format(a, I1)
write(Parsunit,20) "gamma", gamma
!20 format(a, f5.2)
!30 format(a, f5.2, f5.2)
!40 format(a, f5.2, f5.2, f5.2)
close(Parsunit)

write(Rgfile,'(25a)') 'Rgt_Nxb66_',fileend,'.dat'
open(Rgunit, file = Rgfile, form='binary')

write(Rfile,'(25a)') 'Rt_Nxb66_',fileend,'.dat'
open(Runit, file = Rfile, form='binary')

write(Tfile,'(25a)') 'Time_Nxb66_',fileend,'.dat'
open(Tunit,file=Tfile,form='binary')

write(Sfile,'(25a)') 'Source_Nxb66_',fileend,'.dat'
open(Sunit,file=Sfile,form='binary')

write(Vfile,'(25a)') 'V_Nxb66_',fileend,'.dat'
open(Vunit,file=Vfile,form='binary')

write(Wfile,'(25a)') 'W_Nxb66_',fileend,'.dat'   
open(Wunit, file = Wfile, form='binary')

write(Thfile,'(25a)') 'Theta_Nxb66_',fileend,'.dat'
open(Thunit,file=Thfile,form='binary')

write(Bfile,'(25a)') 'B_Nxb66_',fileend,'.dat'
open(Bunit, file = Bfile, form='binary')

write(Thgfile,'(25a)') 'Thetag_Nxb66_',fileend,'.dat'
open(Thgunit,file=Thgfile,form='binary')

write(Bgfile,'(25a)') 'Bg_Nxb66_',fileend,'.dat'
open(Bgunit, file = Bgfile, form='binary')

write(Tsfile,'(25a)') 'Thetas_Nxb66_',fileend,'.dat'
open(Tsunit,file=Tsfile,form='binary')

write(Bsfile,'(25a)') 'Bs_Nxb66_',fileend,'.dat'
open(Bsunit, file = Bsfile, form='binary')

write(Ffile,'(25a)') 'Fib_Nxb66_',fileend,'.dat'
open(Funit, file = Ffile, form='binary')

write(Efile,'(25a)') 'E2_Nxb66_',fileend,'.dat'
open(Eunit, file = Efile, form='binary')

write(Zfile,'(25a)') 'Z_Nxb66_',fileend,'.dat'
open(Z2unit,file=Zfile, form='binary')


do j = 1,Nyb
   write(Bsunit)  Branchs(1:Nxb,j)
   write(Tsunit)  Thetas(1:Nxb,j)
   write(Bunit)   Branch(1:Nxb,j)
   write(Bgunit)  Branchg(1:Nxb,j)
   write(Thunit)  Theta(1:Nxb,j)
   write(Thgunit) Thetag(1:Nxb,j)
   write(Runit)   R(1:Nxb,j)
   write(Rgunit)  Rg(1:Nxb,j)
   write(Eunit)   E2g(1:Nxb,j)
   write(Funit)   Fg(1:Nxb,j)
   write(Z2unit)   Z2(1:Nxb,j)
   write(Sunit)   S10c(1:Nxb,j)
end do

do jz = 1,Nz
   do j = 1,Nyb
       write(Wunit)   W(1:Nxb,j,jz)
       write(Vunit)   V(1:Nxb,j,jz)
    end do
end do
nt = 0
write(Tunit)  to
write(*,*) 'dt=', dt
dogelation = 1
dosource   = 0
!n = 0
do while  (t< tend)
!------------Fluid Solver------------------------
   fbg(:,:,1) = 0.0
   fbg(:,:,2) = 0.0

     ! set the background force - fns initialezed in fluidinit                                                                      
   fns(1:nxb+1,1:nyb+1,:) = 0.0  
   maxu = maxval( abs(u) )
   call new_alpha(asq)
   call fluidstep
!   if (t+dt > tend) then
!      dt = tend-t
!      kx = dt/(2d0*hx**2)
!      ky = dt/(2d0*hy**2)
!   end if
      t = t+dt
!   write(*,*) 't=', t, 'dt=', dt   
   !   write(*,*) 't=', t, 'dt=', dt
if (dosource==1) then
         !----------Z-BC's-----------------------------
   alphat = alpha*(1d0-exp(-nu*(t)))
   !write(*,*) 'alphat=', alphat
   betat   = beta + 3d0/4d0*Z2g(0:Nxb+1,1) -1d0/8d0*Z2g(0:Nxb+1,2)&
           + 3d0/8d0*Z2g(0:Nxb+1,0)
   Az      = (8d0*betat - 6d0*alphat*hy)/(8d0*betat+3d0*alphat*hy)!1d0  
   Dz      = alphat*hy/(8d0*betat + 3d0*alphat*hy)
   BCxRHSz = 2d0*D(1,1:Nyb)*kx*zup
   BCyRHSz = -D(1:Nxb,1)*ky*Bz(1:Nxb)

   ADVz2=0d0
   Rkz = 0d0
!   write(*,*) 'Nxb =',Nxb, 'Nyb=',Nyb,'kx=', kx, 'ky=', ky, 'Zup=', Zup,'BCxz=', BCxz, 'BCxRHSz=', BCxRHSz
!   write(*,*) 'Az=', Az, 'Dz=', Dz, 'Bz=', Bz, 'BCyRHSz=', BCYRHSz, 'D=', D, 'Z2g=', Z2g, 'Z2'
!   write(*,*) Z2, 'ADVz2=', ADVz2, 'Rkz=', Rkz
!   STOP
  
   call diffusion_subr(Nxb,Nyb,kx,ky,Zup,BCxz,BCxRHSz,Az,Dz,Bz,BCyRHSz,D,Z2g,Z2,ADVz2,Rkz)
   call ghost_cell_vec(Z2g,Z2pv, Az,Bz,Dz,Nxb,Nyb,1,Zup,jk)
   call advection(Z2pv,Ulr,Utb,Vlr,Vtb, dt,hx,hy,Nxb,Nyb, ADVz2)

   Z2g = Z2pv(0:Nxb+1,0:Nyb+1)
   Z2g(1:Nxb,1:Nyb) = Z2g(1:Nxb,1:Nyb) + ADVz2
   !write(*,*) Z2g
   call ghost_cells(Z2g,Az,Bz,Dz,Z2p,Nxb,Nyb,1,zup,jk)
   
   !write(*,*) Z2p
   !stop
   Z2g = Z2p
   Z2  = Z2g(1:Nxb,1:Nyb)
 !  write(*,*) 'Z2=', maxval(Z2), minval(Z2)
!   write(*,*) 'Z2=', maxloc(Z2), minloc(Z2)
!   Z2g = 1d0
!   Z2  = 1d0
   !---------------------Advect and Diffusion E2------------
  !---------------E2----------------------------

   Be = hy*alphat*(3d0/4d0*Z2g(0:Nxb+1,1)-1d0/8d0*Z2g(0:Nxb+1,2) &
        + 3d0/8d0*Z2g(0:Nxb+1,0))/betat
   !write(*,*) 'Be=', Be
   BCxRHSe = 2d0*D(1,1:Nyb)*kx*Eup
   BCyRHSe = D(1:Nxb,1)*ky*Be(1:Nxb)
   BCxe    = -1d0
   !------------------------

   RKe = -dt*kat*E2(1:Nxb,1:Nyb)
   ADVe = 0d0
   call ghost_cell_vec(E2g,E2pv, Ae,Be,De,Nxb,Nyb,1,Eup,jk)
   call advection(E2pv,Ulr,Utb,Vlr,Vtb, dt,hx,hy,Nxb,Nyb, ADVe)
   call ghost_cells(E2p,Ae,Be,De,E2g,Nxb,Nyb,1,Eup,jk)                                                      
   call diffusion_subr(Nxb,Nyb,kx,ky,Eup,BCxe,BCxRHSe,Ae,De,Be,BCyRHSe,D,E2g,E2,ADVe,RKe)
   !write(*,*) 'E2=', maxval(E2(1:Nxb,1:Nyb)), minval(E2(1:Nxb,1:Nyb))
   call ghost_cell_vec(E2g,E2pv, Ae,Be,De,Nxb,Nyb,1,Eup,jk)
 ! write(*,*) 'E2=', maxval(E2(1:Nxb,1:Nyb)), minval(E2(1:Nxb,1:Nyb))
   call advection(E2pv,Ulr,Utb,Vlr,Vtb, dt,hx,hy,Nxb,Nyb, ADVe)
   !E2g = E2pv(0:Nxb+1,0:Nyb+1)
   E2g(1:Nxb,1:Nyb) = E2g(1:Nxb,1:Nyb) + ADVe
  call ghost_cells(E2g,Ae,Be,De,E2p,Nxb,Nyb,1,Eup,jk)
   
   E2g = E2p
   E2  = E2g(1:Nxb,1:Nyb)
   write(*,*) 'E2=', maxval(E2), minval(E2)
!  write(*,*) 'E2=', maxloc(E2), minloc(E2)
!  write(*,*) 'E2=', E2(1:Nxb,1)
!  write(*,*) 'E2=', E2g(0:Nxb+1,0:1)
   !---------------------Advect and Diffusion F2------------
   BCxRHSf = 2d0*D(1,1:Nyb)*kx*Fup
   ADVF = 0d0
   RKF  = 0d0
!   write(*,*) 'Fg(35,1)=', Fg(35,1)
   call diffusion_subr(Nxb,Nyb,kx,ky,Fup,BCxf,BCxRHSf,Af,Df,Bf,BCyRHSf,D,Fg,F,ADVf, RKf)!
!   write(*,*) 'maxF=', maxval(Fg(1:Nxb,1:Nyb)), maxloc(Fg(1:Nxb,1:Nyb))
!   write(*,*) 'minF=',minval(Fg(1:Nxb,1:Nyb)),'loc=', minloc(Fg(1:Nxb,1:Nyb))
!   write(*,*) 'Fg(35,1)=', Fg(35,1)
   call ghost_cell_vec(Fg,Fpv, Af,Bf,Df,Nxb,Nyb,1,Fup,jk)
   call advection(Fpv,Ulr,Utb,Vlr,Vtb, dt,hx,hy,Nxb,Nyb, ADVf)

   Fg = Fpv(0:Nxb+1,0:Nyb+1)
   Fg(1:Nxb,1:Nyb) = Fg(1:Nxb,1:Nyb) + ADVf
   
   call ghost_cells(Fg,Af,Bf,Df,Fp1,Nxb,Nyb,1,Fup,jk)
   call RungaKutta(E2g, Fp1, kappaf, kappafs, Nxb, Nyb, dt, Rkf)
   call ghost_cells(Fp1,Af,Bf,Df,Fg,Nxb,Nyb,1,Fup,jk)

   F  = Fg(1:Nxb,1:Nyb)
!   write(*,*) 'F=', maxval(F(1:Nxb,1:Nyb)), minval(F(1:Nxb,1:Nyb))
!   write(*,*) 'Floc =', maxloc(F(1:Nxb,1:Nyb)), minloc(F(1:Nxb,1:Nyb))
 !----------------------------------------------------------------        
   if (minval(F) < -1d-12) then
      write(*,*) 'Error - F is negative', minval(F), 't=', t
      stop
   end if
   !write(*,*) 'Fgend=', Fg(0:Nxb+1,0:1)
S10c = kappaf*E2g*Fg/(kappafs+Fg)
end if
!S10 = 0d0
if (dosource==0) then
   if (t <.1) then
      !write(*,*) 't1=', t
      do ix = 0,Nxb+1

         S10c(ix,1:Nyb) = 64d0*(-5d-1*(tanh(1d2*(xo(ix)-1.5d-1)))+5d-1*(tanh(1d2*(xo(ix)-5d-2))))*(1-exp(-8d0*t))
      end do
   else
!write(*,*) 't2=', t
      do ix = 0,Nxb+1

         S10c(ix,1:Nyb) = 64d0*(-5d-1*(tanh(1d2*(xo(ix)-1.5d-1)))+5d-1*(tanh(1d2*(xo(ix)-5d-2))))*(1-exp(-8d0*.1))*exp(-8d0*(.1-t))
!if (t< .2) then
!S10(ix, 1:Nyb) = !alphao*(5d-1*tanh(100*(xo(ix)-.15))-5d-1*tanh(mu*(xo(ix)-L3)))*(1-exp(-8d0*t)) &
                !*(-2.5d-1*(tanh(1d1*(y(1:Nyb)/Ly-5d-1))-1d0))
!endif 
!if (t>.2) then
!S10(ix, 1:Nyb) = alphao*(5d-1*tanh(mu*(xo(ix)-L2))-5d-1*tanh(mu*(xo(ix)-L3)))*exp(-8d0*t) &
!                *(-2.5d-1*(tanh(1d1*(y(1:Nyb)/Ly-5d-1))-1d0))
!endif
!   if(abs(x(ix)) < 5d-2) then
!      S10(ix,1:Nyb) = 8d0*8d0*exp(-8d0*t)*5d-1*(1+cos(pi*x(ix)/5d-2))
!      end if
      end do
   end if
end if
!S10 = kappaf*E2g*Fg/(kappafs+Fg)
!write(*,*) S10(1:Nxb,1)
!stop
 !----- W Advection in Z-------------------------------------------------
if (dogelation ==1) then
   call b_w(Wold(1:Nxb,1:Nyb,1:Nz),Vold(1:Nxb,1:Nyb,1:Nz), Rold(1:Nxb,1:Nyb),&
            hz, Nxb, Nyb, Nz, dt, z, kl, kb, ADVzW, Funp,1)
   !call LaxWendroff(Wold(1:Nxb,1:Nyb,1:Nz),Nxb,Nyb,Nz,kl,kb,dt,Rold(1:Nxb,1:Nyb),hz,z,LW_z)

   !phi = 1d0
   !phi(1) = 0d0 !5d-1*tanh(2d1*(z-15d-2))+5d-1
   !phi(2) = 0d0
   !phi(3) = 1d-1
   !phi(4) = 2d-1
   !phi(5) = 4d-1
   !phi(6) = 8d-1
   
   !ADVzW(1:Nxb,1:Nyb,1) = LW_z(1:Nxb,1:Nyb,1)
   !do i = 1,Nz-1
   !   ADVzW(1:Nxb,1:Nyb,i) = ADVzW(1:Nxb,1:Nyb,i)*phi(i) + LW_z(1:Nxb,1:Nyb,i)*(1d0-phi(i))
   !end do

!write(*,*) 'B_w: ADVzW=', ADVzW(1:Nxb,1,Nz)
!stop
   !write(*,*) 'B_w: ADVzW=', ADVzW(1:Nxb,1,Nz)
   !---Check CFL----------
   q=maxval(abs(Funp))

   cfl = q*dt/minval(hz)

   if (cfl>1d-10) then
      if (cfl <=2d-1) then
         dt1 = dt*(4.9d-1/cfl)
         dt = min(dt1,dto);
      else if (cfl >=5d-1) then
         dt = dt*(4.9d-1)/cfl
      end if
   else
      dt = dto
   end if
   
   kx = dt/(2d0*hx**2.) 
   ky = dt/(2d0*hy**2.)

   !-----Update W at z=1 ---------------------

      !This line of coded should be added somehow to the RXNW subroutine so I don't have to 
      !change the reaction equation in 2 places

   Wnew(1:Nxb,1:Nyb,Nz) = Wold(1:Nxb,1:Nyb,Nz) -dt*ADVzW(1:Nxb,1:Nyb,Nz) + dt*(-kl*z(Nz)*chi*Rgold(1:Nxb,1:Nyb)**2.) &
                      + dt*(kb/2d0*z(Nz)**2.*Rold(1:Nxb,1:Nyb)**3.) & 
                      - dt*(kb/2d0*(z(Nz)*Rold(1:Nxb,1:Nyb)**3. + 3d0*z(Nz)*chi*Rsold(1:Nxb,1:Nyb)*Rgold(1:Nxb,1:Nyb)**2. &
                      + z(Nz)*chi*Rgold(1:Nxb,1:Nyb)**3.))

  ! write(*,*) 'Wnew=', Wnew(1:Nxb,1,Nz)
   !------Update R using W @ z=1-----------------
   !--Rg n+1--------------
   Rg(1:Nxb,1:Nyb)  = -Wnew(1:Nxb,1:Nyb,Nz)
!   write(*,*) 'Wnew(Nz)=', Wnew(1:Nxb,1,Nz)
!   write(*,*) 'Wold(Nz)=', Wold(1:Nxb,1,Nz)
!   write(*,*) 'ADVzW(1:Nxb,1,Nz)=', ADVzW(1:Nxb,1,Nz)
!   write(*,*) 'Rgold(1:Nxb,1)=', Rgold(1:Nxb,1)
!   write(*,*) 'Rsold(1:Nxb,1)=', Rsold(1:Nxb,1)
!   write(*,*) 'Rold(1:Nxb,1)=', Rold(1:Nxb,1)
!   write(*,*) 'z=', z(Nz)
!   write(*,*) 'kb=', kb
!   write(*,*) 'kl=', kl, 'chi=', chi
!   write(*,*) 't=', t
   !   write(*,*) 'Rold', Rold(33,1), - dt*(kb/2d0*(z(Nz)*Rold(33,1)**3.))
   call ghost_cell_vec(Rsold, Rsv,Ar, Br, Dr, Nxb, Nyb, 1, Rup,jk)
   call advection(Rsv, Ulr, Utb, Vlr, Vtb, dt, hx, hy, Nxb, Nyb, ADVr)
   
   call RXN_R(Nxb,Nyb, Rold, Rgold, Rsold, S10c, kl, kb, chi, dt,rxnr)
   !write(*,*) 'RXNR=', rxnr(1:Nxb,1)
!    Rold(1:Nxb,1:Nyb)= Rold(1:Nxb,1:Nyb) + RXNr(1:Nxb,1:Nyb)
   call ghost_cells(Rold, Ar, Br, Dr, Rph, Nxb, Nyb, 1, Rup,jk)
      !We need a single ghost cell for Rgel (old and new) for the RHS in Diffusion Subroutine
   call ghost_cells(Rg, Ar, Br, Dr, Rgph, Nxb, Nyb,1, Rup,jk)
   Rg = Rgph

   call ghost_cells(Rgold,Ar,Br,Dr,Rgph, Nxb, Nyb,1, Rup,jk)
   Rgold = Rgph

   call diffusion_subr_r(Nxb,Nyb,kx,ky,Rup,BCxR,BCxRHSR,AR,DR,BR,BCyRHSR,D,Rg,Rph,Rgold,RXNr,ADVr,dt,jk)      
!   write(*,*) 'RxnR=', RXNr(1:Nxb,1)

  ! call RXN_R(Nxb,Nyb, Rold, Rgold, Rsold, S10, kl, kb, chi, dt,rxnr)
  ! Rph(1:Nxb,1:Nyb) = Rph(1:Nxb,1:Nyb) + RXNr(1:Nxb,1:Nyb)
   R  = Rph
   Rs = R-Rg
      !write(*,*) 'R(1:Nxb,2)=', R(1:Nxb,1)
      !      write(*,*) 'R=', R(1:Nxb,20)
      !----Diffusion and Advect W in x,y -------
     
   do m= 1,Nz
      Wg(1:Nxb,1:Nyb,m) = Wold(1:Nxb,1:Nyb,m) - z(m)*(Wold(1:Nxb,1:Nyb,Nz)+Wnew(1:Nxb,1:Nyb,Nz))/2d0
   end do

   call ghost_cell_vec(Wg,Wpv,Aw,Dw,Bw,Nxb,Nyb,Nz,Wup,jk)

   do m = 1, Nz-1

      call advection(Wpv(-1:Nxb+2,-1:Nyb+2,m), Ulr, Utb, Vlr,Vtb, dt, hx, hy, Nxb, Nyb, ADVtemp)
         
      ADVw(1:Nxb,1:Nyb,m) = ADVtemp
   end do

   ADVw(1:Nxb,1:Nyb,Nz) = 0d0

   Rgavg = Rgold !(Rg+Rgold)/2d0 !Rgold
   Ravg  = Rold !(R+Rold)/2d0 !Rold
   Rsavg = Rsold !(Rs+Rsold)/2d0 !Rsold
   Rxn_zW = 0d0
   call RXNW(Nxb,Nyb,Nz,chi,Rgavg(1:Nxb,1:Nyb), Ravg(1:Nxb,1:Nyb), Rsavg(1:Nxb,1:Nyb),&
             z, kl, kb, dt,Rxn_zW) 

   Wp(1:Nxb,1:Nyb,1:Nz) = Wold(1:Nxb,1:Nyb,1:Nz)

   call ghost_cells(Wp,  Aw, Bw, Dw,Wg,    Nxb,Nyb,Nz,Wup,jk)
   call ghost_cells(Wnew,Aw, Bw, Dw,Wnewg, Nxb,Nyb,Nz,Wup,jk)
   call diffusion_subr_W(Nxb,Nyb,Nz,kx,ky,Wup,BCxW,BCxRHSW,Aw,Bw, Dw,BCyRHSW,D,Wg,Wp, z,&
                        Wold1, ADVw, Rxn_zW, ADVzW, Wnewg,dt,jk)
      
      !Interpolation at z=0-------------------
      Wg(1:Nxb,1:Nyb,1) = Wg(1:Nxb,1:Nyb,2)*((z(1)-0)*(z(1)-z(3))*(z(1)-z(4)))/((z(2)-0)*(z(2)-z(3))*(z(2)-z(4)))&
                      + Wg(1:Nxb,1:Nyb,3)*((z(1)-0)*(z(1)-z(2))*(z(1)-z(4)))/((z(3)-0)*(z(3)-z(2))*(z(3)-z(4))) &
                      + Wg(1:Nxb,1:Nyb,4) *((z(1)-0)*(z(1)-z(2))*(z(1)-z(3)))/((z(4)-0)*(z(4)-z(2))*(z(4)-z(3)))    
   Wnew = Wg
   W    = Wg(1:Nxb,1:Nyb,1:Nz)
  
      !----------------V equation -------------------------
   call b_w(W,Vold, R(1:Nxb,1:Nyb), &
        hz, Nxb, Nyb, Nz, dt, z, kl, kb, ADVzV, Wz,2)
   Vg(1:Nxb,1:Nyb, 1:Nz) = Vold

   call ghost_cell_vec(Vg,Vpv,AV,DV,BV,Nxb,Nyb,Nz,Vup,jk)

   do m = 1, Nz
      call advection(Vpv(-1:Nxb+2,-1:Nyb+2,m), Ulr, Utb, Vlr,Vtb, dt, hx, hy, Nxb, Nyb, ADVtemp)
      ADVV(1:Nxb,1:Nyb,m) = ADVtemp
      !if(maxval(abs(ADVV)) >1d-12) then
      !   write(*,*) 'Error ADVV'
      !end if
   end do
      
   Vp (1:Nxb,1:Nyb,1:Nz) = Vold
   call ghost_cells(Vp,  AV, BV, DV,Vg, Nxb,Nyb,Nz,Vup,jk)

   call RXN_V(Wnew(1:Nxb,1:Nyb,1:Nz),Wz,z,R(1:Nxb,1:Nyb),S10c(1:Nxb,1:Nyb), kb, Nxb, Nyb, Nz,dt, RxnV)
   call diffusion_subr_V(Nxb,Nyb,Nz,kx,ky,Vup,BCXV,BCxRHSV,AV,BV, DV,BCyRHSV,D,Vg,Vp, &
                           ADVV, RxnV,ADVzV,dt,jk)
   !write(*,*) 'Vg(1:nx,1:Nyb,1) =', Vg(50,1,1)
   !Vg(1:Nxb,1:Nyb,1) = Vg(1:Nxb,1:Nyb,2)*((z(1)-0)*(z(1)-z(3))*(z(1)-z(4)))/((z(2)-0)*(z(2)-z(3))*(z(2)-z(4)))&
   !                + Vg(1:Nxb,1:Nyb,3)*((z(1)-0)*(z(1)-z(2))*(z(1)-z(4)))/((z(3)-0)*(z(3)-z(2))*(z(3)-z(4))) &
   !                + Vg(1:Nxb,1:Nyb,4)*((z(1)-0)*(z(1)-z(2))*(z(1)-z(3)))/((z(4)-0)*(z(4)-z(2))*(z(4)-z(3)))
   V = Vg(1:Nxb,1:Nyb,1:Nz)
   !write(*,*) 'Vg(1:nx,1:Nyb,1) =', Vg(50,1,1)
   !write(*,*) 'max V=', maxval(V)
      !-------------------theta's and B's-------------------------
        !Need to find the integral of V (Vint) and W over z first----
   htrap(1)    = hz(1)/2d0
   htrap(2:Nz) = (hz(1:Nz-1)+hz(2:Nz))/2d0
   htrap(Nz+1) = hz(Nz)/2d0
   !write(*,*) 'htrap=', htrap
   !stop
   do i=1, Nxb
      do j = 1,Nyb
         Vex(1)      = 0d0
         Vex(2:Nz+1) = V(i,j,1:Nz)
         Wex(1)      = 0d0
         Wex(2:Nz+1) = W(i,j,1:Nz)
         Vint(i,j)   = dot_product(htrap, Vex)
         Wint(i,j)   = dot_product(htrap, Wex)
      end do
   end do
   
   thetas  = 0d0
   branchs = 0d0
   branchs(1:Nxb,1:Nyb) = W(1:Nxb,1:Nyb,Nz) - 2d0*Wint
   thetas(1:Nxb,1:Nyb)  = Vint  + 2d0*branchs(1:Nxb,1:Nyb)
   !write(*,*) 'Thetas-diffThetas=', Thetas(1:Nxb,1:Nyb) - (Vint + 2d0*W(1:nx,1:Nyb,Nz) - 2d0*Wint)
   !branchs(1:Nxb,1:Nyb) = W(1:Nxb,1:Nyb,Nz) - 2d0*Wint
   !write(*,*) 'branchs=', branchs(1:Nxb,1)     
   write(*,*) 't=', t
   write(*,*) 'Wint = ', maxval(Wint), minval(Wint)
   write(*,*) 'Vint = ', maxval(Vint), minval(Vint)
   !write(*,*) 'W = ', W(1:Nxb,1,Nz)
   RXNt  = dt*(S10c(1:Nxb,1:Nyb))!+S10old(1:Nxb,1:Nyb))/2d0
   Ravg  = R  !Rold !(R+ Rold)/2d0
   Rsavg = Rs !Rsold!(Rs + Rsold)/2d0
   Rgavg = Rg !Rgold !(Rg+ Rgold)/2d0
   RXNb  = dt*kb/6d0*(Ravg(1:Nxb,1:Nyb)**3. - chi*(3d0*Rsavg(1:Nxb,1:Nyb)*Rgavg(1:Nxb,1:Nyb)**2.+Rgavg(1:Nxb,1:Nyb)**3.))
   Tsp   = 0d0
   Bsp   = 0d0
   Tsp(1:Nxb,1:Nyb)  = (thetas(1:Nxb,1:Nyb))  !  + thetasold)/2d0
   Bsp(1:Nxb,1:Nyb)  = (branchs(1:Nxb,1:Nyb)) !+ branchsold)/2d0

   call ghost_cell_vec(Tsp, Tsv, At, Bt, Dtheta, Nxb, Nyb, 1, Tup,jk)
   call advection(Tsv, Ulr, Utb, Vlr, Vtb, dt, hx, hy, Nxb, Nyb, ADVt)
   !write(*,*) 'ADVt=', minval(ADVt)
   call ghost_cell_vec(Bsp, Bsv, Ab, Bb, Db, Nxb, Nyb, 1, Bup,jk)
   call advection(Bsv, Ulr, Utb, Vlr, Vtb, dt, hx, hy, Nxb, Nyb, ADVb)
   !write(*,*) 'ADVb=', minval(ADVb)
   Tsp = 0d0
   Bsp = 0d0
   
   Tsp(1:Nxb,1:Nyb) = (thetas(1:Nxb,1:Nyb)+thetasold(1:Nxb,1:Nyb))/2d0
   Bsp(1:Nxb,1:Nyb) = (branchs(1:Nxb,1:Nyb))!+branchsold(1:Nxb,1:Nyb))/2d0
   !if (t==dt) then
   !   Thetaold(1:Nxb,1:Nyb) = Thetas(1:Nxb,1:Nyb)
      
   !endif
   !Ghost cells for the Sol theta and Branch--------------
   call ghost_cells(Tsp,  At, Bt, Dtheta,Tsg, Nxb,Nyb,1,Tup,jk)
   call ghost_cells(Bsp,  Ab, Bb, Db,    Bsg, Nxb,Nyb,1,Bup,jk)
   
 !  write(*,*)'Theta'
   call diffusion_subr_tb(Nxb, Nyb, Thetaold,  theta,  Tsg, RXNt,kx,ky,D,ADVt,At,Bt,Dtheta,Tup,jk)
   call diffusion_subr_tb(Nxb, Nyb, Branchold, branch, Bsg, RXNb,kx,ky,D,ADVb,Ab,Bb,Db,    Bup,jk)
!      write(*,*) 'ADVt=', minval(ADVt), maxval(ADVt)
!      write(*,*) 'RXNt=', minval(RXNt), maxval(RXNt)

   Thetag(1:Nxb,1:Nyb)  = Theta(1:Nxb,1:Nyb)  - Thetas(1:Nxb,1:Nyb)
   
   write(*,*) 'Thetag = ', maxval(Thetag), minval(Thetag)
   !if (t.gt.100d0*dt) then
   !   stop
   !end if
   !if (minval(Thetag(1:Nxb,1:Nyb)) < -1d-12) then
   !   write(*,*) 'Thetag is negative', minval(Thetag(1:Nxb,1:Nyb))
      !write(*,*) 'Source rate*dt=', 2d0*S10*dt
   !   write(*,*) ' maxval(ADVt) =', maxval(ADVt) 
   !   STOP
   !end if
   Branchg(1:Nxb,1:Nyb) = Branch(1:Nxb,1:Nyb) - Branchs(1:Nxb,1:Nyb)
   !if (minval(Branchg) <-1d-10) then
   !   write(*,*) ' Branching is going negative!!!!'
   !STOP   
   !end if
   !if (minval(Thetag(1:Nxb,1:Nyb)) < -1d-10) then
   !   write(*,*) ' Theta g is negative!!!!', minval(Thetag(1:nx,1:Nyb))
   !   stop
   !end if
   write(*,*) 'Branchg = ', maxval(Branchg), minval(Branchg)
!   write(*,*) 'R=', minval(R(1:Nxb,1:Nyb)), 'Rs=', minval(Rs(1:Nxb,1:Nyb)), 'Rg=', minval(Rg(1:Nxb,1:Nyb))
!   write(*,*) 'Rgmax=', maxval(Rg(1:Nxb,1:Nyb)), 'min S10=', minval(S10(1:Nxb,1:Nyb)), minloc(S10(1:Nxb,1:Nyb))
!    write(*,*) 't=',t  
!   write(*,*) 'Bg=', minval(Branchg(1:Nxb,1)), 'Bs=', Branchs(23,1), 'B=', Branch(23,1)
!   write(*,*) 'Thetag=', Thetag(23,1), 'Thetas=', Thetas(23,1), 'Theta=', Theta(23,1)
      !----------------------------------------------------------
end if
   if (t > timerecord) then
      do recordidx = 1, Nyb
         write(Tsunit)  Thetas(1:Nxb,recordidx)
         write(Thunit)  Theta(1:Nxb,recordidx)
         write(Thgunit) Thetag(1:Nxb,recordidx)
         write(Bsunit)  Branchs(1:Nxb,recordidx)
         write(Bunit)   Branch(1:Nxb,recordidx)
         write(Bgunit)  Branchg(1:Nxb,recordidx)
         write(Runit)   R(1:Nxb,recordidx)
         write(Rgunit)  Rg(1:Nxb,recordidx)
         write(Sunit)   S10c(1:Nxb,recordidx)
         write(Funit)   Fg(1:Nxb, recordidx)
         write(Eunit)   E2g(1:Nxb,recordidx)
         write(Z2unit)  Z2(1:Nxb,recordidx)
      end do
      do jz =1,Nz
         do recordidx= 1,Nyb
            if (jz == Nz) then
 !           write(*,*) 'jz=', jz, 'recordidx=', recordidx, W(6,recordidx, jz)
            end if
            write(Wunit)   W(1:Nxb,recordidx,jz)
            write(Vunit)   V(1:Nxb,recordidx,jz)
         end do
      end do
      write(Tunit) t
      !write(*,*) 'R=', R(20,1), 'Rs=', Rs(20,1), 'Rg=', Rg(20,1), 'S10=', S10(20,1)
      !write(*,*) 't=',t,  'B=', Branch(20,1), 'Theta=', Theta(20,1)
      !write(*,*) 'vnew=', V(1,1,1:Nz), 'end'
      timerecord = timerecord + 1d-2
   end if

      !-----------------------------------------------------------
   Thetaold  = Theta
   Branchold = Branch
   Wold1  = Wg(0:Nxb+1,0:Nyb+1,Nz)
   Wold   = Wg
   Rold   = R
   Rsold  = Rs
   Rgold  = Rg      
   Vold   = V
   S10old = S10c
   nt     = nt+1
!   t      = t+dt
   Thetasold   = Thetas
   Branchsold  = Branchs
!   n = n+1
enddo
!write(*,*) 'S10=', S10(1:Nxb,1)
!write(*,*) 'Fg=', Fg(1:Nxb,1), 'E2g=', E2g(1:Nxb,1)
!if (mod(nt-1,100) .ne.0) then
   do recordidx = 1, Nyb
      write(Bsunit)  Branchs(1:Nxb,recordidx)
      write(Tsunit)  Thetas(1:Nxb,recordidx)
      write(Thunit)  Theta(1:Nxb,recordidx)
      write(Thgunit) Thetag(1:Nxb,recordidx)
      write(Bunit)   Branch(1:Nxb,recordidx)
      write(Bgunit)  Branchg(1:Nxb,recordidx)
      write(Runit)   R(1:Nxb,recordidx)
      write(Rgunit)  Rg(1:Nxb,recordidx)
      write(Sunit)   S10c(1:Nxb,recordidx)
      write(Funit)   Fg(1:Nxb,recordidx)
      write(Eunit)   E2g(1:Nxb,recordidx)
      write(Z2unit)  Z2g(1:Nxb,recordidx)
   end do
   do jz = 1,Nz
         do recordidx = 1,Nyb
            write(Wunit)   W(1:Nxb,recordidx,jz)
            write(Vunit)   V(1:Nxb,recordidx,jz)
         end do
      end do
      write(tunit) t
!end if
close(Bsunit)
close(Tsunit)
close(Bunit)
close(Bgunit)
close(Thunit)
close(Thgunit)
close(Runit)
close(Rgunit)
close(tunit)
close(Sunit)
close(Vunit)
write(*,*) 'D=', d(1,1), 'dt=', dt, 'hx=', hx, 'kx=', kx
close(Wunit)
close(Z2unit)
write(*,*) 't=',t 
!write(Nyb_char,*) Nyb !CEILING(dt*1e6)
!write(Zfile,'(25a)') 'Z2_N'//trim(adjustl(Nyb_char))//'_'//trim(adjustl(Dtc_char))//'_robin_2.dat'
!open(Zunit, file = Zfile, form='binary')
!do j= 1,Nyb
!   write(Zunit) Z2(1:Nxb,j)
!end do
!write(*,*) 'Z2', Z2(3,1:Nyb)
!close(Zunit)

!write(Efile, '(25a)') 'E2_N'//trim(adjustl(Nyb_char))//'_'//trim(adjustl(Dtc_char))//'_robin_2.dat'
!open(Eunit, file = Efile, form='binary')
!do j= 1,Nyb
!   write(Eunit) E2(1:Nxb,j)
!end do
!write(*,*) 'Z2', Z2(3,1:Nyb)
close(Eunit)

!write(Ffile, '(25a)') 'F_N'//trim(adjustl(Nyb_char))//'_'//trim(adjustl(Dtc_char))//'_robin_2.dat'
!open(Funit, file = Ffile, form='binary')
!do j= 1,Nyb
!   write(Funit) F(1:Nxb,j)
!end do
!write(*,*) 'Z2', Z2(3,1:Nyb)
close(Funit) 
!write(*,*) S10(1:Nxb,1)
!write(*,*) 'nt=', nt
!write(*,*) 'Rg=', Rg(1:Nxb,1:3)
!contains
end subroutine driver
!------------------------------------------------------------------
subroutine b_w(W,V, R, hz, Nxb, Nyb, Nz, dt, z, kl, kb, Fz, fpbw, type)
             
integer:: Nxb, Nyb, Nz, j, m,i, type
double precision, dimension(Nxb, Nyb,1:Nz):: W, Fz, Fpbw, V
double precision, dimension(Nxb,Nyb, 1:Nz):: fo
double precision, dimension(Nxb, Nyb):: R
double precision, dimension(1:Nz):: hz,z
double precision, dimension(1:Nz):: alpha, beta, gamma
double precision:: kb, kl, dt
Fz = 0d0

if (type==1) then
   
   call function_W(W, fo, hz, Nxb, Nyb, Nz,kl,kb,z,R)
   !write(*,*) 'f=', fo(1:Nxb,1,Nz)
   !write(*,*) 'R=', R(1:Nxb,1)
!   write(*,*) 'fo=', fo(20,1,1:Nz)
   !j = -1
   !call function_W(W, j, fm1, hz, Nxb, Nyb, Nz,kl,kb,z,R)

   !j = 0
   !call function_W(W, j, fo, hz, Nxb, Nyb, Nz,kl,kb,z,R)
end if

if (type==2) then
 !j = -2
   call function_V(W, fo, hz, Nxb, Nyb, Nz,kl,kb,z,R,V)

   !j = -1
   !call function_V(W, j, fm1, hz, Nxb, Nyb, Nz,kl,kb,z,R,V)

   !j = 0
   !call function_V(W, j, fo, hz, Nxb, Nyb, Nz,kl,kb,z,R,V)
end if
alpha(2:Nz) = hz(2:Nz)/(hz(1:Nz-1)*(hz(2:Nz)+hz(1:Nz-1)))
alpha(1)    = alpha(2)

beta(2:Nz)   = -(hz(2:Nz)+ hz(1:Nz-1))/(hz(2:Nz)*hz(1:Nz-1))
beta(1)      = beta(2)
gamma        = -(alpha+beta)

do m = 3, Nz
   Fz(1:Nxb,1:Nyb, m)  = alpha(m)*fo(1:Nxb,1:Nyb,m-2) + beta(m)*fo(1:Nxb,1:Nyb,m-1) + gamma(m)*fo(1:Nxb,1:Nyb,m)
end do
Fz(1:Nxb,1:Nyb,2) = beta(2) * fo(1:Nxb,1:Nyb,1) + gamma(2)*fo(1:Nxb,1:Nyb,2)

!First Order Approimztion for fz(1) since W(0) = 0
Fz(1:Nxb, 1:Nyb,1)           = fo(1:Nxb,1:Nyb,1)/hz(1)


!W = W-dt*fz
!write(*,*) 'W=', W(1:Nxb,1:Nyb,1:Nz) 
!write(*,*) 'R=', R**2.
!write(*,*) 'z=', z
if (type==1) then
   do i = 1,Nz
      fpbw(1:Nxb,1:Nyb,i) = kl*W(1:Nxb,1:Nyb,i) + kb/2d0*(W(1:Nxb,1:Nyb,i)**2. + 2d0*z(i)*W(1:Nxb,1:Nyb,i)*R &
                      + z(i)**2.* R**2.) - kb/2d0*z(i)*R**2.
      !write(*,*) 'i=',i,fp(1:Nxb,1:Nyb,i), W(1:Nxb,1:Nyb,i)
   end do
end if


if (type==2) then
   do m = 3, Nz
      fpbw(1:Nxb,1:Nyb, m)  = alpha(m)*W(1:Nxb,1:Nyb,m-2) + beta(m)*W(1:Nxb,1:Nyb,m-1) + gamma(m)*W(1:Nxb,1:Nyb,m)
   end do
   fpbw(1:Nxb,1:Nyb,2)      = beta(2) * W(1:Nxb,1:Nyb,1) + gamma(2)*W(1:Nxb,1:Nyb,2)

   !First Order Approimztion for fz(1) since W(0) = 0    
   fpbw(1:Nxb, 1:Nyb,1)     = W(1:Nxb,1:Nyb,1)/hz(1)

end if
!W(1:Nxb,1:Nyb,1) = W(1:Nxb,1:Nyb,2)*((z(1)-0)*(z(1)-z(3))*(z(1)-z(4)))/((z(2)-0)*(z(2)-z(3))*(z(2)-z(4)))&

!         + W(1:Nxb,1:Nyb,3)*((z(1)-0)*(z(1)-z(2))*(z(1)-z(4)))/((z(3)-0)*(z(3)-z(2))*(z(3)-z(4))) &
!         + W(1:Nxb,1:Nyb,4) *((z(1)-0)*(z(1)-z(2))*(z(1)-z(3)))/((z(4)-0)*(z(4)-z(2))*(z(4)-z(3)))

!write(*,*) 'exit b_w'

!write(*,*) 'fz(20,20,5)', fz(20,20,1:100)

end subroutine b_w

!---------------------------------------------------------------------------
subroutine function_W(W, fji, hz, Nxb, Nyb, Nz,kl,kb,z,R)

integer:: Nxb, Nyb, Nz, m
double precision, dimension(Nxb,Nyb, 1:Nz)::W
double precision, dimension(Nxb, Nyb, 1:Nz)::fji
double precision, dimension(Nxb,Nyb):: Ws, R
double precision, dimension(1:Nz):: hz,z
double precision:: kl, kb

fji = 0

do m  = 1,Nz
   
   Ws =W(1:Nxb,1:Nyb,m)

   fji(1:Nxb,1:Nyb,m) = -1d0*(kl/2d0*Ws**2. - kb/2d0 * z(m)*R**2.*Ws &
                    + kb/2d0*(1d0/3d0*Ws**3.+z(m)*Ws**2.*R +z(m)**2.*R**2.*Ws))


end do


!write(*,*) 'fji', fji(20,3, 1:Nz)
!write(*,*) fji(10,10,1)
end subroutine function_W
!--------------------------------------------------------------------------- 
subroutine function_V(W, fji, hz, Nxb, Nyb, Nz,kl,kb,z,R,V)

integer:: Nxb, Nyb, Nz, m, n, j
double precision, dimension(Nxb,Nyb, 1:Nz)::W,V
double precision, dimension(Nxb, Nyb, 1:Nz)::fji
double precision, dimension(Nxb,Nyb):: Ws, R, Vs
double precision, dimension(1:Nz):: hz,z
double precision:: kl, kb

fji = 0
!n = 3

do m  = 1,Nz

   Ws =W(1:Nxb,1:Nyb,m)
   Vs =V(1:Nxb,1:Nyb,m)

   fji(1:Nxb,1:Nyb,m) = -1d0*((kl*Ws+kb/2d0*(Ws+z(m)*R)**2. - kb/2d0*z(m)*R**2.)*Vs)
 !  n = n+1


end do

end subroutine function_V





!--------------------------------------------------------------------------
subroutine fun_hz
integer:: I
double precision:: hzo, hmax, hend, hmin
double precision, dimension(Nz)::ztemp, hh 
write(*,*) 'Nz=', Nz

hzo = 1d0/(dble(Nz)-1d0)
!write(*,*) Nz-1d0
ztemp = (/(hzo*dble(I),I = 0,Nz-1) /)
!write(*,*) 'ztemp=', ztemp
hmin  = 5d-3
hmax   = hzo
hend  = 1d-3
hh    = 5d-1*(1d0+tanh((7d-1-ztemp)/15d-2))

!hmin + (hmax-hmin) * 0.5 * (1d0+tanh((ztemp-0.1)/0.03)) &
        !     - (hmax-hend) * 0.5 * (1d0+tanh(-(0.7-ztemp)/0.15))

hz(1:Nz)    = hh/sum(hh)
!hz(-1:0)     = hz(1)


end subroutine fun_hz
!------------------------------------------------------------------------------------------
subroutine diffusion_subr_W(Nxb, Nyb, Nz,kx, ky, cup, BCX, BCxRHS, A,B, D, BCyRHS, Diff, G, P,z,&
                           W1old, ADVw, RXN,ADVzW,Wnew,dt,jk)

integer:: Nxb, Nyb, Nz, i, j, m
double precision:: kx,ky,cup,dt
double precision, dimension(0:Nxb+1,1:Nz):: A, B, D
double precision, dimension(0:Nyb+1):: BCX
double precision, dimension(1:Nyb)::   BCxRHS, D22, RHS2y, RHSzW2
double precision, dimension(1:Nxb)::   BCyRHS, D2, RHS2, RHSzW
double precision, dimension(Nxb,Nyb)::  C
double precision, dimension(0:Nxb+1, 0:Nyb+1):: Diff, Gh,Ph, Diff2 
double precision, dimension(0:Nxb+1,0:Nyb+1,1:Nz):: G,P,Wnew
double precision, dimension(1:Nxb,1:Nyb,1:Nz):: ADVw, RXN, ADVzW
double precision, dimension(Nxb-1)::   DL, Du
double precision, dimension(Nyb-1)::   DL2, Du2
double precision, dimension(0:Nxb+1,0:Nyb+1)::zWnew, W1old, zWold
double precision, dimension(Nz):: z
integer:: halfmatinfo, matinfo,jk
Diff2 = Diff
Ph = 0d0
zWnew = 0d0
zWold = 0d0

G(1:Nxb,1:Nyb,1:Nz) = G(1:Nxb,1:Nyb,1:Nz) !+ RXN(1:Nxb,1:Nyb,1:Nz)
if (jk.le.1) then
do m = 1,Nz-1
   zWold = -z(m)*W1old
   zWnew = -z(m)*Wnew(0:Nxb+1,0:Nyb+1,Nz)
   do j = 1, Nyb
     
      call N_mat(Diff2(1:Nxb,j), kx, Nxb, DL, D2, Du, BCX(j),0d0,0)

      call RHS(transpose(Diff2(1:Nxb,j:j+1)), ky, Nxb,&
           transpose(G(1:Nxb,j-1:j+1,m)), 10, RHS2)
      !This is going to need to be rethought through with 
      call RHS_old(Diff2(1:Nxb+1,j:j+1),kx,ky,Nxb,zWold(0:Nxb+1,j-1:j+1),zWnew(0:Nxb+1,j-1:j+1),RHSzW)
           
      
      RHS2 = RHS2 + ADVw(1:Nxb,j,m) + RHSzW - dt*ADVzW(1:Nxb,j,m) + RXN(1:Nxb,j,m)
      
      call DGTSV(Nxb, 1, DL, D2, Du, RHS2, Nxb, halfmatinfo)
      Ph(1:Nxb,j) = RHS2

   end do

   call ghost_cells(Ph, A,B, D, Gh, Nxb, Nyb, 1, cup,jk)
   !write(*,*) 'Ph=', Ph(20,20,1:Nz)

   do i = 1,Nxb
      call N_mat(Diff2(i,1:Nyb), ky, Nyb, DL2, D22, Du2, A(i,m),0d0,0)

      call RHS(Diff2(i:i+1,1:Nyb), kx, Nyb, Gh(i-1:i+1,1:Nyb), i, RHS2y)

      RHS2y(1)  = RHS2y(1) + BCyRHS(i)
  !    RHS2y = RHS2y  - ADVw(i,1:Nyb,m) + RHSzW2 + RXN(i,1:Nyb,m)- dt*ADVzW(i,1:Nyb,m)  
      call DGTSV(Nyb,1,DL2,D22,Du2,RHS2y,Nyb,matinfo)
      C(i,1:Nyb) = RHS2y
   end do

   P(1:Nxb, 1:Nyb,m) = C
end do  

else 
   do m = 1, Nz-1
       zWold = -z(m)*W1old
       zWnew = -z(m)*Wnew(0:Nxb+1,0:Nyb+1,Nz)
      do j = 1, Nyb
         call N_mat(Diff(1:Nxb,j), kx, Nxb, DL, D2, Du, BCX(j),0d0,0)
         !write(*,*) 'Du=', Du
         !write(*,*) 'DL=', Dl
         !write(*,*) 'D2=', D2
!         stop
         RHS2 = (kx*Diff(1:Nxb,j)*G(0:Nxb-1,j,m) + (1d0-2d0*kx*Diff(1:Nxb,j))*G(1:Nxb,j,m) + kx*Diff(1:Nxb,j)*G(2:Nxb+1,j,m))
         !This is going to need to be rethought through with                       
         call RHS_old(Diff(1:Nxb+1,j:j+1),kx,ky,Nxb,zWold(0:Nxb+1,j-1:j+1),zWnew(0:Nxb+1,j-1:j+1),RHSzW)
         !write(*,*) maxval(RHSzW)
         
         RHS2 = RHS2 + ADVw(1:Nxb,j,m) + RHSzW + RXN(1:Nxb,j,m) - dt*ADVzW(1:Nxb,j,m)
         call DGTSV(Nxb, 1, DL, D2, Du, RHS2, Nxb, halfmatinfo)
         C(1:Nxb, j ) = RHS2
      end do
      P(1:Nxb,1:Nyb,m) = C(1:Nxb,1:Nyb)
   end do
endif
P(1:Nxb,1:Nyb,Nz) = Wnew(1:Nxb,1:Nyb,Nz)
!write(*,*) 'P=', P(20,20,1:Nz)
call ghost_cells(P, A,B,D, G, Nxb, Nyb, Nz,Cup,jk)

end subroutine diffusion_subr_W
!------------------------------------------------------------------------------------------
subroutine diffusion_subr_V(Nxb, Nyb, Nz,kx, ky, cup, BCX, BCxRHS, A,B, D, BCyRHS, Diff, G, P, &
                           ADVV, RXNV,ADVzV,dt,jk)

integer:: Nxb, Nyb, Nz, i, j, m
double precision:: kx,ky,cup,dt 
double precision, dimension(0:Nxb+1,1:Nz):: A, B, D
double precision, dimension(0:Nyb+1):: BCX
double precision, dimension(1:Nyb)::   BCxRHS, D22, RHS2y, RHSzW2
double precision, dimension(1:Nxb)::   BCyRHS, D2, RHS2, RHSzW
double precision, dimension(Nxb,Nyb)::  C
double precision, dimension(0:Nxb+1, 0:Nyb+1):: Diff, Gh,Ph 
double precision, dimension(0:Nxb+1,0:Nyb+1,1:Nz):: G,P
double precision, dimension(1:Nxb,1:Nyb,1:Nz):: ADVV, RXNV, ADVzV
double precision, dimension(Nxb-1)::   DL, Du
double precision, dimension(Nyb-1)::   DL2, Du2
integer:: halfmatinfo, matinfo,jk

Ph = 0d0
!write(*,*) 'RXN=', RXNV(1,1,Nz-4:Nz), 'ADVzV', ADVzV(1,1,Nz-4:Nz)
if (jk.le.1) then
do m = 1,Nz
   do j = 1, Nyb
     
      call N_mat(Diff(1:Nxb,j), kx, Nxb, DL, D2, Du, BCX(j),0d0,0)

      call RHS(transpose(Diff(1:Nxb,j:j+1)), ky, Nxb,&
           transpose(G(1:Nxb,j-1:j+1,m)), 10, RHS2)
      !This is going to need to be rethought through with 

      RHS2 = RHS2 + ADVV(1:Nxb,j,m) + RXNV(1:Nxb,j,m)- dt*ADVzV(1:Nxb,j,m)
      
      call DGTSV(Nxb, 1, DL, D2, Du, RHS2, Nxb, halfmatinfo)
      Ph(1:Nxb,j) = RHS2

   end do

   call ghost_cells(Ph, A,B, D, Gh, Nxb, Nyb, 1, cup,jk)
   !write(*,*) 'Ph=', Ph(20,20,1:Nz)

   do i = 1,Nxb
      call N_mat(Diff(i,1:Nyb), ky, Nyb, DL2, D22, Du2, A(i,m),0d0,0)

      call RHS(Diff(i:i+1,1:Nyb), kx, Nyb, Gh(i-1:i+1,1:Nyb), i, RHS2y)
   
      RHS2y(1)  = RHS2y(1) + BCyRHS(i)
  !    RHS2y = RHS2y  - ADVw(i,1:Nyb,m) + RHSzW2 + RXN(i,1:Nyb,m)- dt*ADVzW(i,1:Nyb,m)  
      call DGTSV(Nyb,1,DL2,D22,Du2,RHS2y,Nyb,matinfo)
      C(i,1:Nyb) = RHS2y
   end do

   P(1:Nxb, 1:Nyb,m) = C
end do  
else
   do m = 1, Nz
      do j = 1, Nyb
         call N_mat(Diff(1:Nxb,j), kx, Nxb, DL, D2, Du, BCX(j),0d0,0)
         RHS2 = (kx*Diff(1:Nxb,1)*G(0:Nxb-1,j,m) + (1d0-2d0*kx*Diff(1:Nxb,1))*G(1:Nxb,j,m) + kx*Diff(1:Nxb,j)*G(2:Nxb+1,j,m))
         !This is going to need to be rethought through with                                    
          RHS2 = RHS2 + ADVV(1:Nxb,j,m) + RXNV(1:Nxb,j,m)- dt*ADVzV(1:Nxb,j,m)
         call DGTSV(Nxb, 1, DL, D2, Du, RHS2, Nxb, halfmatinfo)
         C(1:Nxb, j ) = RHS2
      end do
      P(1:Nxb,1:Nyb,m) = C(1:Nxb,1:Nyb)
   end do
endif

!P(1:Nxb,1:Nyb,Nz) = Wnew(1:Nxb,1:Nyb,Nz)
!write(*,*) 'P=', P(20,20,1:Nz)
call ghost_cells(P, A,B,D, G, Nxb, Nyb, Nz,Cup,jk)

end subroutine diffusion_subr_V
!------------------------------------------------------------------------------------------
subroutine diffusion_subr_r(Nxb, Nyb, kx, ky, cup, BCX, BCxRHS, A, D, B, BCyRHS, Diff, Rg, G, Rgold, RXNr, ADVr,dt,jk)
integer:: Nxb, Nyb, i, j,jk
double precision:: kx,ky,cup,dt
double precision, dimension(0:Nxb+1):: A,D,B
double precision, dimension(0:Nyb+1):: BCX
double precision, dimension(1:Nyb)::   BCxRHS, D22, RHS2y,RHSg2
double precision, dimension(1:Nxb)::   BCyRHS, D2, RHS2, RHSg
double precision, dimension(Nxb,Nyb)::  ADVr 
double precision, dimension(0:Nxb+1, 0:Nyb+1):: Rg, Rgold, G, Diff2,Diff, Gh, P, Ph, RXNr 
double precision, dimension(Nxb-1)::   DL, Du
double precision, dimension(Nyb-1)::   DL2, Du2
integer:: halfmatinfo, matinfo
P=0d0
Diff2 = Diff
G(1:Nxb,1:Nyb) = G(1:Nxb,1:Nyb)!+ADVr(1:Nxb,1:Nyb) + RXNr(1:Nxb,1:Nyb)

!do j = 1, Nyb
!   call RHS_old(Diff2(1:Nxb+1,j:j+1),kx,ky,Nxb,-Rg(0:Nxb+1,j-1:j+1),-Rgold(0:Nxb+1,j-1:j+1),RHSg)
!   G(1:Nxb,j) = G(1:Nxb,j) + RHSg
!end do

if (jk.le.1) then

do j = 1, Nyb
     call N_mat(Diff2(1:Nxb,j), kx, Nxb, DL, D2, Du, BCX(j),0d0,1)

     call RHS(transpose(Diff(1:Nxb,j:j+1)), ky, Nxb,&
          transpose(G(1:Nxb,j-1:j+1)), 10, RHS2)

     ! Change the second Rgold to Rg if you want to use a time average
     call RHS_old(Diff2(1:Nxb+1,j:j+1),kx,ky,Nxb,-Rg(0:Nxb+1,j-1:j+1),-Rg(0:Nxb+1,j-1:j+1),RHSg)

!     write(*,*) 'RHSG=', RHSg
    RHS2(1)     = RHS2(1)+BCxRHS(j)
  !write(*,*) 'RHSold=', RHSg
     RHS2        = RHS2 + RHSg +ADVr(1:Nxb,j) + RXNr(1:Nxb,j)

!     write(*,*) 'j=', j, 'RHS2=', RHS2, 'RHSg=', RHSg, 'RXNr=', RXNr(1:Nxb,j)
     call DGTSV(Nxb, 1, DL, D2, Du, RHS2, Nxb, halfmatinfo)
!     write(*,*) 'halfmatinfo=', halfmatinfo
     Ph(1:Nxb,j) = RHS2
!     write(*,*) 'D2=', D2
 !    write(*,*) 'P=', Ph(1:Nxb,j)

  end do

  call ghost_cells(Ph, A, B, D, Gh, Nxb, Nyb, 1,cup,jk)

  do i = 1,Nxb
     call N_mat(Diff2(i,1:Nyb), ky, Nyb, DL2, D22, Du2, A(i),D(i),0)

     call RHS(Diff(i:i+1,1:Nyb), kx, Nyb, Gh(i-1:i+1,1:Nyb), i, RHS2y)
     RHS2y(1)  = RHS2y(1) + BCyRHS(i)

     call DGTSV(Nyb,1,DL2,D22,Du2,RHS2y,Nyb,matinfo)
     P(i,1:Nyb) = RHS2y
  end do
else
   do j = 1,Nyb
         call N_mat(Diff(1:Nxb,j), kx, Nxb, DL, D2, DU, BCX(j),0d0,0)
         RHS2 = (kx*Diff(1:Nxb,1)*G(0:Nxb-1,j) + (1d0-2d0*kx*Diff(1:Nxb,1))*G(1:Nxb,j) + kx*Diff(1:Nxb,j)*G(2:Nxb+1,j))
        
        !This is going to need to be rethought through with   
         call RHS_old(Diff(1:Nxb+1,j:j+1),kx,ky,Nxb,-Rg(0:Nxb+1,j-1:j+1),-Rgold(0:Nxb+1,j-1:j+1),RHSg)
        !if (j==1) then
        !   write(*,*)'RHS=', RHS2
        !end if
          RHS2 = RHS2 + ADVr(1:Nxb,j) + RHSg + RXNr(1:Nxb,j)
         call DGTSV(Nxb, 1, DL, D2, Du, RHS2, Nxb, halfmatinfo)
         P(1:Nxb, j) = RHS2
      end do
endif


!  write(*,*) 'p=', P(1:Nxb,1)
  call ghost_cells(P, A, B, D, G, Nxb, Nyb, 1, Cup,jk)
!  write(*,*) 'Rnew=', G(1:Nxb,1)
end subroutine diffusion_subr_R


!--------------------------------------------<---------------------- 

subroutine ghost_cells(C,A,B, DC, Cnew,Nxb,Nyb,Nz,Cup,jk)

integer:: Nxb, Nyb, Nz,jk
double precision:: Cup
double precision, dimension(0:Nxb+1,0:Nyb+1,1:Nz)::C, Cnew
double precision, dimension(0:Nxb+1,1:Nz):: A, B, DC

Cnew                     = C

Cnew(0:Nxb+1,0,1:Nz)      = A(0:Nxb+1,1:Nz)*C(0:Nxb+1,1,1:Nz)+ DC(0:Nxb+1,1:Nz)*C(0:Nxb+1,2,1:Nz)&
                           + B(0:Nxb+1,1:Nz)
Cnew(0:Nxb+1,Nyb+1,1:Nz)   = C(0:Nxb+1,Nyb,1:Nz)
if (jk.lt.0) then
   Cnew(0,0:Nyb+1,1:Nz)      = 2d0*Cup-C(1,0:Nyb+1,1:Nz) 
else
   Cnew(0,0:Nyb+1,1:Nz)      =  C(1,0:Nyb+1,1:Nz)
endif
Cnew(Nxb+1,0:Nyb+1,1:Nz)   = C(Nxb,0:Nyb+1,1:Nz)
!Cnew(0,0,1:Nz)         = Cnew(1,0,1:Nz)
!Cnew(Nxb+1,0,1:Nz)      = Cnew(Nxb,0,1:Nz)
!Cnew(0,Nyb+1,1:Nz)      = Cnew(0,Nyb,1:Nz)
!Cnew(Nxb+1,Nyb+1,1:Nz)   = Cnew(Nxb+1,Nyb,1:Nz)
!if (maxval(B).gt.0d0) then
!   write(*,*) 'B=', B
!   write(*,*) 'Cnew=', Cnew(0:Nxb+1,0,1)
!endif


end subroutine ghost_cells


!-----------------------------------------------------------------
subroutine N_mat(Diff,k,N,Dl,Diag,Du,BC, BC2,writepar)

implicit none

integer:: N, writepar
double precision::k, BC, BC2
double precision, dimension(N):: Diag
double precision, dimension(N-1):: Dl, Du
double precision, dimension(N):: Diff

!Diagonal    
!write(*,*) 1d+00-BC      
Diag(1)      = 1d+00 + k*(Diff(2)+Diff(1)*(1d0-BC))
Diag(2:N-1)  = 1d+00 + k*(Diff(2:N-1)+Diff(3:N))
Diag(N)      = 1d+00 + k*Diff(N)

!if (writepar==1) then
!write(*,*) 'Diag=', Diag
!end if
!Lower Diagonal 

Dl = -k*Diff(2:N)


!Upper Diagonal 

Du    = -k*Diff(2:N)
Du(1) = -k*(Diff(2)+Diff(1)*BC2)


end subroutine N_mat
!------------------------------------------------------------------------------------------
subroutine RHS(D,k,N,C,i, RHSvec)

integer:: N, i, m
double precision,dimension(2,N):: D
double precision, dimension(3,N):: C
double precision:: k
double precision, dimension(N)::RHSvec

RHSvec = C(2,1:N) + k*(D(2,1:N)*(C(3,1:N)-C(2,1:N)) - D(1,1:N)*(C(2,1:N)-C(1,1:N))) 
!if (i==1) then
!   write(*,*) 'C(3:N)=', C
!   write(*,*) 'k=', k
!   write(*,*) 'D=', D
!end if
end subroutine RHS
!-----------------------------------------------------------------------------------------
subroutine diffusion_subr(Nxb, Nyb, kx, ky, cup, BCX, BCxRHS, A, D, B, BCyRHS, Diff, G, C,ADV,RK)
integer:: Nxb, Nyb, i, j, jk
double precision:: kx,ky,cup
double precision, dimension(0:Nxb+1):: A,D,B
double precision, dimension(0:Nyb+1):: BCX
double precision, dimension(1:Nyb)::   BCxRHS, D22, RHS2y
double precision, dimension(1:Nxb)::   BCyRHS, D2, RHS2
double precision, dimension(Nxb,Nyb)::  C, ADV,Rk
double precision, dimension(0:Nxb+1, 0:Nyb+1):: G, Diff, Gh, P, Ph 
double precision, dimension(Nxb-1)::   DL, Du
double precision, dimension(Nyb-1)::   DL2, Du2
integer:: halfmatinfo, matinfo
!-------------------


!-----------------

!write(*,*) 'G=', G(1:Nxb,1)
P=G
Ph = 0d0
!write(*,*) 'G=', G(35,1:Nyb)
!write(*,*) 'BCx=', BCx
!write(*,*) 'BCxRHS=', BCxRHS
!write(*,*) 'BCyRHS=', BCyRHS
!write(*,*) '--------'
do j = 1, Nyb
     call N_mat(Diff(1:Nxb,j), kx, Nxb, DL, D2, Du, BCX(j),0d0,1)
     !if (j==1) then
     !   write(*,*) 'D2=', D2
     !   write(*,*) 'TDiff=', transpose(Diff(1:Nxb,j:j+1))
     !   write(*,*) 'ky=', ky
     !   write(*,*) 'Nxb=', Nxb
     !   write(*,*) 'TG=', transpose(G(1:Nxb,j-1:j+1))
     !end if
     call RHS(transpose(Diff(1:Nxb,j:j+1)), ky, Nxb,&
          transpose(G(1:Nxb,j-1:j+1)), j, RHS2)
 !    if (maxval(abs(RHS2))/Cup > 1d0) then
 !    write(*,*) 'j=', j, maxval(abs(RHS2))
 !    write(*,*) 'G(1:Nxb,j-1,j+1)=', G(1:Nxb,j-1:j+1)
 !    stop
 !    end if
     !  if (j==1) then
      !  write(*,*) 'RHS2 =', RHS2
      !  end if
     RHS2(1) = RHS2(1)+BCxRHS(j)
     RHS2    = RHS2+ADV(1:Nxb,j)+RK(1:Nxb,j) 
     !if (j ==Nyb) then
        !write(*,*) 'BCxRHS(j)=', BCxRHS(35)
        !write(*,*) 'ADV(1:Nxb,j)=', ADV(35,j)
        !write(*,*) 'RK(1:Nxb,j)=', RK(35,j)
        !write(*,*) 'RHS=', RHS2(1:Nxb)
        !write(*,*) 'Gx=', G(35,j-1:j+1)
        !write(*,*) 'DL = ', Dl, 'D2=', D2,'Du=', Du
     !endif
     call DGTSV(Nxb, 1, DL, D2, Du, RHS2, Nxb, halfmatinfo)
     !write(*,*) 'halfmatinfo=', halfmatinfo
     !write(*,*) 'RHS2=', RHS2
     Ph(1:Nxb,j) = RHS2
     
  end do
  !write(*,*) 'BCx=', BCxRHS
  !write(*,*) 'Ph(1:nx,1)=', Ph(35,1:Nyb)

  call ghost_cells(Ph, A, B, D, Gh, Nxb, Nyb, 1,cup,jk)
  !write(*,*) 'Gh(1:Nxb,1)=', Gh(35,1:Nyb)
  do i = 1,Nxb
     call N_mat(Diff(i,1:Nyb), ky, Nyb, DL2, D22, Du2, A(i),D(i),1)

    call RHS(Diff(i:i+1,1:Nyb), kx, Nyb, Gh(i-1:i+1,1:Nyb), 10, RHS2y)
     RHS2y(1)  = RHS2y(1) + BCyRHS(i)
     call DGTSV(Nyb,1,DL2,D22,Du2,RHS2y,Nyb,matinfo)
     C(i,1:Nyb) = RHS2y
  end do

  P(1:Nxb, 1:Nyb) = C
!  write(*,*) 'P', P(1:Nxb,1)
!write(*,*) 'P=', P(35:66,1:Nyb)   
call ghost_cells(P, A, B, D, G, Nxb, Nyb,1, Cup,jk)
!write(*,*) 'G=', G(35:66,0:Nyb)
end subroutine diffusion_subr

!------------------------------------------------------------------
function philim(theta, casenumber)

integer:: casenumber
double precision:: theta, philim 

!select case (casenumber)

!case (1)
   !minmod                       
!   philim=max(0.0,min(1.0,theta))
!   return
!case (2)
   !superbee                               
!   philim = max(0.0,min(1.0,2.0*theta), min(2.0,theta))
!   return
!case (3)
   !van Leer                          
!   philim = (theta + abs(theta))/(1.0+abs(theta))
!   return
!case (4)
   !monotonized center
   philim = max(0.0, min((1.0+theta)/2.0,2.0,2.0*theta))
!   return
!end select

end function philim


!----------------------------------------------------------------------------
subroutine advection(C, Ulr, Utb, Vlr,Vtb, dt, hx, hy, Nxb, Nyb,ADV)
integer:: Nxb, Nyb
double precision, dimension(-1:Nxb+2, -1:Nyb+2):: F, G
double precision, dimension(-1:Nxb+2, -1:Nyb+2):: C, Cnew, Ulr, Utb, Vlr, Vtb
double precision:: hx, hy,dt, Uij, Vij, Rcij, S, theta, phi
integer:: i,j, casenumber, I2, I3, J2, J3
double precision, dimension(1:Nxb,1:Nyb):: ADV
double precision, dimension(-1:Nyb+2):: Rcvec
!write(*,*) C                     
F = 0.0
G = 0.0
Cnew = 0.0
do j = 0, Nyb+1
   do i = 1, Nxb+1
      Uij  = Ulr(i,j)
      Vij  = Vlr(i,j)
      Rcij = C(i,j) - C(i-1,j)

      if (Uij > 0d0) then
         I2 = i-1
      else
         I2 = i
      end if

      F(i,j) = F(i,j) + Uij*C(I2,j)

      Rcij = C(i,j) - C(i-1,j)
      Uij  = Ulr(i,j)
      Vij  = Vlr(i,j)
      !-----------------------------      
      if (Uij > 0d0) then
         I2 = i
      else
         I2 = i-1
      end if
      if (Vij > 0d0) then
         J2 = j+1
      else
         J2 = j
      end if

      G(I2,J2) = G(I2,J2) - 0.5*dt/hy*Uij*Vij*Rcij

!      write(*,*) I2,J2, G(I2,J2)                                           
      !---------------------------------------       
      !Limited R:   
      if (abs(RCij)==0.0) then
         theta = 0.0
      else
         if (Uij > 0.0) then
            I3 = i-1
         else
            I3 = i+1
         end if

         theta = (C(I3,j) - C(I3-1,j))/RCij


      end if
        casenumber = 4
        phi=philim(theta, casenumber)


      RCij  = RCij*phi
      S     = 0.5*abs(Uij)*(1.0-dt/hx*abs(Uij))*RCij
      F(i,j)= F(i,j) + S
!      if (j==0.or.j==1) then                
!         write(*,*) i, j, F(i,j)               
!      end if                               
      !--------------------------------------- 
      G(i,J2)   = G(i,J2)   + dt/hy*Vij*S
      G(i-1,J2) = G(i-1,J2) - dt/hy*Vij*S
!      write(*,*) i,J2,G(i,J2), G(i-1,J2)           
   end do
end do


do i = 0, Nxb+1
   !Rcvec = C(i,0:Nyb+2) - C(i,-1:Nyb+2)     
   do j = 1, Nyb+1
      Uij  = Utb(i,j)
      Vij  = Vtb(i,j)
!      write(*,*) i,j, Vij, Uij       
      Rcij = C(i,j) - C(i,j-1)
      if (Vij > 0d0) then
         J2 = j-1
      else
         J2 = j
      end if
      G(i,j) = G(i,j) +Vij*C(i,J2)

      !-----------------------------  
      if (Vij > 0d0) then
         J2 = j
      else
         J2 = j-1
      end if
      if (Uij > 0d0) then
         I2 = i+1
      else
         I2 = i
      end if
!      if (I2==1) then                             
!      write(*,*) I2,J2,F(I2,J2), RCij, Vij, Uij      
!      end if                              
      F(I2,J2) = F(I2,J2) - 5d-1*dt/hx*Vij*Uij*Rcij
!      if (I2==1) then                     
!      write(*,*) I2, J2, F(I2,J2)            
!      end if                          
  end do
do j = 0, Nyb+1
      Rcij = C(i,j) - C(i,j-1)
      Uij  = Utb(i,j)
      Vij  = Vtb(i,j)
      if (Uij > 0d0) then
         I2 = i+1
      else
         I2 = i
      end if
      if (Vij > 0d0) then
         J2 = j
      else
         J2 = j-1
      end if
      !--------------------------------------- 
      !Limited                       
      if (abs(RCij)==0.0) then
         theta = 0.0
      else
         if (Vij > 0.0) then
            J3 = j-1
         else
            J3 = j+1
        end if

         theta = (C(i,J3) - C(i,J3-1))/RCij
     end if

      casenumber = 4

      phi=philim(theta, casenumber)
      RCij   = RCij*phi
!      write(*,*) i, RCij     
      S      = 5d-1*abs(Vij)*(1d0-dt/hy*abs(Vij))*RCij

      G(i,j) = G(i,j) + S
!      write(*,*) I2, j, F(I2,j)   
!if (j==0.or.j==1) then       
!         write(*,*) I2,j, F(I2,j), F(I2,j-1)   
!      end if                   
      !---------------------------------------  
      F(I2,j)   = F(I2,j)   + dt/hx*Utb(I2,j)*S
      F(I2,j-1) = F(I2,j-1) - dt/hx*Utb(I2,j-1)*S
!      if (j==0.or.j==1) then  
!         write(*,*) I2,j, F(I2,j), F(I2,j-1) 
!      end if                 
   end do
end do
do i = 1, Nxb
   do j = 1,Nyb
      ADV(i,j) =  - dt/hx*(F(i+1,j)-F(i,j))-dt/hy*(G(i,j+1)-G(i,j))

   end do
end do
!write(*,*) C(1:Nxb,1:Nyb)      
end subroutine advection

!-------------------------------------------------------------------------


subroutine ghost_cell_vec(C,Cnew,A,B,DC,Nxb,Nyb,Nz,Cup,jk)
integer:: Nxb, Nyb, Nz, jk
double precision, dimension(0:Nxb+1,0:Nyb+1,1:Nz) :: C
double precision, dimension(-1:Nxb+2, -1:Nyb+2,1:Nz):: Cnew
double precision:: Cup
double precision, dimension(0:Nxb+1,1:Nz):: A, B, DC

Cnew(1:Nxb, 1:Nyb,1:Nz)    = C(1:Nxb,1:Nyb,1:Nz)
Cnew(0:Nxb+1,0,1:Nz)      = A(0:Nxb+1,1:Nz)*C(0:Nxb+1,1,1:Nz)+ &
                           DC(0:Nxb+1,1:Nz)*C(0:Nxb+1,2,1:Nz) + B(0:Nxb+1,1:Nz)
Cnew(0:Nxb+1,Nyb+1,1:Nz)   = C(0:Nxb+1,Nyb,1:Nz)
if (JK.lt.0) then
   Cnew(0,0:Nyb+1,1:Nz)      = 2d0*Cup-C(1,0:Nyb+1,1:Nz) 
else
   Cnew(0,0:Nyb+1,1:Nz)      = C(1,0:Nyb+1,1:Nz)
endif 
Cnew(Nxb+1,0:Nyb+1,1:Nz)   = C(Nxb,0:Nyb+1,1:Nz)

Cnew(-1:Nxb+2,-1,1:Nz)    = 0d0
Cnew(-1:Nxb+2,Nyb+2,1:Nz)  = 0d0

Cnew(0,-1,1:Nz)   = Cup
Cnew(0,Nyb+2,1:Nz) = Cup

Cnew(-1,-1:Nyb+2,1:Nz)    = Cup

Cnew(Nxb+2, 0:Nyb+1,1:Nz) = C(Nxb,0:Nyb+1,1:Nz)


Cnew(-1,0,1:Nz)=Cup
Cnew(Nxb+2,0,1:Nz) = Cup!0d0       
Cnew(-1,Nyb+1,1:Nz) = 0d0
Cnew(Nxb+2,Nyb+1,1:Nz) = 0d0

Cnew(Nxb+1,-1,1:Nz) =  0d0
Cnew(Nxb+1,Nyb+2,1:Nz)= 0d0

!write(*,*) 'Cnew=', Cnew                                  

end subroutine ghost_cell_vec

!---------------------------------------------------------------------- 
subroutine RungaKutta(E,C,kappaf,kappafs,Nxb,Nyb,dt,Rkf)
integer:: Nxb, Nyb
double precision:: dt, kappaf,kappafs
double precision, dimension(0:Nxb+1,0:Nyb+1):: E, C
double precision, dimension(1:Nxb,1:Nyb):: Rkf,S1
!write(*,*) 5d-0*dt*kf*E(1:33,1:10)*C(1:33,1:10)/(Kfs+C(1:33,1:10))

S1 = C(1:Nxb,1:Nyb)-(5d-1*dt*kappaf*E(1:Nxb,1:Nyb)*C(1:Nxb,1:Nyb))/(kappafs+C(1:Nxb,1:Nyb))
!Cnew  = C-(dt*kf*E*S1)/(kfs+S1)
Rkf = -dt*kappaf*E(1:Nxb,1:Nyb)*S1(1:Nxb,1:Nyb)/(kappafs+S1(1:Nxb,1:Nyb))
C(1:Nxb,1:Nyb) = C(1:Nxb,1:Nyb)+ Rkf
end subroutine RungaKutta
!--------------------------------------------------------------------------------------
subroutine RXN_R(Nxb, Nyb, R, Rg, Rs, S10s, kl, kb, chi,dt, rxnr)

integer:: Nxb, Nyb
double precision:: kb, kl, chi, dt
double precision, dimension(0:Nxb+1,0:Nyb+1)::R, Rg, Rs, Rgavg, rxnr, s10s, mid



rxnr = dt*(-kl*(R**2. - chi*Rg**2.) - kb/2d0*(R**3.- chi*(3d0*Rs*Rg**2.+Rg**3.)) + 2d0*S10s)
!write(*,*) 'R=', R(1:Nxb,1)
!write(*,*) 'Rg=', Rg(1:Nxb,1)
!write(*,*) 'Rs=', Rs(1:Nxb,1)
!write(*,*) 'S10=', S10(1:Nxb,1)
!write(*,*) 'rxnr=', rxnr(1:Nxb,1)
!write(*,*) 'S10dt=', dt*2d0*S10(1:Nxb,1)      
!mid=-kl*(R**2-chi*Rgavg**2)-kb/2d0*(R**3-chi*(3*Rs*Rgavg**2+Rgavg**3))
!write(*,*)'rxnr(1:Nxb,1)=',  rxnr(1:Nxb,1)
!write(*,*) 's10=', 2d0*S10(1:Nxb,20)                                       
!write(*,*) 'dt=', dt                                     

end subroutine RXN_R
!--------------------------------------------------------------------------------------    
subroutine RXNW(Nxb,Nyb,Nz,chi,Rg, R,Rs,z, kl, kb, dt,Rxn_zW)
integer:: Nxb, Nyb,Nz, m
double precision, dimension(1:Nxb,1:Nyb):: R, Rg,Rs
double precision:: kb, kl, chi,dt
double precision, dimension(Nz):: z
double precision, dimension(1:Nxb,1:Nyb,1:Nz):: Rxn_zW
!write(*,*) 'first in RXNW=', Rs(1:Nxb,1)
!write(*,*) 'z=', z, 'kb=', kb, 'kl=', kl, 'chi=', chi
do m = 1, Nz

Rxn_zW(1:Nxb,1:Nyb,m) =dt*(-kl*chi*z(m)*Rg**2. + kb/2d0*z(m)**2.*R**3.) - &
                    dt*(kb/2d0*(z(m)*R**3.+3d0*z(m)*chi*(Rs)*Rg**2.+z(m)*chi*Rg**3.))

end do
!write(*,*) 'RXNW: R(1:Nxb,1)=', Rxn_zW(1:Nxb,1,50)
end subroutine RXNW
!----------------------------------------------------------------------
subroutine diffusion_subr_tb(Nxb,Nyb, Cold, C, Cs, RXN, kx, ky, Diff, ADV, A,B,D,Cup,jk)

integer:: Nxb, Nyb, i ,j,jk
double precision:: kx, ky,  Cup
double precision, dimension(0:Nxb+1,0:Nyb+1):: Cold, C, Cs,G, Diff
double precision, dimension(0:Nxb+1,0:Nyb+1):: CsAr, CsARg
double precision, dimension(1:Nxb,1:Nyb):: RXN, ADV
double precision, dimension(1:Nxb,1:Nyb):: RHS2
double precision, dimension(0:Nxb+1):: A,B,D
C=Cold    

do j = 1,Nyb
   do i = 1, Nxb
      RHS2(i,j) = 2d0*ky*(Diff(i,j+1)*(Cs(i,j+1)-Cs(i,j)) - Diff(i,j)*(Cs(i,j)-Cs(i,j-1))) &
       +2d0*kx*(Diff(i+1,j)*(Cs(i+1,j)-Cs(i,j))-Diff(i,j)*(Cs(i,j)-Cs(i-1,j)))     
      !RHS2(i,j) = 2d0*kx*(Diff(i+1,j)*(Cs(i+1,j)-Cs(i,j))-Diff(i,j)*(Cs(i,j)-Cs(i-1,j)))
   end do
end do
!write(*,*) 'RHS2(max,min) =', maxval(RHS2), minval(RHS2)
!write(*,*) maxval(RHS2), minval(RHS2)
C(1:Nxb,1:Nyb) = Cold(1:Nxb,1:Nyb) + RHS2 + RXN + ADV

call ghost_cells(C, A, B, D, G, Nxb, Nyb, 1, Cup,jk)


C=G

end subroutine diffusion_subr_tb


!------------------------------------------------------------------------------------------ 
subroutine RHS_old(D,kx,ky,Nxb,cold,cnew,RHSC)
integer:: Nxb
double precision, dimension(0:Nxb+1,3):: cold, cnew, C
double precision, dimension(Nxb):: RHSC
double precision, dimension(Nxb+1,2):: D
double precision:: kx, ky


C    = (cold+cnew)/2d0

!Ky is already divided by 2 so multiply it by 2                      
RHSC = 2d0*ky*(D(1:Nxb,2)*(C(1:Nxb,3)-C(1:Nxb,2)) - D(1:Nxb,1)*(C(1:Nxb,2)-C(1:Nxb,1))) &
       +2d0*kx*(D(2:Nxb+1,2)*(C(2:Nxb+1,2)-C(1:Nxb,2))-D(1:Nxb,2)*(C(1:Nxb,2)-C(0:Nxb-1,2)))
!RHSC = 2d0*kx*D(1,1) * (C(2:Nxb+1,2) - 2d0*C(1:Nxb,2) + C(0:Nxb-1,2))
!if ( maxval(abs(2d0*ky*(D(1:Nxb,2)*(C(1:Nxb,3)-C(1:Nxb,2)) - D(1:Nxb,1)*(C(1:Nxb,2)-C(1:Nxb,1))))) >1d-14) then
!write(*,*) 'Error in RHSold, j=', j                                                
!end if                                                           

end subroutine RHS_old
!--------------------------------------------------------------------------------------  

subroutine RXN_V(W,Wz,z,R,S10s, kb, Nxb, Nyb, Nz,dt, RXNv)
integer:: Nxb, Nyb, Nz, m
double precision:: kb, dt
double precision, dimension(Nz):: z
double precision, dimension(1:Nxb,1:Nyb,1:Nz):: W, Wz, RXNv
double precision, dimension(1:Nxb,1:Nyb):: R, S10s

do m = 1, Nz
   RXNv(1:Nxb,1:Nyb,m) =dt*(-kb*(W(1:Nxb,1:Nyb,m)+z(m)*R)**2.*(Wz(1:Nxb,1:Nyb,m)+R)) +dt*2d0*S10s*z(m)
end do


end subroutine RXN_V
!-----------------------------------------------------
!subroutine LaxWendroff(W,Nxb,Nyb,Nz,kl,kb,dt,R,hz,z,LW_z)
!integer:: k, Nxb, Nyb, Nz
!double precision:: dt, kl,kb
!double precision, dimension(Nz):: hz,z
!double precision, dimension(Nxb,Nyb,Nz):: W, Wphalf, Wmhalf, LW_z, fji, fphalf, fmhalf
!double precision, dimension(Nxb,Nyb):: R
!call function_W(W, fji, hz, Nxb, Nyb, Nz,kl,kb,z,R)
!write(*,*) fji
!Wphalf(1:Nxb,1:Nyb,1) = 1d0/2d0*(W(1:Nxb,1:Nyb,2) + W(1:Nxb,1:Nyb,1))   - dt/(2d0*hz(1))*(fji(1:Nxb,1:Nyb,2) - fji(1:Nxb,1:Nyb,1))
!Wmhalf(1:Nxb,1:Nyb,1) = 1d0/2d0*(W(1:Nxb,1:Nyb,1))   - dt/(2d0*hz(1))*(fji(1:Nxb,1:Nyb,1))

!do k = 2,Nz-1
!   Wphalf(1:Nxb,1:Nyb,k) = 1d0/2d0*(W(1:Nxb,1:Nyb,k+1) + W(1:Nxb,1:Nyb,k))   - dt/(2d0*hz(k))*(fji(1:Nxb,1:Nyb,k+1) - fji(1:Nxb,1:Nyb,k))
!   Wmhalf(1:Nxb,1:Nyb,k) = 1d0/2d0*(W(1:Nxb,1:Nyb,k)   + W(1:Nxb,1:Nyb,k-1)) - dt/(2d0*hz(k))*(fji(1:Nxb,1:Nyb,k)   - fji(1:Nxb,1:Nyb,k-1))
!end do

!Wmhalf(1:Nxb,1:Nyb,Nz) = 1d0/2d0*(W(1:Nxb,1:Nyb,Nz)   + W(1:Nxb,1:Nyb,Nz-1)) - dt/(2d0*hz(Nz))*(fji(1:Nxb,1:Nyb,Nz)   - fji(1:Nxb,1:Nyb,Nz-1))
!Wphalf(1:Nxb,1:Nyb,Nz) = Wmhalf(1:Nxb,1:Nyb,Nz) !This is completely wrong, but we dont want to use LW at z=1

!call function_W(Wphalf, fphalf, hz, Nxb, Nyb, Nz,kl,kb,z,R)
!call function_W(Wmhalf, fmhalf, hz, Nxb, Nyb, Nz,kl,kb,z,R)

!write(*,*) 'fph=', fphalf(1:Nxb,1,1)
!write(*,*) 'fmh=', fmhalf(1:Nxb,1,1)

!do k = 1,Nz
!   LW_z(1:Nxb,1:Nyb,k) = 1d0/hz(k)*(fphalf(1:Nxb,1:Nyb,k) - fmhalf(1:Nxb,1:Nyb,k))
!end do
!!write(*,*) LW_z
!end subroutine LaxWendroff

end module gel_mod
!------------------------------------------------
program gelation
  use gel_mod
  use param_mod
  use clot_mod

 
  !integer, parameter:: Nxb =66, Nyb =33, Nz =100   
  call init
  write(*,*) 'calling driver'
!  write(*,*) 'Nx=', Nxb
  call driver

end program gelation
