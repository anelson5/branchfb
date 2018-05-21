module gel_mod
   use grid_mod
   use param_mod
   use bobconst_mod
   use diffusion_mod
   use fluidsolve_mod

   implicit none

   integer,parameter:: Nz = 100
   double precision, dimension(1:Nz):: hz
   double precision:: Nxbd, Nybd, Nzd
   double precision:: Lx, Ly,Lz, hx, hy             !domain parameters
   double precision:: to, tend, t, dt1, dto, timerecord, muc,pic, nu
   double precision:: kx, ky, dt_fluid
   double precision:: DD, dalpha, dbeta, dtc
   double precision:: Kappaon, Kappaoff, kappaoc
   integer::i, j, N1, N2, ix, initialdata, iy, indexrecord, recordidx, dogelation
   integer:: dosource,n, jz, flowon,jk, idx_Bw, bigdomain
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
   double precision, dimension(0:Nxb+1,0:Nyb+1):: Fg, Fpc, Fgh, Fph
   double precision, dimension(-1:Nxb+2,-1:Nyb+2):: Fpv
   !--------------W parameters------------------------------
   integer:: m, nt
   double precision, dimension(Nz):: z,phiLW
   double precision, dimension(1:Nxb,1:Nyb,1:Nz):: W, ADVw, Funp, ADVzW, Rxn_zW, LW_z
   double precision, dimension(0:Nxb+1,0:Nyb+1,1:Nz):: Wg, Wp, Wold, Wnew, Wnewg, Wold1g
   double precision, dimension(-1:Nxb+2,-1:Nyb+2,1:Nz):: Wpv 
   double precision, dimension(0:Nxb+1,1:Nz):: Aw, Dw, Bw
   double precision, dimension(1:Nyb):: BCxRHSw
   double precision, dimension(1:Nxb,1:Nz):: BCyRHSw
   double precision, dimension(0:Nxb+1,0:Nyb+1):: Wold1
   double precision, dimension(0:Nyb+1)::BCxW
   double precision, dimension(1:Nxb,1:Nyb):: ADVtemp
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
   double precision, dimension(1:Nxb,1:Nyb,1:Nz):: Vc, ADVV, ADVzV, Wz, RXNv, Vold, LWV_z
   double precision, dimension(0:Nxb+1,0:Nyb+1,1:Nz):: Vg, Vp, Vnew, Vnewg
   double precision, dimension(-1:Nxb+2,-1:Nyb+2,1:Nz):: Vpv 
   double precision, dimension(0:Nxb+1,1:Nz):: AV, DV, BV
   double precision, dimension(1:Nyb):: BCxRHSV
   double precision, dimension(1:Nxb,1:Nz):: BCyRHSV
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
   double precision, dimension(0:Nxb+1, 0:Nyb+1):: branch,branchs, branchg, branchold, Branchg2, Branch2
   double precision, dimension(0:Nxb+1, 0:Nyb+1):: Bsg, Bp, Bsp, Branchsold
   double precision, dimension(Nxb,Nyb):: ADVb, RXNb
   double precision:: Bup
   double precision, dimension(-1:Nxb+2,-1:Nyb+2):: Bsv
   double precision, dimension(0:Nyb+1)::BCxB
   double precision, dimension(0:Nxb+1)::Ab, Bb, Db
   double precision, dimension(1:Nxb,1:Nyb):: Wintg, Vintg
   double precision, dimension(1:Nz+1):: htrap
   double precision, dimension(1:Nz+1):: Vex,Wex
   !---------------------------------------------
   double precision, dimension(0:Nxb+1,0:Nyb+1):: S10c, S10old
   double precision, dimension(1:Nxb,1:Nyb):: Pore, Diam
   double precision:: FC2, FC1 
   double precision, dimension(1:Nxb,1:Nyb):: KBmat, KLmat
   !--------------------Bump-------------------
   integer:: btype
   double precision, dimension(Nxb, Nyb):: BumpCon

   !---Hindered Variables-----------------------
   double precision, dimension(0:Nxb+1,0:Nyb+1):: phi_f, phi_g, rf, rc, rz
   double precision, dimension(-1:Nxb+2,-1:Nyb+2):: phi_gpv
   double precision, dimension(-1:Nxb+2,-1:Nyb+2):: Huc_lr, Huc_tb, Huz_lr, Huz_tb, Huf_lr, Huf_tb
   double precision, dimension(0:Nxb+1,0:Nyb+1):: Hdc_lr, Hdc_tb, Hdf_lr, Hdf_tb, Hdz_lr, Hdz_tb
   double precision, dimension(-1:Nxb+2,-1:Nyb+2):: uvolc_lr, uvolf_lr, uvolz_lr
   double precision, dimension(-1:Nxb+2,-1:Nyb+2):: uvolc_tb, uvolf_tb, uvolz_tb
   double precision, dimension(-1:Nxb+2,-1:Nyb+2):: vvolc_lr, vvolf_lr, vvolz_lr
   double precision, dimension(-1:Nxb+2,-1:Nyb+2):: vvolc_tb, vvolf_tb, vvolz_tb
   integer, dimension(2):: Tg_minloc
   integer:: hindif_on, fluidon, nobump
   !----Saving files----------------                          
   character(100):: Nyb_char, Dtc_char
   character(6)::fileend
   character(50)::Zfile, Efile, Ffile, Rfile, Rgfile,Tfile, Sfile, Vfile
   character(50):: Wfile,Bfile, Thfile, Bgfile, Thgfile, Bsfile, Tsfile, Parsfile
   character(50):: Velocityufile, Velocityvfile,Alpha2file 
   integer:: Wunit=47, Eunit=48, Vunit = 49, Runit=50, Rgunit=51, Tunit=52, Sunit=53
   integer:: Thunit=54, Bunit=55, Thgunit=56, Bgunit=57, Bsunit = 58, Tsunit=59
   integer:: Z2unit= 61, Funit = 60, Parsunit = 61, Velocityuunit=62, Velocityvunit=66
   integer:: Alpha2unit=67
   real :: asq(1:nxb+1, 1:nyb+1, 2)
   real :: fbg(0:nxb+1,0:nyb+1,2)
   External:: DGTSV
!------------------------------
   real :: mgtol   = 1e-10
   real :: difftol = 1e-7
 CONTAINS
!----------------------------------------------------------------------------------------
subroutine driver
write(*,*) 'Nz=', Nz, 'Nx=', Nxb, 'Ny=', Nyb
write(*,*) 'begin next file'
fileend  ='0312_4'
nobump   = 0  !if nobump = 0 -> yes bump, nobump = 1 -> no bump
btype    = 1
NxbD     = dble(Nxb)
NybD     = dble(Nyb)
NzD      = dble(Nz)
Ls       = pltdiam!1.00d-2
ZupD     = 1d-9 ! The dimensional upstream prothrombic concentration   
DD       = 1d-7 !Dimensional Diffusion   
alphao   = 1d1! 1d2!1d+00 !1d+03 !BC term
dalpha   = alphao*DD*ZupD/(Ls)  ! Dimensional alpha
write(*,*) 'dalpha=', dalpha
alpha    = 0d0 !alphao   
beta     = 5d-2 !5d-4!(alpha(1)*(1d0-1d-1*2d0/pi)-1d-1*(1d0-2d0*1d-1/pi))/
!BC term
dbeta    = beta*ZupD
kcat     = 1d-7
Ro       = dalpha/kcat
write(*,*) 'Ro=', Ro, 'kcat=', kcat, 'koffish=', dbeta*kon-kcat
!koff     = 1d-1 !dbeta*kon - kcat
kon      = 1d6!(kcat+koff)/dbeta
koff     = dbeta*kon-kcat
gamma    = 1d1
Ly       = dheight/xchar !dheight/xchar !2.00d0
Lx       = dble(maspect)*Ly !d-
write(*,*) 'Lx =', Lx, 'Ly=', Ly
!stop
D        = 1d0 !DD/(gamma*LS**2)!1.00d-4/(gamma*LS**2d0)!1.00d-07/(gamma*LS**2)
muc       = 6d0 !Sharpness of Alpha transition between 0 and alpha_o
Kappaon  = ZupD*kon*Ls**2d0/DD !Kappaon2 = 0d0 !kon*Ls*Ro/DD                
Kappaoc  = (koff+kcat)*Ls**2d0/DD
Kappaoff = -koff*Ro*Ls/(ZupD*DD) !1d8*(koff*Ro*Ls/DD)
kld      = 4.8d1  !I have no idea what numbers these should be!!!
kbd      = 4.8d3  ! I have no idea what numbers these should be!!
Lz       = 1d0
Zup      = 1d0
Eup      = 0d0
Fup      = 1d0!/2d1 !1d0
Fupd     = 1d-8
kf       = 1d2
kfs      = FupD

write(*,*) 'FupD=', FupD
!stop
Kappaf   = kf*Zupd*Ls**2/(DD*Fupd)
Kappafs  = Kfs/FupD
kat      = 1d-1*Ls**2/DD
nu       = 3d0!4d0  
kl       = 1d4!5d-1!1d0 !Ls**2.*Fupd*kld/DD!1d0!4.8d3*Ls**2/DD
kb       = kl!/1d1!/2d0!/1d1!/2d0 !<1d1!Ls**2.*Fupd**2.*kbd/DD!5d-1 ! 4.8d1*Ls**2/DD
kbMat    = kb
klMat    = kl

!------Pore/diam size parameters -------------
FC1 = 6.02e-15!FupD*6.02d23 
FC2 = 4.275d-19*6.02d23*FupD !Volume of Fibrin monomer xavogadrosxFupD


write(*,*) 'kb=', kb, 'kappaf=', kappaf, 'kappafs', kappafs
!stop
!write(*,*) 'Lx=', Lx, 'Ly=', Ly
!stop
!Initial Values
W   = 0d0
!R = 0.0
chi = 1d0
!-----------FK Parameters------------------------------------ 
!Ly = 2d-1
!Lx = 2d-1
!D  = 1d-2
JK = -1 !-1=cheryl's bcs and diffusion, 1=FK bcs and CZ Diff, 2=FK BC and Diff     
!The boundary conditions also have to be altered in both Ghost_cells 
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
fluidon    = 1
bigdomain = 0
do i = 75, Nxb
   if (bigdomain ==1) then
      kbMat(i,1:Nyb)=kb*5d-1*(1 - tanh(10d0*(x(i) - 750d-2)))  !325
      klMat(i,1:Nyb)=kl*5d-1*(1 - tanh(10d0*(x(i) - 750d-2)))  !325
   else
      kbMat(i,1:Nyb)=kb*5d-1*(1 - tanh(10d0*(x(i) - 325d-2)))  !325
      klMat(i,1:Nyb)=kl*5d-1*(1 - tanh(10d0*(x(i) - 325d-2)))  !325
   end if
end do
write(*,*) 'kbmat(1:Nxb,1)=', kbmat(1:Nxb,1)



to = 0d0
write(*,*) 'dt = ', dt, min(hy/1d5,minval(hz)/8d0,hx/1d5)
dt_fluid = 1d0*hb/(3d0*umax/uchar)
if (fluidon ==1) then
   dt = min(dt_fluid,hy/1d4,minval(hz)/8d0,hx/1d4) !1d0*hb/(4d0*umax/uchar) !min(hy/5d4,minval(hz)/8d0,hx/5d4)!minval(hz)/8d0 !min(hy/1d3,minval(hz)/8d0,hx/1d3) !minval(hz)/4d0
else 
   dt = min(hy/5d3,minval(hz)/8d0,hx/5d3)
end if
tend = 1d0  !dble(tfinal)!1d-1 !0.8!1d-1!2.5d-1!2d0*dt!8d-1 !1d0!5d0*dt !1d2* dt !5d-2
!tend = tend+dt/
timerecord= 0d0 !1d-1
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
write(*,*) 'hx=', hx, 'dt=', dt
BCxz  = -1d0

L2 =  1.6!Lx/2d0 - Lx/1d1
L3 =  2.4!Lx/2d0 + Lx/1d1
write(*,*) 'L2=', L2, 'L3=', L3, 'Lx=', Lx, 'Ly=', Ly
write(*,*) 'hx = ', hx, 'x=', maxval(x), 'y=', maxval(y)
!stop

!alpha = alphao*(5d-1*tanh(muc*(xo-L2))-5d-1*tanh(muc*(xo-L3)))
!write(*,*) 'alphao=', alphao, 'alpha=', maxval(alpha)              
!stop
call alpha_function(btype)
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
Fpc = Fup
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
thetasold = 0d0
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
branchsold=0d0
S10old    = 0d0
uvolz_lr     = 0d0
vvolz_tb     = 0d0
!----------------------------------------------
Hdc_lr=1d0
Hdc_tb=1d0
Hdf_lr=1d0
Hdf_tb=1d0
Hdz_lr=1d0
Hdz_tb=1d0
Huc_lr=1d0
Huc_tb=1d0
Huz_lr=1d0
Huz_tb=1d0
Huf_lr=1d0
Huf_tb=1d0
uvolc_lr=0d0
uvolf_lr=0d0
uvolz_lr=0d0
uvolc_tb=0d0
uvolf_tb=0d0
uvolz_tb=0d0
vvolc_lr=0d0
vvolf_lr=0d0
vvolz_lr=0d0
vvolc_tb=0d0
vvolf_tb=0d0
vvolz_tb=0d0
call ghost_cells(Z2p, Az, Bz,Dz, Z2g, Nxb, Nyb,1, zup,jk,uvolz_lr(1,0:Nyb+1), vvolz_tb(0:Nxb+1,1))
write(Dtc_char,*) dtc
!write(*,*) 'asq=', asq
!fluidon    = 0 !1
if (fluidon==0) then
   Ulr =  0d0 
else
   Ulr= u(-1:nxb+2,-1:nyb+2,1)!0d0
endif
!Ulr = 0d0
!Utb = 0d0
!Vlr = 0d0
!Vtb = 0d0
!if (flowon==1) then
!   do i = -1, Nxb+2
!      Ulr(i,1:Nyb) = Ls**2/DD*gamma*y(1:Nyb)*(1d0-y(1:Nyb)/Ly)!y*(4d0-8d0*y) !1-(1-yv(0:Ny+1))**2
!   end do
!   do i = -1, Nxb+2
!      Utb(i,1:Nyb) = Ls**2/DD*gamma*(y(1:Nyb)-hy/2)*(1d0-(y(1:Nyb)-hy/2)/Ly)!1-(1-(yv(0:Ny+1)-hy/2))**2
!   end do
!end if

Utb(0:Nxb+1,0:Nyb+1) = 25d-2*(Ulr(0:nxb+1,0:nyb+1) + Ulr(1:nxb+2,0:nyb+1) + Ulr(1:nxb+2,-1:nyb) + Ulr(0:nxb+1,-1:nyb)) !0d0
Utb(-1,0:Nyb+1)      = Utb(0,0:Nyb+1)
Utb(Nxb+2,0:Nyb+1)   = Utb(Nxb+1, 0:Nyb+1)
Utb(-1:Nxb+2,-1)     = 0d0
Utb(-1:Nxb+2,Nyb+2)  = 0d0
write(*,*) 'Ls=', Ls, 'gamma=', gamma, 'Dd=', DD, 'Ly=', Ly, 'y=', y
write(*,*)'Ulr=', Ulr(1, 0:Nyb+1)
write(*,*) 'Utb=', Utb(1,1:Nyb+1)

Vtb = u(-1:nxb+2,-1:nyb+2,2)!0d0
!write(*,*) 'Vtb=', maxval(Vtb)
Vlr(0:Nxb+1,0:Nyb+1) = 25d-2*(Vtb(0:nxb+1,0:nyb+1) + Vtb(-1:nxb,0:nyb+1) + Vtb(0:nxb+1,1:Nyb+2) + Vtb(-1:nxb,1:nyb+2)) !0d0
Vlr(-1,0:Nyb+1)      = 0d0
Vlr(Nxb+2,0:Nyb+1)   = Vlr(Nxb+1,0:Nyb+1)
Vlr(-1:Nxb+2,-1)     = 0d0
Vlr(-1:Nxb+2,Nyb+2)  = 0d0
write(*,*) 'velocity,u=', maxval(Utb), maxval(Ulr), minval(Utb), minval(Ulr)
write(*,*) 'velocity,V =', maxval(Vtb), maxval(Vlr)
asq=0.0
t = to
dto = dt
write(*,*) TRIM(fileend), '.dat'

write(Parsfile, '(25a)') 'Parm_',TRIM(fileend),'.dat'
open(Parsunit, file=Parsfile)

write(Parsunit,22) "shear rate", 4.0*umax/(ymax*1d-2) 
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

write(Velocityufile,'(25a)') 'U_Nxb66_',fileend,'.dat'
open(Velocityuunit,file=Velocityufile, form='binary')

write(Velocityvfile,'(25a)') 'Vvel_Nxb66_',fileend,'.dat'
open(Velocityvunit,file=Velocityvfile, form='binary')

write(Alpha2file,'(25a)') 'Alpha2_Nxb66_',fileend,'.dat'
open(Alpha2unit,file=Alpha2file, form='binary')
do j = 1,Nyb
   write(Bsunit)  Branchs(0:Nxb+1,j)
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
   write(Velocityuunit)  Ulr(1:Nxb,j)  
   write(Velocityvunit)  u(1:Nxb,j,2)!Vtb(1:Nxb,j)
   write(Alpha2unit) asq(1:Nxb,j,1)
end do
do jz = 1,Nz
   do j = 1,Nyb
       write(Wunit)   W(1:Nxb,j,jz)
       write(Vunit)   Vnew(1:Nxb,j,jz)
    end do
end do
nt = 0
write(Tunit)  to
write(*,*) 'dt=', dt
dogelation = 1
dosource   = 1
hindif_on  = 0 !1 means on

!n = 0
 !call new_alpha(asq)
fbg(0:nxb+1,0:nyb+1,1) = 0.0
fbg(0:nxb+1,0:nyb+1,2) = 0.0
fns(1:nxb+1,1:nyb+1,:) = 0.0
if (fluidon ==1) then
   dt = min(dt_fluid,hy/1d4,minval(hz)/8d0,hx/1d4) 
else
   dt = min(hy/5d3,minval(hz)/8d0,hx/5d3)
end if
call Bump_type(btype)
call alpha_sq(t)
!stop
!--------------------Begin Time Steps-----------------------------------
write(*,*) 'pre time step=', dt
do while  (t< tend)
   write(*,*) 'dt = ', dt, 'CFL=', CFL
!stop
   t = t + dt
   !write(*,*) 'kx=', kx
   if (fluidon==1) then
   call alpha_sq(t)
   write(*,*) 'maxval(asq)', maxval(asq)
   call new_alpha(asq)
   call fluidstep!(asq)
   Ulr = u(-1:nxb+2,-1:nyb+2,1)!0d0
   Utb(0:Nxb+1,0:Nyb+1) = 25d-2*(Ulr(0:nxb+1,0:nyb+1) + Ulr(1:nxb+2,0:nyb+1) + Ulr(1:nxb+2,-1:nyb) + Ulr(0:nxb+1,-1:nyb)) !0d0
   Utb(-1,0:Nyb+1)      = Utb(0,0:nyb+1)
   Utb(Nxb+2,0:Nyb+1)   = Utb(Nxb+1, 0:Nyb+1)
   Utb(-1:Nxb+2,-1)     = 0d0
   Utb(-1:Nxb+2,Nyb+2)  = 0d0

   Vtb = u(-1:nxb+2,-1:nyb+2,2)!0d0
   Vlr(0:Nxb+1,0:Nyb+1) = 25d-2*(Vtb(0:nxb+1,0:nyb+1) + Vtb(-1:nxb,0:nyb+1) + Vtb(0:nxb+1,1:Nyb+2) + Vtb(-1:nxb,1:nyb+2)) !0d0
   Vlr(-1,0:Nyb+1)      = 0d0
   Vlr(Nxb+2,0:Nyb+1)   = Vlr(Nxb+1,0:Nyb+1)
    Vlr(-1:Nxb+2,-1)     = 0d0
   Vlr(-1:Nxb+2,Nyb+2)  = 0d0
   else 
      Ulr = 0d0
      Utb = 0d0
      Vlr = 0d0
      Vtb = 0d0
      end if
  !do j = 1,Nyb
   !   write(Velocityuunit)  Ulr(1:Nxb,j)
   !   write(Velocityvunit)  Vtb(1:Nxb,j)
   !end do
   !Need to do the ghost_cells for phi_g that was we can find the hindered diffusion terms.
   phi_g = 0.0!Thetag*FC2
   phi_f = 1.0-phi_g


    if (hindif_on==1) then
      call ghost_cells_phig(phi_g)
      phi_gpv=0.0
      phi_gpv(0:Nxb+1,0:Nyb+1) = phi_g
      call hinddifcoef
   end if

!   write(*,*) 'max=', maxval(phi_g)
!   write(*,*) 'maxf=', maxval(phi_f)
!   write(*,*) 'maxHu=', maxval(Huc_lr), maxval(Huc_TB)
   !write(*,*) 'maxHd=', maxval(Hdc_lr), maxval(Hdc_TB), minval(Hdc_lr), minval(Hdc_TB)
!   write(*,*) '------------------------------------------------------------------'
   if (dosource==1) then
      if (fluidon==0) then
         if (t.gt.0d0) then !1d-1) then
            alphat = alpha*(1d0-exp(-nu*(t)))!-1d-1)))
            !write(*,*) 'max=', maxval(alphat)
         else
            alphat = 0d0
         end if
      else
         if (t.gt.5d-2) then   
            alphat = alpha*(1d0-exp(-nu*(t-5d-2)))
            !write(*,*) 'max=', maxval(alphat
         else
            alphat = 0d0
         end if

   end if
      betat   = phi_f(0:Nxb+1,0)*beta + 3d0/4d0*Z2g(0:Nxb+1,1) -1d0/8d0*Z2g(0:Nxb+1,2)&
           + 3d0/8d0*Z2g(0:Nxb+1,0)
      Az      = (8d0*betat - 6d0*alphat*hy)/(8d0*betat+3d0*alphat*hy)!1d0 
      Dz      = alphat*hy/(8d0*betat + 3d0*alphat*hy)
      Bz      = -hy*vvolz_tb(0:Nxb+1,1)*(15.0/8.0*Z2g(0:Nxb+1,1) + 10.0/8.0*Z2g(0:Nxb+1,2) + 3d0/8d0*Z2g(0:Nxb+1,3))*betat&
              /(betat+3d0/8d0*alphat*hy)

      BCyRHSz = D(1:Nxb,1)*Hdz_tb(1:Nxb,1)*ky*Bz(1:Nxb)
      BCxRHSz = 2d0*D(1,1:Nyb)*Hdz_lr(1,1:Nyb)*kx*zup

      Rkz = 0d0

      call ghost_cell_vec(Z2g,Z2pv, Az,Bz,Dz,Nxb,Nyb,1,Zup,jk,uvolz_lr(1,0:nyb+1),vvolz_tb(0:Nxb+1,1))
      call advection(Z2pv,Huz_lr*Ulr+uvolz_Lr,Huz_tb*Utb+uvolz_tb,Huz_lr*Vlr+vvolz_lr,Huz_tb*Vtb+vvolz_tb, dt,hx,hy,Nxb,Nyb, ADVz2)
      call ghost_cells(Z2g,Az,Bz,Dz,Z2p,Nxb,Nyb,1,zup,jk,uvolz_lr(1,0:nyb+1),vvolz_tb(0:Nxb+1,1))
      Z2g = Z2p

      call diffusion_subr(Nxb,Nyb,kx,ky,Zup,BCxz,BCxRHSz,Az,Dz,Bz,BCyRHSz,Hdz_TB*D,Hdz_LR*D,Z2g,Z2,ADVz2,Rkz,uvolz_lr(1,0:Nyb+1), vvolz_tb(0:Nxb+1,1))
      Z2  = Z2g(1:Nxb,1:Nyb)

      Be = hy*alphat*(3d0/4d0*Z2g(0:Nxb+1,1)-1d0/8d0*Z2g(0:Nxb+1,2) &
           + 3d0/8d0*Z2g(0:Nxb+1,0))/(betat) + vvolz_tb(0:Nxb+1,1)*(15.0/8.0*E2g(0:Nxb+1,1) + 10.0/8.0*E2g(0:Nxb+1,2) &
           + 3d0/8d0*E2g(0:Nxb+1,3))
      

      BCxRHSe = 2d0*D(1,1:Nyb)*Hdz_lr(1,1:Nyb)*kx*Eup
      BCyRHSe = D(1:Nxb,1)*Hdz_tb(1:Nxb,1)*ky*Be(1:Nxb)

      BCxe    = -1d0
      !------------------------                                          
      RKe = -dt*kat*(E2(1:Nxb,1:Nyb) - 5d-1*dt*kat*E2(1:Nxb,1:Nyb))
      ADVe = 0d0

      call ghost_cell_vec(E2g,E2pv, Ae,Be,De,Nxb,Nyb,1,Eup,jk,uvolz_lr(1,0:nyb+1),vvolz_tb(0:Nxb+1,1))
      call advection(E2pv,Huz_lr*Ulr+uvolz_Lr,Huz_tb*Utb+uvolz_tb,Huz_lr*Vlr+vvolz_lr,Huz_tb*Vtb+vvolz_tb,dt,hx,hy,Nxb,Nyb, ADVe)
      call ghost_cells(E2p,Ae,Be,De,E2g,Nxb,Nyb,1,Eup,jk,uvolz_lr(1,0:nyb+1),vvolz_tb(0:Nxb+1,1))      
      call diffusion_subr(Nxb,Nyb,kx,ky,Eup,BCxe,BCxRHSe,Ae,De,Be,BCyRHSe,Hdz_TB*D,Hdz_LR*D,E2g,E2,ADVe,RKe,uvolz_lr(1,0:Nyb+1), vvolz_tb(0:Nxb+1,1))
      call ghost_cells(E2g,Ae,Be,De,E2p,Nxb,Nyb,1,Eup,jk,uvolz_lr(1,0:nyb+1),vvolz_tb(0:Nxb+1,1))

      E2g = E2p
      E2  = E2g(1:Nxb,1:Nyb)
!      write(*,*) 'E2=', maxval(E2), minval(E2)
!     write(*,*) 'end E'
      !---------------------Advect and Diffusion F2------------
      Bf =vvolf_tb(0:Nxb+1,1)*(15.0/8.0*Fg(0:Nxb+1,1) + 10.0/8.0*Fg(0:Nxb+1,2) &
           + 3d0/8d0*Fg(0:Nxb+1,3))
      BCxRHSf = 2d0*D(1,1:Nyb)*Hdf_lr(1,1:Nyb)*kx*Fup
      BCyRHSf = D(1:Nxb,1)*Hdz_tb(1:Nxb,1)*ky*Bf(1:Nxb)
      ADVF = 0d0
      RKF  = 0d0

      call ghost_cell_vec(Fg,Fpv, Af,Bf,Df,Nxb,Nyb,1,Fup,jk,uvolf_lr(1,0:nyb+1),vvolf_tb(0:Nxb+1,1))
      call advection(Fpv,Huf_lr*Ulr+uvolf_Lr,Huf_tb*Utb+uvolf_tb,Huf_lr*Vlr+vvolf_lr,Huf_tb*Vtb+vvolf_tb, dt,hx,hy,Nxb,Nyb, ADVf)

      Fg = Fpv(0:Nxb+1,0:Nyb+1)


      call ghost_cells(Fg,Af,Bf,Df,Fpc,Nxb,Nyb,1,Fup,jk,uvolf_lr(1,0:nyb+1),vvolf_tb(0:Nxb+1,1))
      call RungaKutta(E2g, Fpc, kappaf, kappafs, Nxb, Nyb, dt, Rkf)
      call ghost_cells(Fpc,Af,Bf,Df,Fg,Nxb,Nyb,1,Fup,jk,uvolf_lr(1,0:nyb+1),vvolf_tb(0:Nxb+1,1))
      call diffusion_subr(Nxb,Nyb,kx,ky,Fup,BCxf,BCxRHSf,Af,Df,Bf,BCyRHSf,Hdf_TB*D,Hdf_LR*D,Fg,F,ADVf, RKf,uvolf_lr(0:Nxb+1,1), vvolf_tb(1,0:Nyb+1))! 
      F  = Fg(1:Nxb,1:Nyb)
     write(*,*) 'max F=', maxval(F), minval(F)
      S10c = kappaf*E2g*Fg/(kappafs+Fg)
   end if

   !---------------------Gelation Step-----------------------------------------------
   if (dogelation ==1) then
      call beam_warming(ADVzW,Funp,1)
      call Lax_wendroff(LW_z)
      
      phiLW = 5d-1*tanh(2d1*(z-25d-2))+5d-1 
      ADVzW(1:Nxb,1:Nyb,1) = LW_z(1:Nxb,1:Nyb,1)                                  
         do i = 2,Nz-1
            ADVzW(1:Nxb,1:Nyb,i) = ADVzW(1:Nxb,1:Nyb,i)*phiLW(i) + LW_z(1:Nxb,1:Nyb,i)*(1d0-phiLW(i))
         end do
         ADVzW(1:Nxb,1:Nyb,Nz)= ADVzW(1:Nxb,1:Nyb,Nz)
      !---Check CFL----------                                                                              
      q=maxval(abs(Funp))

      cfl = q*dt/minval(hz)
      !dt_fluid = 3d0*hx/(4d0*max(maxval(Ulr),maxval(Utb), maxval(Vlr),maxval(Vtb))/uchar)
      !write(*,*) 'dt_fluid=', dt_fluid
      if (cfl>1d-12) then
         if (cfl <=2d-1) then
            dt1 = dt*(4.9d-1/cfl)
            dt = min(dt1,dto,dt_fluid);
         else if (cfl >=5d-1) then
            dt = dt*(4.9d-1)/cfl
         end if
      !else
       !  dt = !dto
      end if

      kx = dt/(2d0*hx**2.)
      ky = dt/(2d0*hy**2.)
      !-----Update W at z=1 ---------------------      
      !This line of coded should be added somehow to the RXNW subroutine so I don't have to
      !change the reaction equation in 2 places

      Wnew(1:Nxb,1:Nyb,Nz) = Wold(1:Nxb,1:Nyb,Nz) -dt*ADVzW(1:Nxb,1:Nyb,Nz) + dt*(-klMat/phi_f(1:Nxb,1:Nyb)*z(Nz)*chi*Rgold(1:Nxb,1:Nyb)**2.) &
                          + dt*(kbMat/(2d0*phi_f(1:Nxb,1:Nyb)**2)*z(Nz)**2.*Rold(1:Nxb,1:Nyb)**3.) &
                          - dt*kbMAT/(2d0*phi_f(1:Nxb,1:Nyb)**2)*z(Nz)*(-Rold(1:Nxb,1:Nyb)**3. + 3d0*chi*Rsold(1:Nxb,1:Nyb)*Rgold(1:Nxb,1:Nyb)**2. &
                          + chi*Rgold(1:Nxb,1:Nyb)**3.)
      !------Update R using W @ z=1-----------------
      !--Rg n+1--------------   
      Rg(1:Nxb,1:Nyb)  = -Wnew(1:Nxb,1:Nyb,Nz)

      Br = vvolc_tb(0:Nxb+1,1)*(15.0/8.0*Rsold(0:Nxb+1,1) + 10.0/8.0*Rsold(0:Nxb+1,2) &
           + 3d0/8d0*Rsold(0:Nxb+1,3)) 
      BCxRHSr = 2d0*D(1,1:Nyb)*Hdc_lr(1,1:Nyb)*kx*0d0
      BCyRHSr = D(1:Nxb,1)*Hdc_tb(1:Nxb,1)*ky*Br(1:Nxb)
      call ghost_cell_vec(Rsold, Rsv,Ar, Br, Dr, Nxb, Nyb, 1, Rup,jk,uvolc_lr(1,0:nyb+1),vvolc_tb(0:Nxb+1,1))
      call advection(Rsv, Huc_lr*Ulr+uvolc_Lr,Huc_tb*Utb+uvolc_tb,Huc_lr*Vlr+vvolc_lr,Huc_tb*Vtb+vvolc_tb, dt, hx, hy, Nxb, Nyb, ADVr)
      call RXN_R(rxnr)
      
      !!!You can't use Ghost_cells subroutine for Rg - Rg doesn't have BC its must go through interpolation
      call gel_ghost_cells(Rg,Rgph, Nxb, Nyb)
      Rg = Rgph
      call gel_ghost_cells(Rgold,Rgph, Nxb, Nyb) 
      Rgold = Rgph
      
      
      Br = vvolc_tb(0:Nxb+1,1)*(15.0/8.0*Rsold(0:Nxb+1,1) + 10.0/8.0*Rsold(0:Nxb+1,2)) +  (2d0*Rg(0:Nxb+1,1)-3d0*Rg(0:Nxb+1,2)+Rg(0:Nxb+1,3))
           !+ 3d0/8d0*Rsold(0:Nxb+1,3)) + (71d0*Rgold(0:Nxb+1,1) - 138d0*Rgold(0:Nxb+1,2) + 88d0*Rgold(0:Nxb+1,3) - 21d0*Rgold(0:Nxb+1,4))/96d0  
      call ghost_cells(Rold, Ar, Br, Dr, Rph, Nxb, Nyb, 1, Rup,jk,uvolc_lr(1,0:nyb+1),vvolc_tb(0:Nxb+1,1)) 
      BCyRHSr = D(1:Nxb,1)*Hdc_tb(1:Nxb,1)*ky*Br(1:Nxb)
      call diffusion_subr_r(Nxb,Nyb,kx,ky,Rup,BCxR,BCxRHSR,AR,DR,BR,BCyRHSR,Hdc_TB*D,Hdc_LR*D,Rg,Rph,Rgold,RXNr,ADVr,dt,jk,uvolc_lr(1,0:Nyb+1), vvolc_tb(0:Nxb+1,1))
      R  = Rph
      Rs = R-Rg
      

      !----Diffusion and Advect W in x,y -------        
      do idx_Bw=1,Nz
      Bw(0:Nxb+1,idx_Bw) = vvolc_tb(0:Nxb+1,1)*(15.0/8.0*W(0:Nxb+1,1,idx_BW) + 10.0/8.0*W(0:Nxb+1,2,idx_Bw) &
           + 3d0/8d0*W(0:Nxb+1,3,idx_Bw)) 
      end do
      BCxRHSw = 2d0*D(1,1:Nyb)*Hdc_lr(1,1:Nyb)*kx*0d0
      do idx_Bw=1,Nz
      BCyRHSw(1:Nxb,idx_Bw) = D(1:Nxb,1)*Hdc_tb(1:Nxb,1)*ky*Bw(1:Nxb,idx_Bw)
      end do
      do m= 1,Nz
         Wg(1:Nxb,1:Nyb,m) = Wold(1:Nxb,1:Nyb,m) - z(m)*(Wold(1:Nxb,1:Nyb,Nz)+Wnew(1:Nxb,1:Nyb,Nz))/2d0
      end do

      call ghost_cell_vec(Wg,Wpv,Aw,Dw,Bw,Nxb,Nyb,Nz,Wup,jk,uvolc_lr(1,0:nyb+1),vvolc_tb(0:Nxb+1,1))

      do m = 1, Nz-1

         call advection(Wpv(-1:Nxb+2,-1:Nyb+2,m), Huc_lr*Ulr+uvolc_Lr,Huz_tb*Utb+uvolc_tb,Huc_lr*Vlr+vvolc_lr,Huc_tb*Vtb+vvolc_tb,&
                       dt, hx, hy, Nxb, Nyb, ADVtemp)
         
         ADVw(1:Nxb,1:Nyb,m) = ADVtemp
      end do
      
      ADVw(1:Nxb,1:Nyb,Nz) = 0d0
      
      Rgavg = Rgold!(Rg+Rgold)/2d0 !Rgold
      Ravg  = Rold !(R+Rold)/2d0 !Rold    
      Rsavg = Rsold!(Rs+Rsold)/2d0 !Rsold
      Rxn_zW = 0d0
      call RXNW(Rxn_zW)
      !-------To extend the domain of W must us BC => W(y=-1) = W(y=1) + z h_y dRg/dy(at y=0)---------------------
      do idx_Bw=1,Nz
      Bw(0:Nxb+1,idx_Bw) = vvolc_tb(0:Nxb+1,1)*(15.0/8.0*W(0:Nxb+1,1,idx_BW) + 10.0/8.0*W(0:Nxb+1,2,idx_Bw) &
           + 3d0/8d0*W(0:Nxb+1,3,idx_Bw)) - z(idx_Bw)*(2d0*Rg(0:Nxb+1,1)-3d0*Rg(0:Nxb+1,2)+Rg(0:Nxb+1,3))
      !*(71d0*Rgold(0:Nxb+1,1) -138d0*Rgold(0:Nxb+1,2) + 88d0*Rgold(0:Nxb+1,3) - 21d0*Rgold(0:Nxb+1,4))/96d0
      end do
      BCxRHSw = 2d0*D(1,1:Nyb)*Hdc_lr(1,1:Nyb)*kx*0d0
      do idx_Bw=1,Nz
      BCyRHSw(1:Nxb,idx_Bw) = D(1:Nxb,1)*Hdc_tb(1:Nxb,1)*ky*Bw(1:Nxb,idx_Bw)
      end do
      
      Wp(1:Nxb,1:Nyb,1:Nz) = Wold(1:Nxb,1:Nyb,1:Nz)
      Wold1g(1:Nxb,1:Nyb,1:Nz) = Wold(1:Nxb,1:Nyb,1:Nz)
      Wnewg(1:Nxb,1:Nyb,1:Nz) = Wnew(1:Nxb,1:Nyb,1:Nz)
      call gel_ghost_cells(Wold(0:Nxb+1,0:Nyb+1,Nz),Wold1g(0:Nxb+1,0:Nyb+1,Nz),Nxb,Nyb)
      call gel_ghost_cells(Wnew(0:Nxb+1,0:Nyb+1,Nz),Wnewg(0:Nxb+1,0:Nyb+1,Nz),Nxb,Nyb)
      call ghost_cells(Wp,  Aw, Bw, Dw,Wg,    Nxb,Nyb,Nz,Wup,jk,uvolc_lr(1,0:nyb+1),vvolc_tb(0:Nxb+1,1))
      
      
      do idx_Bw = 1,Nz
         Wg(0:Nxb+1,0,idx_Bw) = Wpv(0:Nxb+1,0,idx_Bw) - z(idx_Bw)*Rgold(0:Nxb+1,0) 
     end do
      call diffusion_subr_W(Nxb,Nyb,Nz,kx,ky,Wup,BCxW,BCxRHSW,Aw,Bw, Dw,BCyRHSW,Hdc_TB*D,Hdc_LR*D,Wg,z,&
                            Wold1g(0:Nxb+1,0:Nyb+1,Nz),ADVw, Rxn_zW, ADVzW, Wnewg,dt,jk,uvolc_lr(1,0:Nyb+1), vvolc_tb(0:Nxb+1,1))
      !Interpolation at z=0-------------------                     
      Wg(1:Nxb,1:Nyb,1) = Wg(1:Nxb,1:Nyb,2)*((z(1)-0)*(z(1)-z(3))*(z(1)-z(4)))/((z(2)-0)*(z(2)-z(3))*(z(2)-z(4)))&
                        + Wg(1:Nxb,1:Nyb,3)*((z(1)-0)*(z(1)-z(2))*(z(1)-z(4)))/((z(3)-0)*(z(3)-z(2))*(z(3)-z(4))) &
                        + Wg(1:Nxb,1:Nyb,4)*((z(1)-0)*(z(1)-z(2))*(z(1)-z(3)))/((z(4)-0)*(z(4)-z(2))*(z(4)-z(3)))
      Wnew = Wg
      W    = Wg(1:Nxb,1:Nyb,1:Nz)
       !----------------V equation -------------------------
      call beam_warming(ADVzV,Wz,2)
      call Lax_WendroffV(LWV_z)
      
      ADVzV(1:Nxb,1:Nyb,1) = LWV_z(1:Nxb,1:Nyb,1)
      do i = 2,Nz-1
         ADVzV(1:Nxb,1:Nyb,i) = ADVzV(1:Nxb,1:Nyb,i)*phiLW(i) + LWV_z(1:Nxb,1:Nyb,i)*(1d0-phiLW(i))
      end do
     
      Vg(1:Nxb,1:Nyb, 1:Nz) = Vold
      do idx_Bw=1,Nz
      BV(0:Nxb+1,idx_Bw) = vvolc_tb(0:Nxb+1,1)*(15.0/8.0*Vg(0:Nxb+1,1,idx_Bw) + 10.0/8.0*Vg(0:Nxb+1,2,idx_Bw) &
           + 3d0/8d0*Vg(0:Nxb+1,3,idx_Bw))
      end do
      BCxRHSV = 2d0*D(1,1:Nyb)*Hdc_lr(1,1:Nyb)*kx*0d0
      do idx_Bw=1,Nz
      BCyRHSV(1:Nxb,idx_Bw) = D(1:Nxb,1)*Hdc_tb(1:Nxb,1)*ky*Bv(1:Nxb,idx_Bw)
      end do
      call ghost_cell_vec(Vg,Vpv,AV,DV,BV,Nxb,Nyb,Nz,Vup,jk,uvolc_lr(1,0:nyb+1),vvolc_tb(0:Nxb+1,1))

      do m = 1, Nz
         call advection(Vpv(-1:Nxb+2,-1:Nyb+2,m), Huc_lr*Ulr+uvolc_Lr,Huc_tb*Utb+uvolc_tb,Huc_lr*Vlr+vvolc_lr,Huc_tb*Vtb+vvolc_tb,&
                       dt, hx, hy, Nxb, Nyb, ADVtemp)
         ADVV(1:Nxb,1:Nyb,m) = ADVtemp
         !if(maxval(abs(ADVV)) >1d-12) then
         !   write(*,*) 'Error ADVV'
         !end if
      end do

      Vp(1:Nxb,1:Nyb,1:Nz) = Vold
      call ghost_cells(Vp,  AV, BV, DV,Vg, Nxb,Nyb,Nz,Vup,jk,uvolc_lr(1,0:nyb+1),vvolc_tb(0:Nxb+1,1))

      call RXN_V(RxnV)
      call diffusion_subr_V(Nxb,Nyb,Nz,kx,ky,Vup,BCXV,BCxRHSV,AV,BV, DV,BCyRHSV,Hdc_TB*D,Hdc_LR*D,Vg,Vp, &
           ADVV, RxnV,ADVzV,dt,jk,uvolc_lr(1,0:Nyb+1), vvolc_tb(0:Nxb+1,1))
      
      Vg(1:Nxb,1:Nyb,1) = Vg(1:Nxb,1:Nyb,2)*((z(1)-0)*(z(1)-z(3))*(z(1)-z(4)))/((z(2)-0)*(z(2)-z(3))*(z(2)-z(4)))&
                        + Vg(1:Nxb,1:Nyb,3)*((z(1)-0)*(z(1)-z(2))*(z(1)-z(4)))/((z(3)-0)*(z(3)-z(2))*(z(3)-z(4))) &  
                        + Vg(1:Nxb,1:Nyb,4)*((z(1)-0)*(z(1)-z(2))*(z(1)-z(3)))/((z(4)-0)*(z(4)-z(2))*(z(4)-z(3)))
      Vc = Vg(1:Nxb,1:Nyb,1:Nz)
       !-------------------theta's and B's-------------------------
        !Need to find the integral of V (Vint) and W over z first----
      htrap(1)    = hz(1)/2d0
      htrap(2:Nz) = (hz(1:Nz-1)+hz(2:Nz))/2d0
      htrap(Nz+1) = hz(Nz)/2d0
      
      do i = 1,Nxb
         do j = 1,Nyb
            Vex(1)      = 0d0
            Vex(2:Nz+1) = Vc(i,j,1:Nz)
            Wex(1)      = 0d0
            Wex(2:Nz+1) = Wnew(i,j,1:Nz)
            Vintg(i,j)  = dot_product(htrap, Vex)
            Wintg(i,j)  = dot_product(htrap, Wex)
         end do
      end do

      thetas  = 0d0
      branchs = 0d0
      branchs(1:Nxb,1:Nyb) = -Rg(1:Nxb,1:Nyb) - 2d0*Wintg
      thetas(1:Nxb,1:Nyb)  = Vintg  + 2d0*branchs(1:Nxb,1:Nyb)
      
      RXNt  = dt*S10old(1:Nxb,1:Nyb)!(S10old(1:Nxb,1:Nyb)+S10old(1:Nxb,1:Nyb))/2d0
      Ravg  = R!(R+Rold)/2d0 !R
      Rsavg = Rs!(Rs+Rsold)/2d0 !Rsold!(Rs + Rsold)/2d0 !Rs
      Rgavg = Rg!(Rg+Rgold)/2d0 !Rgold !(Rg+ Rgold)/2d0 !Rg
      RXNb  = dt*kbMAT/(6d0*phi_f(1:Nxb,1:Nyb)**2)*(Ravg(1:Nxb,1:Nyb)**3. - chi*(3d0*Rsavg(1:Nxb,1:Nyb)*Rgavg(1:Nxb,1:Nyb)**2.+Rgavg(1:Nxb,1:Nyb)**3.))
      Tsp   = 0d0
      Bsp   = 0d0
      Tsp(1:Nxb,1:Nyb)  = (thetasold(1:Nxb,1:Nyb))!  + thetasold(1:Nxb,1:Nyb))/2d0
      Bsp(1:Nxb,1:Nyb)  = (branchsold(1:Nxb,1:Nyb))! + branchsold(1:Nxb,1:Nyb))/2d0 
      
      Bt = vvolc_tb(0:Nxb+1,1)*(15.0/8.0*Thetas(0:Nxb+1,1) + 10.0/8.0*Thetas(0:Nxb+1,2) &
           + 3d0/8d0*Thetas(0:Nxb+1,3))
      Bb = vvolc_tb(0:Nxb+1,1)*(15.0/8.0*Branchs(0:Nxb+1,1) + 10.0/8.0*Branchs(0:Nxb+1,2) &
           + 3d0/8d0*Branchs(0:Nxb+1,3))
      call ghost_cell_vec(Tsp, Tsv, At, Bt, Dtheta, Nxb, Nyb, 1, Tup,jk,uvolc_lr(1,0:nyb+1),vvolc_tb(0:Nxb+1,1))

      call advection(Tsv, Huc_lr*Ulr+uvolc_Lr,Huc_tb*Utb+uvolc_tb,Huc_lr*Vlr+vvolc_lr,Huc_tb*Vtb+vvolc_tb, dt, hx, hy, Nxb, Nyb, ADVt)


      call ghost_cell_vec(Bsp, Bsv, Ab, Bb, Db, Nxb, Nyb, 1, Bup,jk,uvolc_lr(1,0:nyb+1),vvolc_tb(0:Nxb+1,1))
      call advection(Bsv, Huc_lr*Ulr+uvolc_Lr,Huc_tb*Utb+uvolc_tb,Huc_lr*Vlr+vvolc_lr,Huc_tb*Vtb+vvolc_tb, dt, hx, hy, Nxb, Nyb, ADVb)

      Tsp = 0d0
      Bsp = 0d0

      Tsp(1:Nxb,1:Nyb) = (thetas(1:Nxb,1:Nyb))!old(1:Nxb,1:Nyb)+thetas(1:Nxb,1:Nyb))/2d0
      Bsp(1:Nxb,1:Nyb) = (branchsold(1:Nxb,1:Nyb)+branchs(1:Nxb,1:Nyb))/2d0
      
      call ghost_cells(Tsp,  At, Bt, Dtheta,Tsg, Nxb,Nyb,1,Tup,jk,uvolc_lr(1,0:nyb+1),vvolc_tb(0:Nxb+1,1))
      call ghost_cells(Bsp,  Ab, Bb, Db,    Bsg, Nxb,Nyb,1,Bup,jk,uvolc_lr(1,0:nyb+1),vvolc_tb(0:Nxb+1,1)) 

      write(*,*) 't=', t

      call diffusion_subr_tb(Nxb, Nyb, Thetaold,  theta,  Tsg, RXNt,kx,ky,Hdc_TB*D,Hdc_LR*D,ADVt,At,Bt,Dtheta,Tup,jk,uvolc_lr(1,0:Nyb+1), vvolc_tb(0:Nxb+1,1))
      call diffusion_subr_tb(Nxb, Nyb, Branchold, branch, Bsg, RXNb,kx,ky,Hdc_TB*D,Hdc_LR*D,ADVb,Ab,Bb,Db,    Bup,jk,uvolc_lr(1,0:Nyb+1), vvolc_tb(0:Nxb+1,1))
      

      Branchg(1:Nxb,1:Nyb) = Rg(1:Nxb,1:Nyb)
      Branch(1:Nxb,1:Nyb)  = Branchg(1:Nxb,1:Nyb) + Branchs(1:Nxb,1:Nyb)

      Thetag(1:Nxb,1:Nyb)  = Theta(1:Nxb,1:Nyb)  - Thetas(1:Nxb,1:Nyb)
      Tg_minloc = minloc(Thetag)
      write(*,*) 'Thetas=', minval(Thetas), minloc(Thetas), maxval(Thetas), maxloc(Thetas)

      call poresize(Pore, Diam)
   end if
    Thetaold  = Theta
    Branchold = Branch
    Wold1  = Wg(0:Nxb+1,0:Nyb+1,Nz)
    Wold   = Wg
    Rold   = R
    Rsold  = Rs
    Rgold  = Rg
    Vold   = Vc
    S10old = S10c
    nt     = nt+1

    Thetasold   = Thetas
    Branchsold  = Branchs

    if (t > timerecord) then
      do recordidx = 1, Nyb
         write(Tsunit)  Thetas(1:Nxb,recordidx)
         write(Thunit)  Theta(1:Nxb,recordidx)
         write(Thgunit) Thetag(1:Nxb,recordidx)
         write(Bsunit)  Branchs(0:Nxb+1,recordidx)
         write(Bunit)   Branch(1:Nxb,recordidx)
         write(Bgunit)  Branchg(1:Nxb,recordidx)
         write(Runit)   R(1:Nxb,recordidx)
         write(Rgunit)  Rg(1:Nxb,recordidx)
         write(Sunit)   S10c(1:Nxb,recordidx)
         write(Funit)   Fg(1:Nxb, recordidx)
         write(Eunit)   E2g(1:Nxb,recordidx)
         write(Z2unit)  Z2(1:Nxb,recordidx)
         write(Velocityuunit)  Ulr(1:Nxb,recordidx)
         write(Velocityvunit)   u(1:Nxb,recordidx,2)
         write(Alpha2unit) asq(1:Nxb,recordidx,1)
      end do
      write(*,*) 't=', t

      do jz =1,Nz
         do recordidx= 1,Nyb
            if (jz == Nz) then
            end if
            write(Wunit)   W(1:Nxb,recordidx,jz)
            write(Vunit)   Vc(1:Nxb,recordidx,jz)
         end do
      end do
      write(Tunit) t

         timerecord = timerecord + 1d-3

   end if

end do


do recordidx = 1, Nyb
   write(Bsunit)  Branchs(0:Nxb+1,recordidx)
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
   write(Velocityuunit)  Ulr(1:Nxb,recordidx)
   write(Velocityvunit)  u(1:Nxb,recordidx,2) !Vtb(1:Nxb,recordidx)
   write(Alpha2unit) asq(1:Nxb,recordidx,1)
end do
do jz = 1,Nz
   do recordidx = 1,Nyb
      write(Wunit)   W(1:Nxb,recordidx,jz)
      write(Vunit)   Vc(1:Nxb,recordidx,jz)
   end do
end do
write(tunit) t
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
!close(Vunit)
write(*,*) 'D=', d(1,1), 'dt=', dt, 'hx=', hx, 'kx=', kx
close(Wunit)
close(Z2unit)
write(*,*) 't=',t
close(Eunit)
close(Funit)
close(Velocityuunit)
close(Velocityvunit)
end subroutine driver
!--------------------------------------------------------------------------
subroutine fun_hz
integer:: I
double precision:: hzo, hmax, hend, hmin
double precision, dimension(Nz)::ztemp, hh 
write(*,*) 'Nz=', Nz

hzo = 1d0/(dble(Nz)-1d0)
ztemp = (/(hzo*dble(I),I = 0,Nz-1) /)


hmin  = 5d-3
hmax   = hzo
hend  = 1d-3
hh    = 5d-1*(1d0+tanh((7d-1-ztemp)/15d-2))

!hmin + (hmax-hmin) * 0.5 * (1d0+tanh((ztemp-0.1)/0.03)) &
!     - (hmax-hend) * 0.5 * (1d0+tanh(-(0.7-ztemp)/0.15))

hz(1:Nz)    = hh/sum(hh)


end subroutine fun_hz
!---------------------------------------------------------------------
subroutine ghost_cells(C,A,B, DC, Cnew,Nxb,Nyb,Nz,Cup,jk,uvol,vvol)

integer:: Nxb, Nyb, Nz,jk
double precision:: Cup
double precision, dimension(0:Nxb+1,0:Nyb+1,1:Nz)::C, Cnew
double precision, dimension(0:Nxb+1,1:Nz):: A, B, DC
double precision, dimension(0:Nxb+1):: vvol
double precision, dimension(0:Nyb+1):: uvol
Cnew    = C

Cnew(0:Nxb+1,0,1:Nz)       = A(0:Nxb+1,1:Nz)*C(0:Nxb+1,1,1:Nz)+ DC(0:Nxb+1,1:Nz)*C(0:Nxb+1,2,1:Nz)&
                           + B(0:Nxb+1,1:Nz)

do m = 1,Nz
Cnew(0:Nxb+1,Nyb+1,m)   = C(0:Nxb+1,Nyb,m)+ hx/8d0*(15d0*C(0:Nxb+1,Nyb,m)+10d0*C(0:Nxb+1,Nyb-1,m)+3d0*C(0:Nxb+1,Nyb-2,m)) &
                           * vvol(0:Nxb+1)
end do
if (jk.lt.0) then
   Cnew(0,0:Nyb+1,1:Nz)    = 2d0*Cup-C(1,0:Nyb+1,1:Nz) 
else
   Cnew(0,0:Nyb+1,1:Nz)    =  C(1,0:Nyb+1,1:Nz)
endif
do m = 1,Nz
Cnew(Nxb+1,0:Nyb+1,m)   = C(Nxb,0:Nyb+1,m) + hx/8d0*(15d0*C(Nxb,0:Nyb+1,m)+10d0*C(Nxb-1,0:Nyb+1,m)+3d0*C(Nxb-2,0:Nyb+1,m)) &
                           * uvol(0:Nyb+1)
end do


end subroutine ghost_cells
!-------------------------------------------------------------------------
subroutine ghost_cells_phig(phi_g)
double precision, dimension(0:Nxb+1,0:Nyb+1):: phi_g
phi_g(0:Nxb+1,Nyb+1)   = 0.5*(phi_g(0:Nxb+1,Nyb)-2.0*phi_g(0:Nxb+1,Nyb-1) + phi_g(0:Nxb+1,Nyb-2))/hy**2.*(-hy/2.)**2.0 &
                      + (-hy/2.0)*(-2.0*phi_g(0:Nxb+1,Nyb)-3.0*phi_g(0:Nxb+1,Nyb-1)+phi_g(0:Nxb+1,Nyb-2))/hy &
                      + (15.0*phi_g(0:Nxb+1,Nyb) -10.0*phi_g(0:Nxb+1,Nyb-1) + 3.0*phi_g(0:Nxb+1,Nyb-2))/8.0
phi_g(0:Nxb+1,0)       = 0.5*(phi_g(0:Nxb+1,1)-2.0*phi_g(0:Nxb+1,2) + phi_g(0:Nxb+1,3))/hy**2.*(-hy/2.)**2.0 &
                      + (-hy/2.0)*(-2.0*phi_g(0:Nxb+1,1)-3.0*phi_g(0:Nxb+1,2)+phi_g(0:Nxb+1,3))/hy &
                      + (15.0*phi_g(0:Nxb+1,1) -10.0*phi_g(0:Nxb+1,2) + 3.0*phi_g(0:Nxb+1,3))/8.0
phi_g(0,0:Nyb+1)       = 0.5*(phi_g(1,0:Nyb+1)-2.0*phi_g(2,0:Nyb+1) + phi_g(3,0:Nyb+1))/hx**2.*(-hx/2.)**2.0 &
                      + (-hx/2.0)*(-2.0*phi_g(1,0:Nyb+1)-3.0*phi_g(2,0:Nyb+1)+phi_g(3,0:Nyb+1))/hx &
                      + (15.0*phi_g(1,0:Nyb+1) -10.0*phi_g(2,0:Nyb+1) + 3.0*phi_g(3,0:Nyb+1))/8.0
phi_g(Nxb+1,0:Nyb+1)   = 0.5*(phi_g(Nxb,0:Nyb+1)-2.0*phi_g(Nxb-1,0:Nyb+1) + phi_g(Nxb-2,0:Nyb+1))/hx**2.*(-hx/2.)**2.0 &
                      + (-hx/2.0)*(-2.0*phi_g(Nxb,0:Nyb+1)-3.0*phi_g(Nxb-1,0:Nyb+1)+phi_g(Nxb-2,0:Nyb+1))/hx &
                      + (15.0*phi_g(Nxb,0:Nyb+1) -10.0*phi_g(Nxb-1,0:Nyb+1) + 3.0*phi_g(Nxb-2,0:Nyb+1))/8.0
end subroutine ghost_cells_phig
!--------------------------------------------------------------------------
subroutine Bump_type(btypeA)
integer:: btypeA
double precision:: A_space, B_space, L2_1, L2_2, L3_1, L3_2
double precision, dimension(Nxb, Nyb):: A_bump, B_bump, C_bump, Bump_tot
A_space = 21.33/100
B_space = 7.5/100;

L2_1 = L2+A_space
L2_2 = L2_1+B_space
L3_1 = L2_2+A_space
L3_2 = L3_1+B_space
do i = 1, Nxb
   do j = 1,Nyb
      if (btypeA > 1) then
         A_bump(i,j) = ((.5*tanh(25*(x(i)-L2))-5e-1*tanh(25*(x(i)-L2_1))))*(-.5*tanh(25.*(y(j) - .05))+0.5);
         B_bump(i,j) = ((.5*tanh(25*(x(i)-L2_2))-5e-1*tanh(25*(x(i)-L3_1))))*(-.5*tanh(25.*(y(j) - .05))+0.5);
         C_bump(i,j) = ((.5*tanh(25*(x(i)-L3_2))-5e-1*tanh(25d0*(x(i)-L3))))*(-.5*tanh(25.*(y(j) - .05))+0.5);
      elseif (btypeA ==1) then
         Bump_tot(i,j) = ((5d-1*tanh(5d0*(x(i)-L2))-5d-1*tanh(5d0*(x(i)-L3))))*(-5d-1*tanh(5d1*(y(j) - 0.05))+5d-1);
      end if
   end do
end do
if (btypeA==1) then
   Bumpcon = Bump_tot
elseif (btypeA ==2) then
   BumpCon = A_bump+B_bump+C_bump
elseif(btypeA==3) then
   BumpCon = A_bump+C_bump
end if



end subroutine Bump_type



!--------------------------------------------------------------------------
subroutine diffusion_subr(Nxb, Nyb, kx, ky, cup, BCX, BCxRHS, A, D, B, BCyRHS, Diff_TB, Diff_LR, G, C,ADV,RK,uvol, vvol)
integer:: Nxb, Nyb, i, j, jk
double precision:: kx,ky,cup
double precision, dimension(0:Nxb+1):: A,D,B
double precision, dimension(0:Nyb+1):: BCX
double precision, dimension(1:Nyb)::   BCxRHS, D22, RHS2y
double precision, dimension(1:Nxb)::   BCyRHS, D2, RHS2
double precision, dimension(Nxb,Nyb)::  C, ADV,Rk
double precision, dimension(0:Nxb+1, 0:Nyb+1):: G, Diff_LR, Gh, P, Ph, Diff_TB 
double precision, dimension(Nxb-1)::   DL, Du
double precision, dimension(Nyb-1)::   DL2, Du2
double precision, dimension(0:Nyb+1):: uvol
double precision, dimension(0:Nxb+1):: vvol
integer:: halfmatinfo, matinfo
!------------------


!----------------- 
!write(*,*) 'G=', G(1:Nxb,1)
P=G
Ph = 0d0
do j = 1, Nyb
     call N_mat(Diff_LR(1:Nxb,j), kx, Nxb, DL, D2, Du, BCX(j),0d0,1)
     call RHS(transpose(Diff_TB(1:Nxb,j:j+1)), ky, Nxb,&
     transpose(G(1:Nxb,j-1:j+1)), j, RHS2)
      RHS2(1) = RHS2(1)+BCxRHS(j)
     RHS2    = RHS2+ADV(1:Nxb,j)+RK(1:Nxb,j)
     call DGTSV(Nxb, 1, DL, D2, Du, RHS2, Nxb, halfmatinfo)
     Ph(1:Nxb,j) = RHS2

  end do
  call ghost_cells(Ph, A, B, D, Gh, Nxb, Nyb, 1,cup,jk,uvol,vvol)
  do i = 1,Nxb
     call N_mat(Diff_LR(i,1:Nyb), ky, Nyb, DL2, D22, Du2, A(i),D(i),1)

    call RHS(Diff_TB(i:i+1,1:Nyb), kx, Nyb, Gh(i-1:i+1,1:Nyb), 10, RHS2y)
     RHS2y(1)  = RHS2y(1) + BCyRHS(i)
     call DGTSV(Nyb,1,DL2,D22,Du2,RHS2y,Nyb,matinfo)
     C(i,1:Nyb) = RHS2y
  end do

  P(1:Nxb, 1:Nyb) = C
call ghost_cells(P, A, B, D, G, Nxb, Nyb,1, Cup,jk,uvol,vvol)

end subroutine diffusion_subr

!--------------------------------------------------------------------------
subroutine N_mat(Diff,k,N,Dl,Diag,Du,BC, BC2,writepar)

implicit none
integer:: N, writepar
double precision::k, BC, BC2
double precision, dimension(N):: Diag
double precision, dimension(N-1):: Dl, Du
double precision, dimension(N):: Diff

!Diagonal   

Diag(1)      = 1d+00 + k*(Diff(2)+Diff(1)*(1d0-BC))
Diag(2:N-1)  = 1d+00 + k*(Diff(2:N-1)+Diff(3:N))
Diag(N)      = 1d+00 + k*Diff(N)

!Lower Diagonal

Dl = -k*Diff(2:N)


!Upper Diagonal 

Du    = -k*Diff(2:N)
Du(1) = -k*(Diff(2)+Diff(1)*BC2)


end subroutine N_mat
!-----------------------------------------------------------------------------------------
subroutine RHS(D,k,N,C,i, RHSvec)

integer:: N, i, m
double precision,dimension(2,N):: D
double precision, dimension(3,N):: C
double precision:: k
double precision, dimension(N)::RHSvec

RHSvec = C(2,1:N) + k*(D(2,1:N)*(C(3,1:N)-C(2,1:N)) - D(1,1:N)*(C(2,1:N)-C(1,1:N)))

end subroutine RHS

!-------------------------------------------------------------------------        


subroutine ghost_cell_vec(C,Cnew,A,B,DC,Nxb,Nyb,Nz,Cup,jk,uvol,vvol1)
integer:: Nxb, Nyb, Nz, jk, m
double precision, dimension(0:Nxb+1,0:Nyb+1,1:Nz) :: C
double precision, dimension(-1:Nxb+2, -1:Nyb+2,1:Nz):: Cnew
double precision:: Cup
double precision, dimension(0:Nxb+1,1:Nz):: A, B, DC
double precision, dimension(0:Nxb+1):: vvol1
double precision, dimension(0:Nyb+1):: uvol


Cnew(1:Nxb, 1:Nyb,1:Nz)    = C(1:Nxb,1:Nyb,1:Nz)
Cnew(0:Nxb+1,0,1:Nz)      = A(0:Nxb+1,1:Nz)*C(0:Nxb+1,1,1:Nz)+ &
                           DC(0:Nxb+1,1:Nz)*C(0:Nxb+1,2,1:Nz) + B(0:Nxb+1,1:Nz)

if (Nz==1) then
   Cnew(0:Nxb+1,Nyb+1,1)   = C(0:Nxb+1,Nyb,1)+ hx/8d0*(15d0*C(0:Nxb+1,Nyb,1)+1d1*C(0:Nxb+1,Nyb-1,1)+3d0*C(0:Nxb+1,Nyb-2,1)) &
                           * vvol1(0:Nxb+1)
else
   do m = 1,Nz
      Cnew(0:Nxb+1,Nyb+1,m)   = C(0:Nxb+1,Nyb,m)+ hx/8d0*(15d0*C(0:Nxb+1,Nyb,m)+10d0*C(0:Nxb+1,Nyb-1,m)+3d0*C(0:Nxb+1,Nyb-2,m)) &
                           * vvol1(0:Nxb+1)
end do
end if

if (JK.lt.0) then
   Cnew(0,0:Nyb+1,1:Nz)      = 2d0*Cup-C(1,0:Nyb+1,1:Nz)
else
   Cnew(0,0:Nyb+1,1:Nz)      = C(1,0:Nyb+1,1:Nz)
endif

do m = 1,Nz
Cnew(Nxb+1,0:Nyb+1,m)   = C(Nxb,0:Nyb+1,m) + hx/8d0*(15d0*C(Nxb,0:Nyb+1,m)+10d0*C(Nxb-1,0:Nyb+1,m)+3d0*C(Nxb-2,0:Nyb+1,m)) &
                          *uvol(0:Nyb+1)
end do


Cnew(-1:Nxb+2,-1,1:Nz)    = 0d0
Cnew(-1:Nxb+2,Nyb+2,1:Nz)  = 0d0

Cnew(0,-1,1:Nz)   = Cup
Cnew(0,Nyb+2,1:Nz) = Cup
Cnew(-1,-1:Nyb+2,1:Nz)    = Cup

Cnew(Nxb+2, 0:Nyb+1,1:Nz) = C(Nxb,0:Nyb+1,1:Nz)


Cnew(-1,0,1:Nz)=Cup
Cnew(Nxb+2,0,1:Nz) = 0d0                                                      
Cnew(-1,Nyb+1,1:Nz) = 0d0
Cnew(Nxb+2,Nyb+1,1:Nz) =  C(Nxb,Nyb+1,1:Nz)
Cnew(Nxb+1,-1,1:Nz)   =  C(Nxb,0,1:Nz)
Cnew(Nxb+1,Nyb+2,1:Nz) = C(Nxb,Nyb+1,1:Nz) !0d0



end subroutine ghost_cell_vec
!--------------------------------------------------------------------------
subroutine advection(C, Ulr, Utb, Vlr,Vtb, dt, hx, hy, Nxb, Nyb,ADV)                   
integer:: Nxb, Nyb
double precision, dimension(-1:Nxb+2, -1:Nyb+2):: F, G
double precision, dimension(-1:Nxb+2, -1:Nyb+2):: C, Cnew, Ulr, Utb, Vlr, Vtb
double precision:: hx, hy,dt, Uij, Vij, Rcij, S, thetaa, phi
integer:: i,j, casenumber, I2, I3, J2, J3
double precision, dimension(1:Nxb,1:Nyb):: ADV
double precision, dimension(-1:Nyb+2):: Rcvec

F = 0.0
G = 0.0
Cnew = 0.0
do j = 0, Nyb+1
   do i = 0, Nxb+1
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

      G(I2,J2) = G(I2,J2) - 5d-1*dt/hy*Uij*Vij*Rcij


      !---------------------------------------                                    
      !Limited R:                                                                 
      if (abs(RCij)==0.0) then
         thetaa = 0.0
      else
         if (Uij > 0.0) then
            I3 = i-1
         else
            I3 = i+1
         end if

         thetaa = (C(I3,j) - C(I3-1,j))/RCij
         end if
        casenumber = 4
        phi=philim(thetaa, casenumber)


      RCij  = RCij*phi
      S     = 0.5*abs(Uij)*(1.0-dt/hx*abs(Uij))*RCij
      F(i,j)= F(i,j) + S

      !---------------------------------------                                    
      G(i,J2)   = G(i,J2)   + dt/hy*Vij*S
      G(i-1,J2) = G(i-1,J2) - dt/hy*Vij*S

   end do
end do


do i = 0, Nxb+1

   do j = 0, Nyb+1
      Uij  = Utb(i,j)
      Vij  = Vtb(i,j)

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

      F(I2,J2) = F(I2,J2) - 5d-1*dt/hx*Vij*Uij*Rcij

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
         thetaa = 0.0
      else
         if (Vij > 0.0) then
            J3 = j-1
            else
            J3 = j+1
        end if

         thetaa = (C(i,J3) - C(i,J3-1))/RCij
     end if

      casenumber = 4

      phi=philim(thetaa, casenumber)
      RCij   = RCij*phi

      S      = 5d-1*abs(Vij)*(1d0-dt/hy*abs(Vij))*RCij

      G(i,j) = G(i,j) + S

      !---------------------------------------                                    
      F(I2,j)   = F(I2,j)   + dt/hx*Utb(I2,j)*S
      F(I2,j-1) = F(I2,j-1) - dt/hx*Utb(I2,j-1)*S

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
!------------------------------------------------------------------               
function philim(thetaa, casenumber)

integer:: casenumber
double precision:: thetaa, philim 

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
   philim = max(0.0, min((1.0+thetaa)/2.0,2.0,2.0*thetaa))
   !   return                                                                        
!end select                                                                       

end function philim
!----------------------------------------------------------------------
subroutine RungaKutta(E,C,kappaf,kappafs,Nxb,Nyb,dt,Rkf)
integer:: Nxb, Nyb
double precision:: dt, kappaf,kappafs
double precision, dimension(0:Nxb+1,0:Nyb+1):: E, C
double precision, dimension(1:Nxb,1:Nyb):: Rkf,S1

S1 = C(1:Nxb,1:Nyb)-(5d-1*dt*kappaf*E(1:Nxb,1:Nyb)*C(1:Nxb,1:Nyb))/(phi_f(1:Nxb,1:Nyb)*kappafs+C(1:Nxb,1:Nyb))

Rkf = -dt*kappaf*E(1:Nxb,1:Nyb)*S1(1:Nxb,1:Nyb)/(phi_f(1:Nxb,1:Nyb)*kappafs+S1(1:Nxb,1:Nyb))

end subroutine RungaKutta
!------------------------------------------------------------------
subroutine beam_warming(Fz,Fpw,type)

integer:: type, m
double precision, dimension(Nxb, Nyb,1:Nz):: Fz, Fpw,fun_o
double precision, dimension(1:Nz):: alpha_w, beta_w, gamma_w
alpha_w(2:Nz) = hz(2:Nz)/(hz(1:Nz-1)*(hz(2:Nz)+hz(1:Nz-1)))
alpha_w(1)    = alpha_w(2)
beta_w(2:Nz)   = -(hz(2:Nz)+ hz(1:Nz-1))/(hz(2:Nz)*hz(1:Nz-1))
beta_w(1)      = beta_w(2)        
gamma_w        = -(alpha_w+beta_w)
fun_o = 0.0
if (type==1) then
   call function_W(fun_o)
   
end if
if (type==2) then
   call function_V(fun_o)
end if
do m = 3, Nz                                                                                   
   Fz(1:Nxb,1:Nyb, m)  = alpha_w(m)*fun_o(1:Nxb,1:Nyb,m-2) + beta_w(m)*fun_o(1:Nxb,1:Nyb,m-1) &
                       + gamma_w(m)*fun_o(1:Nxb,1:Nyb,m)
end do                                                                                         
Fz(1:Nxb,1:Nyb,2)      = beta_w(2) * fun_o(1:Nxb,1:Nyb,1) + gamma_w(2)*fun_o(1:Nxb,1:Nyb,2)               
Fz(1:Nxb, 1:Nyb,1)     = fun_o(1:Nxb,1:Nyb,1)/hz(1)


if (type==1) then
   do i = 1,Nz
      Fpw(1:Nxb,1:Nyb,i) = klMat/phi_f(1:Nxb,1:Nyb)*Wold(1:Nxb,1:Nyb,i) + kbMAT/(2d0*phi_f(1:Nxb,1:Nyb)**2)*(Wold(1:Nxb,1:Nyb,i)**2.&
                         + 2d0*z(i)*Wold(1:Nxb,1:Nyb,i)*Rold(1:Nxb,1:Nyb) &
                         + z(i)**2.*Rold(1:Nxb,1:Nyb)**2.) - kbMAT/(2d0*phi_f(1:Nxb,1:Nyb)**2)*z(i)*(Rold(1:Nxb,1:Nyb)**2.+ chi*Rgold(1:Nxb,1:Nyb)**2.)

   end do
end if


if (type==2) then
   do m = 3, Nz
      Fpw(1:Nxb,1:Nyb, m)  = alpha_w(m)*W(1:Nxb,1:Nyb,m-2) + beta_w(m)*W(1:Nxb,1:Nyb,m-1) + gamma_w(m)*W(1:Nxb,1:Nyb,m)
   end do
   Fpw(1:Nxb,1:Nyb,2)      = beta_w(2) * W(1:Nxb,1:Nyb,1) + gamma_w(2)*W(1:Nxb,1:Nyb,2)
   Fpw(1:Nxb,1:Nyb,1)      = W(1:Nxb,1:Nyb,1)/hz(1)
end if
end subroutine beam_warming
!--------------------------------------------------------------------------------------
subroutine function_W(fji)

integer:: midx
double precision, dimension(Nxb, Nyb, 1:Nz)::fji
double precision, dimension(Nxb,Nyb):: Ws,Rshort
double precision:: zi
fji = 0
Rshort = Rold(1:Nxb,1:Nyb)

do midx  = 1,Nz
   zi = z(midx)
   Ws = Wold(1:Nxb,1:Nyb,midx)
   fji(1:Nxb,1:Nyb,midx) = -1d0*(klMAT/(2d0*phi_f(1:Nxb,1:Nyb))*Ws**2. - kbMAT/(2d0*phi_f(1:Nxb,1:Nyb)**2)*zi*Ws*(Rshort**2.-chi*Rgold(1:Nxb,1:Nyb)**2.)&
                         - kbMAT/(2d0*phi_f(1:Nxb,1:Nyb)**2)*zi**2.*(Rshort**3.-chi*Rshort*Rgold(1:Nxb,1:Nyb)**2.) &
                         + kbMAT/(2d0*phi_f(1:Nxb,1:Nyb)**2)*(1d0/3d0*Ws**3.+zi*Ws**2.*Rshort +zi**2.*Rshort**2.*Ws))


end do



end subroutine function_W

!---------------------------------------------------------------------------                            
subroutine function_V(fji)

integer:: midx
double precision, dimension(Nxb, Nyb, 1:Nz)::fji
double precision, dimension(Nxb,Nyb):: Ws, Rshort, Vs

fji = 0

Rshort = (R(1:Nxb,1:Nyb)+Rold(1:Nxb,1:Nyb))/2.0
do midx  = 1,Nz
   Ws = W(1:Nxb,1:Nyb,midx)
   Vs = Vold(1:Nxb,1:Nyb,midx)

   fji(1:Nxb,1:Nyb,midx) = -1d0*(((klMAT/phi_f(1:Nxb,1:Nyb))*Ws+kbMAT/(2d0*phi_f(1:Nxb,1:Nyb)**2)*(Ws+z(midx)*Rshort)**2. &
                           - kbMAT/(2d0*phi_f(1:Nxb,1:Nyb)**2)*z(midx)*(Rshort**2.-chi*Rg(1:Nxb,1:Nyb)**2.))*Vs)



end do

end subroutine function_V
!--------------------------------------------------------------------------------------
subroutine RXN_R(rxnr)


double precision, dimension(0:Nxb+1,0:Nyb+1)::rxnr, klmat2, kbmat2
klmat2(1:Nxb,1:Nyb) = klmat
kbmat2(1:Nxb,1:Nyb) = kbmat

klmat2(1:Nxb,0) = klmat2(1:Nxb,1)
kbmat2(1:Nxb,0) = kbmat2(1:Nxb,1)

klmat2(1:Nxb,Nyb) = klmat2(1:Nxb,Nyb)
kbmat2(1:Nxb,Nyb) = kbmat2(1:Nxb,Nyb)

klmat2(0,0:Nyb+1) = kl
kbmat2(0,0:Nyb+1) = kb

kbmat2(Nxb+1,0:Nyb+1) = 0
klmat2(Nxb+1,0:Nyb+1) = 0 


rxnr = dt*(-klmat2/phi_f(0:Nxb+1,0:Nyb+1)*(Rold**2. - chi*Rgold**2.) - kbmat2/(2d0*phi_f(0:Nxb+1,0:Nyb+1)**2)*(Rold**3.- chi*(3d0*Rsold*Rgold**2.+Rgold**3.)) + 2d0*S10c)

end subroutine RXN_R
!------------------------------------------------------------------------------------------
subroutine diffusion_subr_r(Nxb, Nyb, kx, ky, cup, BCX, BCxRHS, A, D, B, BCyRHS, Diff_TB,Diff_LR, Rg, G, Rgold, RXNr, ADVr,dt,jk,uvol, vvol)
integer:: Nxb, Nyb, i, j,jk
double precision:: kx,ky,cup,dt
double precision, dimension(0:Nxb+1):: A,D,B
double precision, dimension(0:Nyb+1):: BCX
double precision, dimension(1:Nyb)::   BCxRHS, D22, RHS2y,RHSg2
double precision, dimension(1:Nxb)::   BCyRHS, D2, RHS2, RHSg
double precision, dimension(Nxb,Nyb)::  ADVr 
double precision, dimension(0:Nxb+1, 0:Nyb+1):: Rg, Rgold, G, Diff2,Diff_TB, DIff_LR, Gh, P, Ph, RXNr 
double precision, dimension(Nxb-1)::   DL, Du
double precision, dimension(Nyb-1)::   DL2, Du2
double precision, dimension(0:Nyb+1)::uvol
double precision, dimension(0:Nxb+1):: vvol
integer:: halfmatinfo, matinfo
P=0d0

G(1:Nxb,1:Nyb) = G(1:Nxb,1:Nyb)!+ADVr(1:Nxb,1:Nyb) + RXNr(1:Nxb,1:Nyb)


if (jk.le.1) then

   do j = 1, Nyb
      call N_mat(Diff_TB(1:Nxb,j), kx, Nxb, DL, D2, Du, BCX(j),0d0,1)
      call RHS(transpose(Diff_LR(1:Nxb,j:j+1)), ky, Nxb,&
           transpose(G(1:Nxb,j-1:j+1)), 10, RHS2)

      ! Change the second Rgold to Rg if you want to use a time average
      call RHS_old(Diff_TB(1:Nxb+1,j:j+1),Diff_LR(1:Nxb,j:j+1), &
                   kx,ky,Nxb,-Rgold(0:Nxb+1,j-1:j+1),-Rg(0:Nxb+1,j-1:j+1),RHSg)

      RHS2(1)     = RHS2(1)+BCxRHS(j)
      RHS2        = RHS2 + RHSg +ADVr(1:Nxb,j) + RXNr(1:Nxb,j)
     
      call DGTSV(Nxb, 1, DL, D2, Du, RHS2, Nxb, halfmatinfo)

      Ph(1:Nxb,j) = RHS2
   end do
   call ghost_cells(Ph, A, B, D, Gh, Nxb, Nyb, 1,cup,jk,uvol,vvol)
   
   do i = 1,Nxb
      call N_mat(Diff_TB(i,1:Nyb), ky, Nyb, DL2, D22, Du2, A(i),D(i),0)

      call RHS(Diff_LR(i:i+1,1:Nyb), kx, Nyb, Gh(i-1:i+1,1:Nyb), i, RHS2y)
      RHS2y(1)  = RHS2y(1) + BCyRHS(i)
      call DGTSV(Nyb,1,DL2,D22,Du2,RHS2y,Nyb,matinfo)
      P(i,1:Nyb) = RHS2y
   end do
else
   do j = 1,Nyb
      call N_mat(Diff_TB(1:Nxb,j), kx, Nxb, DL, D2, DU, BCX(j),0d0,0)
      RHS2 = (kx*Diff_LR(1:Nxb,1)*G(0:Nxb-1,j) &
           + (1d0-2d0*kx*Diff_LR(1:Nxb,1))*G(1:Nxb,j) + kx*Diff_LR(1:Nxb,j)*G(2:Nxb+1,j))
      
      !This is going to need to be rethought through with
      call RHS_old(Diff_TB(1:Nxb+1,j:j+1),Diff_LR(1:nxb+1,j:j+1), &
                kx,ky,Nxb,-Rg(0:Nxb+1,j-1:j+1),-Rgold(0:Nxb+1,j-1:j+1),RHSg)

      RHS2 = RHS2 + ADVr(1:Nxb,j) + RHSg + RXNr(1:Nxb,j)
      call DGTSV(Nxb, 1, DL, D2, Du, RHS2, Nxb, halfmatinfo)
      P(1:Nxb, j) = RHS2
   end do
endif
call ghost_cells(P, A, B, D, G, Nxb, Nyb, 1, Cup,jk,uvol,vvol)

end subroutine diffusion_subr_R
!------------------------------------------------------------------------------------------
subroutine RHS_old(D_TB, D_LR,kx,ky,Nxb,cold,cnew,RHSC)
integer:: Nxb
double precision, dimension(0:Nxb+1,3):: cold, cnew, C
double precision, dimension(Nxb):: RHSC
double precision, dimension(Nxb+1,2):: D_TB, D_LR
double precision:: kx, ky


C    = (cold+cnew)/2d0

!Ky is already divided by 2 so multiply it by 2
RHSC = 2d0*ky*(D_TB(1:Nxb,2)*(C(1:Nxb,3)-C(1:Nxb,2)) - D_TB(1:Nxb,1)*(C(1:Nxb,2)-C(1:Nxb,1))) &
     + 2d0*kx*(D_LR(2:Nxb+1,2)*(C(2:Nxb+1,2)-C(1:Nxb,2))-D_LR(1:Nxb,2)*(C(1:Nxb,2)-C(0:Nxb-1,2)))


end subroutine RHS_old
!------------------------------------------------------------------------------------------
subroutine RHS_old2(D_TB,D_LR,kx,ky,Nyb,cold,cnew,RHSC)
integer:: Nyb
double precision, dimension(3,0:Nyb+1):: cold, cnew, C
double precision, dimension(Nyb):: RHSC
double precision, dimension(2,Nyb+1):: D_TB, D_LR
double precision:: kx, ky


C    = (cold+cnew)/2d0

!Ky is already divided by 2 so multiply it by 2
RHSC = 2d0*kx*(D_LR(2,1:Nyb)  * (C(3,1:Nyb)-C(2,1:Nyb))  - D_LR(1,1:Nyb)*(C(2,1:Nyb)-C(1,1:Nyb))) &
     + 2d0*ky*(D_TB(2,2:Nyb+1)* (C(2,2:Nyb+1)-C(2,1:Nyb))- D_TB(2,1:Nyb)*(C(2,1:Nyb)-C(2,0:Nyb-1)))

end subroutine RHS_old2
!--------------------------------------------------------------------------------------
subroutine RXNW(Rxn_zW)
integer:: m
double precision, dimension(1:Nxb,1:Nyb,1:Nz):: Rxn_zW

do m = 1, Nz

Rxn_zW(1:Nxb,1:Nyb,m) = dt*(-klmat/phi_f(1:Nxb,1:Nyb)*chi*z(m)*Rgavg(1:Nxb,1:Nyb)**2. + kbmat/(2d0*phi_f(1:Nxb,1:Nyb)**2)*z(m)**2.*Ravg(1:Nxb,1:Nyb)**3.)  &
                      - dt*kbmat/(2d0*phi_f(1:Nxb,1:Nyb)**2)*(-z(m)*Ravg(1:Nxb,1:Nyb)**3.+3d0*z(m)*chi*Rsavg(1:Nxb,1:Nyb)*Rgavg(1:Nxb,1:Nyb)**2.&
                      + z(m)*chi*Rgavg(1:Nxb,1:Nyb)**3.) 

end do

end subroutine RXNW
!---------------------------------------------------------------------------
subroutine diffusion_subr_W(Nxb, Nyb, Nz,kx, ky, cup, BCX, BCxRHS, A,B, D, BCyRHS, Diff_TB, Diff_LR, G,z,&
                           Wold1, ADVw, RXN_zW,ADVzW,Wnewg,dt,jk,uvol, vvol)
integer:: Nxb, Nyb, Nz, i, j, m
double precision:: kx,ky,cup,dt
double precision, dimension(0:Nxb+1,1:Nz):: A, B, D
double precision, dimension(0:Nyb+1):: BCX
double precision, dimension(1:Nyb)::   BCxRHS, D22, RHS2y, RHSzW2
double precision, dimension(1:Nxb)::   D2, RHS2, RHSzW
double precision, dimension(1:Nxb,1:Nz):: BCyRHS
double precision, dimension(Nxb,Nyb)::  C
double precision, dimension(0:Nxb+1, 0:Nyb+1):: Diff_LR, Diff_TB, Gh,Ph, Diff2 
double precision, dimension(0:Nxb+1,0:Nyb+1,1:Nz):: G,P,Wnewg
double precision, dimension(1:Nxb,1:Nyb,1:Nz):: ADVw, RXN_zW, ADVzW
double precision, dimension(Nxb-1)::   DL, Du
double precision, dimension(Nyb-1)::   DL2, Du2
double precision, dimension(0:Nxb+1,0:Nyb+1)::zWnew, Wold1, zWold
double precision, dimension(Nz):: z
double precision, dimension(0:Nyb+1):: uvol
double precision, dimension(0:Nxb+1):: vvol
integer:: halfmatinfo, matinfo,jk

Ph = 0d0
zWnew = 0d0
zWold = 0d0


if (jk.le.1) then
   do m = 1,Nz-1
      zWold = -z(m)*Wold1
      zWnew = -z(m)*Wnewg(0:Nxb+1,0:Nyb+1,Nz)

      do j = 1, Nyb
         call N_mat(Diff_LR(1:Nxb,j), kx, Nxb, DL, D2, Du, BCX(j),0d0,0)
         
         call RHS(transpose(Diff_TB(1:Nxb,j:j+1)), ky, Nxb,&
              transpose(G(1:Nxb,j-1:j+1,m)), 10, RHS2)
       !This is going to need to be rethought through with 
         call RHS_old(Diff_TB(1:Nxb+1,j:j+1),Diff_LR(1:Nxb,j:j+1), kx,ky,Nxb,zWold(0:Nxb+1,j-1:j+1),zWnew(0:Nxb+1,j-1:j+1),RHSzW)
         RHS2 = RHS2 + RHSzW  + ADVw(1:Nxb,j,m)- dt*ADVzW(1:Nxb,j,m) + RXN_zW(1:Nxb,j,m)
         
         call DGTSV(Nxb, 1, DL, D2, Du, RHS2, Nxb, halfmatinfo)
         Ph(1:Nxb,j) = RHS2

      end do

      call ghost_cells(Ph, A(0:Nxb+1,m),B(0:Nxb+1,m), D(0:Nxb+1,m), Gh, Nxb, Nyb, 1, cup,jk,uvol,vvol)

      
      do i = 1,Nxb
         call N_mat(Diff_TB(i,1:Nyb), ky, Nyb, DL2, D22, Du2, A(i,m),0d0,0)

         call RHS(Diff_LR(i:i+1,1:Nyb), kx, Nyb, Gh(i-1:i+1,1:Nyb), i, RHS2y)
        
         RHS2y(1)  = RHS2y(1) + BCyRHS(i,m)

         call DGTSV(Nyb,1,DL2,D22,Du2,RHS2y,Nyb,matinfo)
         C(i,1:Nyb) = RHS2y
      end do

      P(1:Nxb, 1:Nyb,m) = C
   end do
else
   do m = 1, Nz-1
      zWold = -z(m)*Wold1
      zWnew = -z(m)*Wnew(0:Nxb+1,0:Nyb+1,Nz)
      do j = 1, Nyb
         call N_mat(Diff_LR(1:Nxb,j), kx, Nxb, DL, D2, Du, BCX(j),0d0,0)

         RHS2 = (kx*Diff_LR(1:Nxb,j)*G(0:Nxb-1,j,m) + (1d0-2d0*kx*Diff_LR(1:Nxb,j))*G(1:Nxb,j,m) + kx*Diff_LR(1:Nxb,j)*G(2:Nxb+1,j,m))
         !This is going to need to be rethought through with
         call RHS_old(Diff_TB(1:Nxb+1,j:j+1),Diff_LR(1:Nxb+1,j:j+1), &
                     kx,ky,Nxb,zWold(0:Nxb+1,j-1:j+1),zWnew(0:Nxb+1,j-1:j+1),RHSzW)


         RHS2 = RHS2 + ADVw(1:Nxb,j,m) + RHSzW + RXN_zW(1:Nxb,j,m) - dt*ADVzW(1:Nxb,j,m)
         call DGTSV(Nxb, 1, DL, D2, Du, RHS2, Nxb, halfmatinfo)
         C(1:Nxb, j ) = RHS2
      end do
      P(1:Nxb,1:Nyb,m) = C(1:Nxb,1:Nyb)
   end do
endif
P(1:Nxb,1:Nyb,Nz) = Wnew(1:Nxb,1:Nyb,Nz)
call ghost_cells(P, A,B,D, G, Nxb, Nyb, Nz,Cup,jk,uvol,vvol)

end subroutine diffusion_subr_W
!--------------------------------------------------------------------------------------
subroutine Lax_Wendroff(LW_z)
double precision, dimension(1:Nxb,1:Nyb,1:Nz):: LW_z
integer:: i
double precision, dimension(Nxb,Nyb,Nz):: Wh, Fpw, fLw
double precision, dimension(Nxb,Nyb):: flwpi, flwmi,flwi, Fpwpi, Fpwi
call function_W(fLw)

Wh(1:Nxb,1:Nyb,1)    = Wold(1:Nxb,1:Nyb,1)/2d0
Wh(1:Nxb,1:Nyb,2:Nz) = (Wold(1:Nxb,1:Nyb,2:Nz) + Wold(1:Nxb,1:Nyb,1:Nz-1))/2d0
call function_WLW(Wh,Fpw)
do i = 2,Nz-1
   flwpi = flw(1:Nxb,1:Nyb,i+1)
   flwmi = flw(1:Nxb,1:Nyb,i-1)
   flwi  = flw(1:Nxb,1:Nyb,i)
   Fpwpi = Fpw(1:Nxb,1:Nyb,i+1)
   Fpwi  = Fpw(1:Nxb,1:Nyb,i)
   LW_z(1:Nxb,1:Nyb,i) = (flwpi-flwmi)/(hz(i+1)+hz(i)) - dt/(hz(i+1)+hz(i))*(Fpwpi*(Flwpi-Flwi)/hz(i+1)-Fpwi*(Flwi-Flwmi)/hz(i))
end do
LW_z(1:Nxb,1:Nyb,1) = (flw(1:Nxb,1:Nyb,2)-0d0)/(hz(2)+hz(1)) &
                    - dt/(hz(2)+hz(1))*(Fpw(1:Nxb,1:Nyb,2)*(flw(1:Nxb,1:Nyb,2)-flw(1:Nxb,1:Nyb,1))/hz(2) &
                    - Fpw(1:Nxb,1:Nyb,1)*(Flw(1:Nxb,1:Nyb,1)-0d0)/hz(1))
LW_z(1:Nxb,1:Nyb,Nz) = 0d0
end subroutine Lax_Wendroff
!--------------------------------------------------------------------------------------  
subroutine function_WLW(Whalf,Fphalf)

integer:: midx
double precision, dimension(Nxb, Nyb, 1:Nz)::Fphalf, Whalf
double precision, dimension(Nxb,Nyb):: Ws,Rshort


Rshort = Rold(1:Nxb,1:Nyb)
do midx  = 1,Nz
   Ws = Whalf(1:Nxb,1:Nyb,midx)
   Fphalf(1:Nxb,1:Nyb,midx) = -1d0*(klmat/phi_f(1:Nxb,1:Nyb)*Ws + kbmat/(2d0*phi_f(1:Nxb,1:Nyb)**2)*(Ws**2. + 2d0*z(midx)*Ws*Rshort &
                            + z(midx)**2.*Rshort**2.) - kbmat/(2d0*phi_f(1:Nxb,1:Nyb)**2)*z(midx)*(Rshort**2.- chi*Rgold(1:Nxb,1:Nyb)**2.))
end do



end subroutine function_WLW

!--------------------------------------------------------------------------------------
subroutine Lax_WendroffV(LWV_z)
double precision, dimension(1:Nxb,1:Nyb,1:Nz):: LWV_z
integer:: kidx
double precision, dimension(Nxb,Nyb,Nz):: Vhalf, FpV,flv
double precision, dimension(Nxb,Nyb):: flvpi,flpi,flvmi,flvi,Fpvpi,Fpvi
Vhalf(1:Nxb,1:Nyb,1)    = Vold(1:Nxb,1:Nyb,1)/2d0
Vhalf(1:Nxb,1:Nyb,2:Nz) = (Vold(1:Nxb,1:Nyb,2:Nz) + Vold(1:nxb,1:Nyb,1:Nz-1))/2d0
call function_V(flV)

call function_VLW(Vhalf, Fpv)


do kidx = 2,Nz-1
   flvpi = flv(1:Nxb,1:Nyb,kidx+1)
   flvmi = flv(1:Nxb,1:Nyb,kidx-1)
   flvi  = flv(1:Nxb,1:Nyb,kidx)
   Fpvpi = Fpv(1:Nxb,1:Nyb,kidx+1)
   Fpvi  = Fpv(1:Nxb,1:Nyb,kidx)
   LWV_z(1:Nxb,1:Nyb,kidx) = (flvpi-flvmi)/(hz(kidx+1)+hz(kidx)) -dt/(hz(kidx+1)+hz(kidx))*(Fpvpi*(Flvpi-Flvi)/hz(kidx+1)-Fpvi*(Flvi-Flvmi)/hz(kidx))
end do
LWV_z(1:Nxb,1:Nyb,1) = (flv(1:Nxb,1:Nyb,2)-0d0)/(hz(2)+hz(1)) &
                     - dt/(hz(2)+hz(1))*(Fpv(1:Nxb,1:Nyb,2)*(flv(1:Nxb,1:Nyb,2)-flv(1:Nxb,1:Nyb,1))/hz(2) &
                     - Fpv(1:Nxb,1:Nyb,1)*(Flv(1:Nxb,1:Nyb,1)-0)/hz(1))
LWV_z(1:Nxb,1:Nyb,Nz) = 0d0





!write(*,*) 'fph=', fphalf(1:Nx,1,1)
!write(*,*) 'fmh=', fmhalf(1:Nx,1,1)

end subroutine Lax_WendroffV
!--------------------------------------------------------------------------------------
subroutine function_VLW(Vhalf,fhalf)

integer:: midx
double precision, dimension(Nxb, Nyb, 1:Nz)::fhalf, Vhalf
double precision, dimension(Nxb,Nyb):: Vs,Rshort, Ws

fhalf = 0

do midx  = 1,Nz
   Ws = W(1:Nxb,1:Nyb,midx)
   Rshort = R(1:Nxb,1:Nyb)
   fhalf(1:Nxb,1:Nyb,midx) = -1d0*((klmat/phi_f(1:Nxb,1:Nyb)*Ws+kbmat/(2d0*phi_f(1:Nxb,1:Nyb)**2)*(Ws+ z(midx)*Rshort)**2. &
                            - kbmat/(2d0*phi_f(1:Nxb,1:Nyb)**2)*z(midx)*(Rshort**2.-chi*Rg(1:Nxb,1:Nyb)**2.)))
 
end do


end subroutine function_VLW
!--------------------------------------------------------------------------------------

subroutine RXN_V(RXNv)
integer:: midx
double precision, dimension(1:Nxb,1:Nyb,1:Nz):: RXNv
double precision, dimension(1:Nxb,1:Nyb):: Ravg
Ravg = (Rold(1:Nxb,1:Nyb) + R(1:Nxb,1:Nyb))/2d0

do midx = 1, Nz
   RXNv(1:Nxb,1:Nyb,midx) =dt*(-kbmat/phi_f(1:Nxb,1:Nyb)**2*(Wnew(1:Nxb,1:Nyb,midx)+z(midx)*Ravg)**2.*(Wz(1:Nxb,1:Nyb,midx)+Ravg)) &
                       + dt*2d0*S10c(1:Nxb,1:Nyb)*z(midx)
end do


end subroutine RXN_V
!------------------------------------------------------------------------------------------ 
subroutine diffusion_subr_V(Nxb, Nyb, Nz,kx, ky, cup, BCX, BCxRHS, A,B, D, BCyRHS, Diff_TB, Diff_LR, G, P, &
                           ADVV, RXNV,ADVzV,dt,jk,uvol, vvol)

integer:: Nxb, Nyb, Nz, i, j, m
double precision:: kx,ky,cup,dt 
double precision, dimension(0:Nxb+1,1:Nz):: A, B, D
double precision, dimension(0:Nyb+1):: BCX
double precision, dimension(1:Nyb)::   BCxRHS, D22, RHS2y, RHSzW2
double precision, dimension(1:Nxb)::   D2, RHS2, RHSzW
double precision, dimension(1:Nxb,1:Nz):: BCyRHS
double precision, dimension(Nxb,Nyb)::  C
double precision, dimension(0:Nxb+1, 0:Nyb+1):: Diff_TB, Diff_LR, Gh,Ph 
double precision, dimension(0:Nxb+1,0:Nyb+1,1:Nz):: G,P
double precision, dimension(1:Nxb,1:Nyb,1:Nz):: ADVV, RXNV, ADVzV
double precision, dimension(Nxb-1)::   DL, Du
double precision, dimension(Nyb-1)::   DL2, Du2
double precision, dimension(0:Nyb+1):: uvol
double precision, dimension(0:Nxb+1):: vvol
integer:: halfmatinfo, matinfo,jk

Ph = 0d0

if (jk.le.1) then
do m = 1,Nz


 do j = 1, Nyb
      call N_mat(Diff_LR(1:Nxb,j), kx, Nxb, DL, D2, Du, BCX(j),0d0,0)
      call RHS(transpose(Diff_TB(1:Nxb,j:j+1)), ky, Nxb,&
           transpose(G(1:Nxb,j-1:j+1,m)), 10, RHS2)
      !This is going to need to be rethought through with                                   
      RHS2(1) = RHS2(1)+BCxRHS(j)
      RHS2 = RHS2 + ADVV(1:Nxb,j,m) + RXNV(1:Nxb,j,m)- dt*ADVzV(1:Nxb,j,m)                                                                                         
      call DGTSV(Nxb, 1, DL, D2, Du, RHS2, Nxb, halfmatinfo)
      PH(1:Nxb,j) = RHS2
   end do
 call ghost_cells(Ph, A(0:Nxb+1,m),B(0:Nxb+1,m), D(0:Nxb+1,m), Gh, Nxb, Nyb, 1, cup,jk,uvol,vvol)
 do i = 1,Nxb
      call N_mat(Diff_TB(i,1:Nyb), ky, Nyb, DL2, D22, Du2, A(i,m),0d0,0)
      call RHS(Diff_LR(i:i+1,1:Nyb), kx, Nyb, Gh(i-1:i+1,1:Nyb), i, RHS2y)
      RHS2y(1)  = RHS2y(1) + BCyRHS(i,m)
      !RHS2y = RHS2y  + ADVV(i,1:Nyb,m) + RXNV(i,1:Nyb,m)- dt*ADVzV(i,1:Nyb,m)
      call DGTSV(Nyb,1,DL2,D22,Du2,RHS2y,Nyb,matinfo)
      C(i,1:Nyb) = RHS2y
   end do

!------------------------------------


   P(1:Nxb, 1:Nyb,m) = C
end do
else
   do m = 1, Nz
      do j = 1, Nyb
         call N_mat(Diff_LR(1:Nxb,j), kx, Nxb, DL, D2, Du, BCX(j),0d0,0)
         RHS2 = (kx*Diff_LR(1:Nxb,1)*G(0:Nxb-1,j,m) + (1d0-2d0*kx*Diff_LR(1:Nxb,1))*G(1:Nxb,j,m) + kx*Diff_LR(1:Nxb,j)*G(2:Nxb+1,j,m))
         !This is going to need to be rethought through with
          RHS2 = RHS2 + ADVV(1:Nxb,j,m) + RXNV(1:Nxb,j,m)- dt*ADVzV(1:Nxb,j,m)
         call DGTSV(Nxb, 1, DL, D2, Du, RHS2, Nxb, halfmatinfo)
         C(1:Nxb, j ) = RHS2
      end do
      P(1:Nxb,1:Nyb,m) = C(1:Nxb,1:Nyb)
   end do
endif


call ghost_cells(P, A,B,D, G, Nxb, Nyb, Nz,Cup,jk,uvol,vvol)
end subroutine diffusion_subr_V
!----------------------------------------------------------------------
subroutine diffusion_subr_tb(Nxb,Nyb, Cold, C, Cs, RXN, kx, ky, Diff_TB,Diff_LR, ADV, A,B,D,Cup,jk,uvol, vvol)

integer:: Nxb, Nyb, i ,j,jk
double precision:: kx, ky,  Cup
double precision, dimension(0:Nxb+1,0:Nyb+1):: Cold, C, Cs,G, Diff_TB, Diff_LR
double precision, dimension(1:Nxb,1:Nyb):: RHS2,RXN, ADV
double precision, dimension(0:Nxb+1):: A,B,D
double precision, dimension(0:Nyb+1):: uvol
double precision, dimension(0:Nxb+1):: vvol
C=Cold

do j = 1,Nyb
   do i = 1, Nxb
      RHS2(i,j) = 2d0*ky*(Diff_TB(i,j+1)*(Cs(i,j+1)-Cs(i,j)) - Diff_LR(i,j)*(Cs(i,j)-Cs(i,j-1))) &
                + 2d0*kx*(Diff_TB(i+1,j)*(Cs(i+1,j)-Cs(i,j)) - Diff_LR(i,j)*(Cs(i,j)-Cs(i-1,j)))
   end do
end do

C(1:Nxb,1:Nyb) = Cold(1:Nxb,1:Nyb) + RHS2 + RXN + ADV

call ghost_cells(C, A, B, D, G, Nxb, Nyb, 1, Cup,jk,uvol,vvol)
C=G

end subroutine diffusion_subr_tb
!--------------------------------------------------------------------------
subroutine alpha_sq(t)
double precision:: t
double precision, dimension(1:Nxb,1:Nyb):: VolFac, Permeability
double precision:: Ftheta, Funbg
integer:: test_1
asq = 0.0

if (nobump == 0) then
asq(1:Nxb,1:Nyb,1) = (5000.)**2d0* BumpCon*(1.0-exp(-200.0*t))
asq(1:Nxb,1:Nyb,2) = (5000.)**2d0* BumpCon*(1.0-exp(-200.0*t))

!   do i = 1,nxb                                                                                      
!      do j = 1,nyb
         !asq(i,j,1) = (5000.)**2d0*(5d-1*tanh(10.*(x(i)-(L2-L2/10.))) - 5d-1*tanh(10.*(x(i)-(L2+L2/10.))))&
         !          *(-5d-1*tanh(50.*(y(j) - 0.05))+0.5)*(1.0-exp(-200.0*t))
         !(100.0*tanh((i-1.0)*xmax/nxb-rzL) - 100.0*tanh((i-1.0)*xmax/nxb-rzR)) &   
         !*(-0.5*tanh((j-0.5)*ymax/nyb-ymax/100.0-2.0) + 0.5)*(1.0-exp(-5.0*(ttotal)))   
        ! asq(i,j,2) =(5000.)**2*(5d-1*tanh(10.*(x(i)-(L2-L2/10.)))- 5d-1*tanh(10.*(x(i)-(L2+L2/10.))))&
          !       *(-5d-1*tanh(50.*(y(j) - 0.05))+0.5)*(1.0-exp(-200.0*t))
         !(100.0*tanh((i-0.5)*xmax/nxb-rzL) - 100.0*tanh((i-0.5)*xmax/nxb-rzR)) &        
         !*(-0.5*tanh((j-1)*ymax/nyb-ymax/100.0-2.0) + 0.5)*(1.0-exp(-5.0*(ttotal)))
      !end do
  ! end do
end if

test_1 = 1

do i = 1,nxb
   do j = 1,nyb

      
      if (thetag(i,j) > 5d-1) then

         VolFac(i,j) = thetag(i,j)*FC2
         
         if (VolFac(i,j) .ge. 1.0 ) then
            VolFac(i,j) = 0.999
         end if
         Ftheta = -(3.0*thetag(i,j)*FC2+sqrt(3.0)*sqrt(FC2*3.0**(3.0/4.0)*thetag(i,j)))/(-3.0**(3.0/4.0)+3.0*thetag(i,j)*FC2);
         Funbg = 2.0/3.0*(9.0**(1.0/3.0))*((FC1*3.0**(1.0/4.0)*Branchg(i,j)**2.0)**(1.0/3.0))/Branchg(i,j);
         Pore(i,j) = Funbg/(1+Ftheta);
         Diam(i,j)=  Ftheta*Pore(i,j);
         Permeability(i,j) = Diam(i,j)**2.0/Ls**2.0*1d0/(4d0*VolFac(i,j)) &
              *(-log(VolFac(i,j)) - 1.5 + 2.0*VolFac(i,j)-VolFac(i,j)**2.0/2.0)
         
         asq(i,j,1) = 1.0/Permeability(i,j) + asq(i,j,1)
         asq(i,j,2) = 1.0/Permeability(i,j) + asq(i,j,2)
         if (test_1 .gt. 0) then
            write(*,*) 'Thetag=', thetag(i,j),'Perm', Permeability(i,j),'Diam', Diam(i,j), 'asq(imj) =', asq(i,j,1)   
         test_1 = -2
         end if
         
      end if
   end do
end do

end subroutine alpha_sq
!-------------------------------------------------------------------------  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! FUNCTION: TRIMBT -- return the length of a string ignoring the first    
!                     blank or tab note: the ascii code for <tab> is 009
!
function trimbt( instring )
  integer      :: trimbt
  character(*) :: instring
  integer      :: slen
  trimbt = 1
  do
     if( (instring(trimbt:trimbt).eq.' '     ).or.  &
         (instring(trimbt:trimbt).eq.achar(9))      &
       )then
        trimbt = trimbt - 1
        exit
     end if

     if( slen.eq.len(instring) ) then
        exit
     end if

     trimbt = trimbt + 1
  end do

  return

end function trimbt
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! SUBROUTINE: INIT --- intialize things              
!
subroutine init
  character(40) :: filename    ! filename for writing
  integer       :: strlen      ! length of runname string
  integer       :: k
  integer       :: ycount
  real          :: hym,ytmp


  ! load the parameters from a file               
  !
write(*,*) 'loadparams'
  call loadparams

  ! put the rest of the data in the Data directory
  !                       
  runname = './Data/'//runname
  strlen = trimbt(runname)
call setparams
call fluidinit
end subroutine init
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  SUBROUTINE: FLUIDINIT --- initialize the fluid solver
!             
subroutine fluidinit
  real :: x,y
  integer :: i,j
call alpha_sq(0.0)
write(*,*) 'fluidin max=', maxval(asq)
call diff_init( asq )
call gc_init(gc_phidt,nxb,nyb,hb,0.0,-1.0,NMN,DIR,NMN,NMN)

  ! set the tolerance of the projection and momentum solver
  !

  call MG_settol( mgtol )

  call diff_settol( difftol )


  ! initialize with a Poiseuille flow with max vel as umax 
  ! zero outside the domain         

 !ADDED THIS FOR TESTING CONVERGENCE
write(*,*) 'umax=', umax
u(:,:,1)=0.0

  do j=1,nyb
     y = ymin + hb*(j-0.5)
     u(:,j,1) = umax/(ymax**2)*4.0*(y-ymin)*(ymax - y) !this is dimensional!
     
     ! note using the 0.01 instead of 4 since the ymax = 20 and not 1
  end do
write(*,*) 'wall shear rate=', umax/ymax**2.* 4d0*(ymax+ymin) !!CZ-this looks dimensional!!!
write(*,*) 'CZ wall shear rate=', 4.0*umax/(ymax*1d-2) 

write(*,*) 'umax=', umax
write(*,*) 'ymax=', ymax
write(*,*) 'ymin=', ymin
write(*,*) 'u - max=', maxval(u)
write(*,*) 'umax/uchar=', umax/uchar



!write(*,*) 'ymax=', u(0,:,1)**2 -u(2,:,1)**2
!write(*,*) 'ymin=', ymin
  ! set the vertical velocity and pressure to zero
  !
  u(:,:,2) = 0.0
  p(:,:) = 0.0


  ! nondimensionalize the velocities
  write(*,*) 'uchar=', uchar
  u(:,:,1) = u(:,:,1)/uchar
write(*,*) 'utildemax=', maxval(u(:,:,1))
!stop

end subroutine fluidinit
!-------------------------------------------------------------------------
subroutine poresize(Pore,Diam)
double precision, dimension(Nxb,Nyb):: Pore, Diam   !Pore size (inner length of hexagon), Fiber Diameter
integer:: i,j
double precision:: Ftheta, Funbg !function of theta and bg check maple code poresize for more info


!FC1 = 6.02d15
!FC2 = 4.275d-19*6.02d23*1d-8
Pore = 0.0 
Diam = 0.0 
do i = 1,Nxb
   do j = 1,Nyb
      if (thetag(i,j).gt.5d0) then
         Ftheta = -(3.0*thetag(i,j)*FC2+sqrt(3.0)*sqrt(FC2*3.0**(3.0/4.0)*thetag(i,j)))/(-3.0**(3.0/4.0)+3.0*thetag(i,j)*FC2);
            Funbg = 2.0/3.0*(9.0**(1.0/3.0))*((FC1*3.0**(1.0/4.0)*Branchg(i,j)**2.0)**(1.0/3.0))/Branchg(i,j);
            Pore(i,j) = Funbg/(1+Ftheta);
            Diam(i,j)=  Ftheta*Pore(i,j);
         end if
      end do
end do
end subroutine poresize
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!------Subroutine for gel componenets - ie interpolation must be used to filled in ghost cells because gel doesn't move-------
!-----This is an interpolation using 3 inner points (there for its a quadratic interpolation)---------------------------------
subroutine gel_ghost_cells(Gel_initial,Gel_ghost, Nxb, Nyb) 
integer:: Nxb, Nyb
double precision, dimension(0:Nxb+1,0:Nyb+1):: Gel_initial, Gel_ghost

Gel_ghost = Gel_initial
Gel_ghost(0,0:Nyb+1)     = 3d0*Gel_initial(1,0:Nyb+1)-3d0*Gel_initial(2,0:Nyb+1) + Gel_initial(3,0:Nyb+1)
Gel_ghost(Nxb+1,0:Nyb+1) = 3d0*Gel_initial(Nxb,0:Nyb+1)-3d0*Gel_initial(Nxb-1,0:Nyb+1) + Gel_initial(Nxb-2,0:Nyb+1)
Gel_ghost(0:Nxb+1,0)     = 3d0*Gel_initial(0:Nxb+1,1)-3d0*Gel_initial(0:Nxb+1,2) + Gel_initial(0:Nyb+1,3)
Gel_ghost(0:Nxb+1,Nyb+1) = 3d0*Gel_initial(0:Nxb+1,Nyb)-3d0*Gel_initial(0:Nxb+1,Nyb-1) + Gel_initial(0:Nxb+1,Nyb-2)


end subroutine gel_ghost_cells
!-------Hindered Diffusion Coefficients---------------------------
subroutine hinddifcoef
double precision, dimension(-1:Nxb+2,-1:Nyb+2):: phig_lr, phig_tb
integer:: phi_id
double precision, dimension(-1:Nxb+2,-1:Nyb+2):: rctb, rclr, rftb, rflr, rztb, rzlr
double precision, dimension(-1:nxb+2,-1:Nyb+2):: grad_phi_u_lr, grad_phi_u_tb,grad_phi_v_lr, grad_phi_v_tb
phig_lr = 0d0
phig_tb = 0d0

do phi_id = -1, Nxb+1
   phig_lr(phi_id+1,-1:Nyb+2) = 0.5*(phi_gpv(phi_id+1,-1:Nyb+2) + phi_gpv(phi_id,-1:Nyb+2))
end do
do phi_id = -1,Nyb+1
  phig_tb(-1:Nxb+2,phi_id+1) = 0.5*(phi_gpv(-1:Nxb+2,phi_id+1) + phi_gpv(-1:Nxb+2,phi_id))
end do

phig_lr(1,-1:Nyb+2)     = (15.0*phi_gpv(2,-1:Nyb+2) - 10.0*phi_gpv(3,-1:Nyb+2) + 3.0*phi_gpv(4,-1:Nyb+2))/8.0
phig_lr(Nxb+1,-1:Nyb+2) = (15.0*phi_gpv(Nxb,-1:Nyb+2) -10.0*phi_gpv(Nxb-1,-1:Nyb+2) + 3.0*phi_gpv(Nxb-2,-1:Nyb+2))/8.0
phig_tb(-1:Nxb+2,1)     = (15.0*phi_gpv(-1:Nxb+2,2) -10.0*phi_gpv(-1:Nxb+2,3) + 3.0*phi_gpv(-1:Nxb+2,4))/8.0
phig_tb(-1:Nxb+2,Nyb+1) = (15.0*phi_gpv(-1:Nxb+2,Nyb) -10.0*phi_gpv(-1:Nxb+2,Nyb-1) + 3.0*phi_gpv(-1:Nxb+2,Nyb-2))/8.0
!There needs to be a better correction for the (1:Nxb,1) location.  


grad_phi_u_lr = 0d0
grad_phi_v_tb = 0d0
grad_phi_u_lr(0:Nxb+2,-1:Nyb+2) = (1d0-phi_gpv(0:Nxb+2,-1:Nyb+2)) - (1d0-phi_gpv(-1:Nxb+1,-1:Nyb+2))
grad_phi_v_tb(-1:Nxb+2,0:Nyb+2) = (1d0-phi_gpv(-1:Nxb+2,0:Nyb+2)) -(1d0-phi_gpv(-1:Nxb+2,-1:Nyb+1))

!This is assuming a cubic form of interpolation (which goes down an order for the derivative)
grad_phi_v_tb(-1:Nxb+2,1) = 1.0/(24.0*hy)*(71.0*grad_phi_v_tb(-1:Nxb+2,2) +141.0*grad_phi_v_tb(-1:Nxb+2,3)&
                          + 93.0*grad_phi_v_tb(-1:Nxb+2,4))
grad_phi_v_tb(-1:Nxb+2,Nyb+1) = -1.0/(24.0*hy)*(71.0*grad_phi_v_tb(-1:Nxb+2,Nyb) +141.0*grad_phi_v_tb(-1:Nxb+2,Nyb-1)&
                          + 93.0*grad_phi_v_tb(-1:Nxb+2,Nyb-2))
!---------------------------------------------------------------------------------------------------------------------

grad_phi_v_lr = 0d0
grad_phi_u_tb = 0d0
grad_phi_v_lr(0:Nxb+1,0:Nyb+1) = 0.25*(grad_phi_v_tb(0:Nxb+1,0:Nyb+1) + grad_phi_v_tb(1:Nxb+2,0:Nyb+1)&
                               + grad_phi_v_tb(0:Nxb+1,-1:Nyb) + grad_phi_v_tb(1:Nxb+2,-1:Nyb))
grad_phi_u_tb(0:Nxb+1,0:Nyb+1) = 0.25*(grad_phi_u_lr(0:Nxb+1,0:Nyb+1) + grad_phi_u_lr(1:Nxb+2,0:Nyb+1)&
                               + grad_phi_u_lr(0:Nxb+1,-1:Nyb) + grad_phi_u_lr(1:Nxb+2,-1:Nyb))


rc = 10.1*phi_gpv(0:Nxb+1,0:Nyb+1) + 1d0
rf = 10.1*phi_gpv(0:Nxb+1,0:Nyb+1) + 1d0
rz = 30.56*phi_gpv(0:Nxb+1,0:Nyb+1) + 1d0

rctb = 10.1*phig_tb + 1d0
rclr = 10.1*phig_lr + 1d0

rftb = 10.1*phig_tb + 1d0
rflr = 10.1*phig_lr + 1d0

rztb = 30.56*phig_tb + 1d0
rzlr = 30.56*phig_lr + 1d0

Hdc_lr = 1d0/(1d0+(rclr(0:nxb+1,0:Nyb+1) - 1d0)*phig_lr(0:nxb+1,0:Nyb+1))
Hdz_lr = 1d0/(1d0+(rzlr(0:nxb+1,0:Nyb+1)  - 1d0)*phig_lr(0:nxb+1,0:Nyb+1))
Hdf_lr = 1d0/(1d0+(rflr(0:nxb+1,0:Nyb+1)  - 1d0)*phig_lr(0:nxb+1,0:Nyb+1))

Hdc_tb = 1d0/(1d0+(rctb(0:nxb+1,0:Nyb+1)  - 1d0)*phig_tb(0:nxb+1,0:Nyb+1))
Hdz_tb = 1d0/(1d0+(rztb(0:nxb+1,0:Nyb+1)  - 1d0)*phig_tb(0:nxb+1,0:Nyb+1))
Hdf_tb = 1d0/(1d0+(rftb(0:nxb+1,0:Nyb+1)  - 1d0)*phig_tb(0:nxb+1,0:Nyb+1))

Huc_tb = (1d0-phig_tb)/(1d0 + (rctb-1d0)*phig_tb) 
Huc_lr = (1d0-phig_lr)/(1d0 + (rclr-1d0)*phig_lr)
Huz_tb = (1d0-phig_tb)/(1d0 + (rztb-1d0)*phig_tb)
Huz_lr = (1d0-phig_lr)/(1d0 + (rzlr-1d0)*phig_lr)
Huf_tb = (1d0-phig_tb)/(1d0 + (rftb-1d0)*phig_tb)
Huf_lr = (1d0-phig_lr)/(1d0 + (rflr-1d0)*phig_lr)

uvolc_tb = grad_phi_u_tb/(1d0-phig_tb)
uvolc_lr = grad_phi_u_lr/(1d0-phig_lr)
uvolf_tb = grad_phi_u_tb/(1d0-phig_tb)
uvolf_lr = grad_phi_u_lr/(1d0-phig_lr)
uvolz_tb = grad_phi_u_tb/(1d0-phig_tb)
uvolz_lr = grad_phi_u_lr/(1d0-phig_lr)
vvolc_tb = grad_phi_v_tb/(1d0-phig_tb)
vvolc_lr = grad_phi_v_lr/(1d0-phig_lr)
vvolf_tb = grad_phi_v_tb/(1d0-phig_tb)
vvolf_lr = grad_phi_v_lr/(1d0-phig_lr)
vvolz_tb = grad_phi_v_tb/(1d0-phig_tb)
vvolz_lr = grad_phi_v_lr/(1d0-phig_lr)

end subroutine hinddifcoef
!--------------------------------------------------------------------------
subroutine alpha_function(btypeA)
double precision, dimension(0:Nxb+1):: alphafun
double precision:: Diffa
integer:: idxc1, idxc2,indicator, btypeA
double precision, dimension(0:Nxb+1):: A_line, B_line, C_line, Line_tot
double precision:: A_space, B_space, L2_1, L2_2, L3_1, L3_2
A_space = 21.33/100
B_space = 7.5/100;

L2_1 = L2+A_space
L2_2 = L2_1+B_space
L3_1 = L2_2+A_space
L3_2 = L3_1+B_space

if (btype ==0) then
alphafun = 10*((5e-1*tanh(5*(x-L2))-5e-1*tanh(5*(x-L3))));
indicator = 0  
i = 1
do while (indicator < 1)
     if (x(i)>L2) then
         idxc1 = i
         indicator = 2
     else
        i = i+1

     end if
   
  end do

indicator = 0
i = Nxb/2
do while (indicator < 1)
      if (x(i) > L3) then
         idxc2 = i
         indicator = 2
      else
         i = i + 1
      end if

end do



    Diffa = alphao - alphafun(idxc1);
    Alpha(0:idxc1) = alphafun(0:idxc1);
    Alpha(idxc2:Nxb+1) = alphafun(idxc2:Nxb+1);
    Alpha(idxc1+8:idxc2-8) = alphao;
    Alpha(idxc1+1) = Alpha(idxc1) + Diffa/32;
    Alpha(idxc1+2) = Alpha(idxc1+1) + Diffa/16;
    Alpha(idxc1+3) = Alpha(idxc1+2) + Diffa/8;
    Alpha(idxc1+4) = Alpha(idxc1+3) + Diffa/4;
    Alpha(idxc1+5) = Alpha(idxc1+4) + Diffa/4;
    Alpha(idxc1+6) = Alpha(idxc1+5) + Diffa/4;
    Alpha(idxc1+7) = Alpha(idxc1+6) + Diffa/32;
    
    do j = 1,8
        Alpha(idxc2-8+j) = Alpha(idxc1+8-j);
    end do
else
   do i = 0, Nxb+1

         if (btypeA > 1) then
            A_line(i) = ((.5*tanh(25*(x(i)-L2))-5e-1*tanh(25*(x(i)-L2_1))))!*(-.5*tanh(25.*(y(j) - .05))+0.5);
            B_line(i) = ((.5*tanh(25*(x(i)-L2_2))-5e-1*tanh(25*(x(i)-L3_1))))!*(-.5*tanh(25.*(y(j) - .05))+0.5);
            C_line(i) = ((.5*tanh(25*(x(i)-L3_2))-5e-1*tanh(25d0*(x(i)-L3))))!*(-.5*tanh(25.*(y(j) - .05))+0.5);
         elseif (btypeA ==1) then
            line_tot(i) = ((5d-1*tanh(5d0*(x(i)-L2))-5d-1*tanh(5d0*(x(i)-L3))))!*(-5d-1*tanh(5d1*(y(j) - 0.05))+5d-1);
         end if
      !end do
   end do
if (btypeA==1) then
   Alpha = alphao*line_tot
elseif (btypeA ==2) then
   Alpha = alphao*(A_line+B_line+C_line)
elseif(btypeA==3) then
   Alpha = alphao*(A_line+C_line)
end if
end if




end subroutine 



end module gel_mod
!--------------------------------------------------------------------------    
!--------------------------------------------------------------------------    



! FINALLY, THE MAIN PROGRAM                                                    
!                                                                              
program clot
  use gel_mod
  use param_mod
!  use clot_mod
write(*,*) 'start driver'
call init
call driver
write(*,*) 'end driver'
end program clot
