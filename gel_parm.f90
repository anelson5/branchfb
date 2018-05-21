module gel_parm

use grid_mod
use param_mod
implicit none
integer,parameter:: Nz = 100
double precision, dimension(1:Nz):: hz
double precision:: Lx, Ly,Lz, hx, hy             !domain parameters
double precision:: to, tend, t, dt1, dto, timerecord, mu_c,pi, nu
double precision:: kx, ky, Nxbd, Nybd, Nzd
double precision:: DD, dalpha, dbeta, dtc
double precision:: Kappaon, Kappaoff, kappaoc
integer::i, j, N1, N2, ix,jk,initialdata, iy, indexrecord, recordidx, dogelation
integer:: dosource,n, jz, flowon
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
double precision:: Zup,ke, ZupD
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
double precision, dimension(Nz):: z,phi
double precision, dimension(1:Nxb,1:Nyb,1:Nz):: W, ADVw, Funp, ADVzW, Rxn_zW,LW_z
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
!----Saving files----------------   
!character(100):: Nyb_char, Dtc_char
!character(6)::fileend
!character(50)::Zfile, Efile, Ffile, Rfile, Rgfile,Tfile, Sfile, Vfile
!character(50):: Wfile,Bfile, Thfile, Bgfile, Thgfile, Bsfile, Tsfile, Parsfile
!integer:: Wunit=47, Eunit=48, Vunit = 49, Runit=50, Rgunit=51, Tunit=52, Sunit=53
!integer:: Thunit=54, Bunit=55, Thgunit=56, Bgunit=57, Bsunit = 58, Tsunit=59
!integer:: Z2unit= 61, Funit = 60, Parsunit = 61

CONTAINS

!-----------SETPARAMS-GEL-----------------------------------
!This subroutine sets the values for the gelation parameters
!-----------------------------------------------------------

subroutine setparams_gel
implicit none
write(*,*) 'setparams_gel'
!Nxb,Nyb
Ly  = dheight
Lx  = 4d0*Ly
NxbD = dble(Nxb)
NybD = dble(Nyb)
NzD = dble(Nz)
!Grid parameters
!Lx = 6d0
!Ly = 2.5d-1
hx = Lx/NxbD
hy = Ly/NybD
!Injury Site
L2 = 2.5d0
L3 = 3.5d0
!x,y
xo = (/ (hx*(ix) - hx/2d0, ix = 0,Nxb+1)/)
y  = (/ (hy*(iy) - hy/2d0, iy = 0,Nyb+1)/)

call fun_hz(hz,Nz,z)

NxbD      = dble(Nxb)
NybD      = dble(Nyb)
NzD      = dble(Nz)
pi       = ACOS(-1d0)
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
mu_c       = 6d0 !Sharpness of Alpha transition between 0 and alpha_o              
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



JK = -1
end subroutine setparams_gel

!-------Fun_hz------------------------
!This subroutine makes the variable grid spacing for z.
!
subroutine fun_hz(hz,Nz,z)
integer:: i, Nz
double precision:: hzo, hmax, hend, hmin
double precision, dimension(Nz)::ztemp, hh ,hz, z

hzo      = 1d0/(dble(Nz)-1d0)
ztemp    = (/(hzo*dble(I),I = 0,Nz-1) /)
hmin     = 5d-3
hmax     = hzo
hend     = 1d-3
hh       = 5d-1*(1d0+tanh((7d-1-ztemp)/15d-2))
hz(1:Nz) = hh/sum(hh)

end subroutine fun_hz



end module gel_parm
