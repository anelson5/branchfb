!
! MODULE: DIFFUSION_MOD --- this module handles updating the heat equation 
!                           on an irregular domain.  This is written for 
!                           the case of XXXX
!                            
module diffusion_mod
  use grid_mod
  use param_mod
  implicit none

  PRIVATE

  ! indicies to loop over --- used for different boundary conditions
  !
  integer :: iminU, imaxU
  integer :: jminU, jmaxU
  integer :: iminV, imaxV
  integer :: jminV, jmaxV

  ! stencil for the irregular laplacian
  !
  real :: lapstenU(1:nxb+1, 1:nyb+1, 5)
  real :: lapstenV(1:nxb+1, 1:nyb+1, 5)
  
  ! alpha^2 values at the edges
  !
  real :: asq(1:nxb+1, 1:nyb+1, 2)

  ! max number of iterations
  !
  integer :: maxiter = 1000
  real    :: tol     = 1e-10

  ! list of public routines and variables
  !
  PUBLIC :: diff_init, step_diffusion, diff_settol, new_alpha
 
  
CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! SUBROUTINE: DIFF_INIT --- compute the stencils and the masks
!
subroutine diff_init( asq_in )
  real :: asq_in(1:nxb+1, 1:nyb+1, 2)
    write(*,*) 'asq_in (diff_init) = ', maxval(asq_in)
  ! set regular stencils everywhere
  !
  lapstenU(:,:, 1 ) = -4.0
  lapstenU(:,:,2:5) =  1.0

  lapstenV(:,:, 1 ) = -4.0
  lapstenV(:,:,2:5) =  1.0

  ! record alpha values
  !
  asq = asq_in


  ! modify the stencil near the boundary 
  ! & record the number of unknowns
  !
  call stencilbc

end subroutine diff_init
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! SUBROUTINE: NEW_ALPHA get the new alpha values
!
subroutine new_alpha( asq_new )
  real :: asq_new(1:nxb+1, 1:nyb+1, 2)

  ! record alpha values
  !
  asq = asq_new

end subroutine new_alpha
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! SUBROUTINE: DIFF_SETTOL --- set the tolerance of the diffusion solver
!
subroutine diff_settol( newtol )
  real :: newtol

  tol = newtol
  
end subroutine diff_settol
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! SUBROUTINE: STENCILBC --- modify the stencil near the boundaries
!                           and decide how many unknowns
!                           
!                           This is written for LR boundaries periodic and
!                           TB boundaries homogeneous Dirichlet
!
subroutine stencilbc

  ! limits for the unknowns
  !
  iminU = 2;  imaxU = nxb+1  ! now have dirichlet L, neumann R
  jminU = 1;  jmaxU = nyb


  ! bottom boundary
  !
  !lapstenU(iminU:imaxU, jminU, 1) = lapstenU(iminU:imaxU, jminU, 1) - 1.0
  !lapstenU(iminU:imaxU, jminU, 4) =  0.0

  ! top boundary
  !
  !lapstenU(iminU:imaxU, jmaxU, 1) = lapstenU(iminU:imaxU, jmaxU, 1) - 1.0
  !lapstenU(iminU:imaxU, jmaxU, 2) =  0.0

  !!!!!
  !!!!! now the same for v
  !!!!!

  ! limits for the unknowns
  !
  iminV = 1;  imaxV = nxb
  jminV = 2;  jmaxV = nyb 

  !!!
  !
  ! no need to modify top/bottom boundary as long as values are set
  ! to dirichlet condtions for j=1 and j=nyb+1
  !
  !!!

  !!!
  !
  ! other boundaries are periodic, so they will be handled with ghost cells
  !
  !!!

end subroutine stencilbc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! SUBROUTINE: STEP_DIFFUSION --- take a time step of size dt of the 
!                                diffusion equation
!                                   u   = coeff*Lu + f
!                                    t
!                                using Crank-Nicholson time step
!
subroutine step_diffusion(u,f,coeff,dt,UVflag)
  real    :: u(-1:nxb+2,-1:nyb+2)        ! begining and end solution
  real    :: f(-1:nxb+2,-1:nyb+2)        ! forcing term
  real    :: coeff                       ! diffusion coefficient 
  real    :: dt                          ! time step
  integer :: UVflag                      ! set to 1 for u 
                                         ! set to 2 for v
 
  real    :: d                          ! coefficient

  real :: rhs(1:nxb+1,1:nyb+1)          ! rhs for solve in CN
  integer :: imin,imax,jmin,jmax        ! range of unknowns
  real ::  st(1:nxb+1,1:nyb+1,5)        ! stencil for (I-d*L)

  ! coefficient in front of laplacian
  !
  d = 0.5*dt*coeff/(hb*hb)
  
  ! set imin,imax,jmin,jmax
  !
  if( UVflag.eq.1 ) then
 
     ! record number of unknowns
     !
     imin = iminU; imax = imaxU
     jmin = jminU; jmax = jmaxU

     ! set up the right side
     !
     call setrhs(rhs,u,f,d,dt,lapstenU,imin,imax,jmin,jmax,UVflag)
     !write(*,*) 'setrhs = ', f
     !write(*,*) 'lapstenu=', lapstenU
     !write(*,*) 'u=', u(1,:)
     ! set up the stencil
     !
!     write(*,*) 'asq diffuionmod=', maxval(asq(:,:,1)), maxval(asq(:,:,2))
     st(:,:,1  ) = 1.0 - d/Re*lapstenU(:,:, 1  ) + dt/Re*asq(:,:,1)
     st(:,:,2:5) =     - d/Re*lapstenU(:,:, 2:5)
     !write(*,*) UVflag
     ! solve
     !
     call solve(u,rhs,st,imin,imax,jmin,jmax,UVflag)
     
  else

     ! record number of unknowns
     !
     imin = iminV; imax = imaxV
     jmin = jminV; jmax = jmaxV
 
     ! set up the right side
     !
     call setrhs(rhs,u,f,d,dt,lapstenV,imin,imax,jmin,jmax,UVflag)
     
     ! set up the stencil
     !
     st(:,:,1  ) = 1.0 - d/Re*lapstenV(:,:, 1  ) + dt/Re*asq(:,:,2)
     st(:,:,2:5) =     - d/Re*lapstenV(:,:, 2:5)
 !    write(*,*) maxval(lapstenV), minval(lapstenV)
     ! solve
     !
     call solve(u,rhs,st,imin,imax,jmin,jmax,UVflag)
  end if

end subroutine step_diffusion
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! SUBROUTINE: SETRHS --- compute the right side
!
subroutine setrhs(rhs,u,f,d,dt,lapsten,imin,imax,jmin,jmax,UVflag)

  real    :: rhs(1:nxb+1,1:nyb+1)          ! rhs for solve in CN
  real    :: u(-1:nxb+2,-1:nyb+2)          ! solution from input
  real    :: f(-1:nxb+2,-1:nyb+2)          ! forcing term
  real    :: d                             ! coefficient
  real    :: dt                            ! time step
  real    :: lapsten(1:nxb+1, 1:nyb+1, 5)  ! stencil

  integer :: imin,imax,jmin,jmax           ! range of unknowns
  integer :: UVflag			   ! 1 for u 2 for v
  integer :: i,j                           ! loop counters

  ! fill ghost cells
  !

  if (UVflag.eq.1) then
     call fillgcu(u)
  else
     call fillgcv(u)
  end if

  ! set the rhs
  !
!write(*,*) 'f=', maxval(abs(f))
  do i=imin,imax
     do j=jmin,jmax
        rhs(i,j) = dt*f(i,j) + (1.0 + d*lapsten(i,j,1))*u(i,j)         &
             + d*( lapsten(i,j,2)*u(i,j+1) + lapsten(i,j,3)*u(i+1,j) + &
                   lapsten(i,j,4)*u(i,j-1) + lapsten(i,j,5)*u(i-1,j))
     end do
  end do

end subroutine setrhs
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! SUBROUTINE: SOLVE --- solve equation
!                         (I - d*L)u = rhs
!                       by an iterative scheme
!
subroutine solve(u,rhs,st,imin,imax,jmin,jmax,UVflag)

  real    ::    u(-1:nxb+2,-1:nyb+2)        ! solution
  real    ::    uold(-1:nxb+2,-1:nyb+2)     ! solution
  real    ::  rhs( 1:nxb+1, 1:nyb+1)        ! right hand side

  integer :: imin,imax,jmin,jmax           ! range of unknowns
  
  integer :: UVflag			   ! 1 for u 2 for v  

  integer :: i,j

  real ::  st(1:nxb+1,1:nyb+1,5)   ! stencil for (I-d*L)
  real :: res(1:nxb+1,1:nyb+1)     ! residual
  real :: maxres                   ! max of the residual

  integer :: k                 ! counter for number of iterations
  integer :: itsor

  ! begin iterations
  !
  uold=u
  k = 0
  do
     ! relax
     !
     call gsrelax(u,rhs,st,imin,imax,jmin,jmax,UVflag)
    !write(*,*) 'max =', minval(abs(u))
     
     ! compute residual
     !
     !call residual(res,u,rhs,st,imin,imax,jmin,jmax,UVflag)

     res(imin:imax,jmin:jmax)=u(imin:imax,jmin:jmax)-uold(imin:imax,jmin:jmax)
     k = k+1

     ! check if we have converged
     !
     maxres = maxval(abs(res(imin:imax,jmin:jmax)))

     if( maxres.lt.tol ) then
        !write(*,*) 'converged in ',k
        ! write(50,*) k
        exit
     end if

     ! check if we have exceeded the maximum number of iterations
     !
     if( k.ge.maxiter ) then
        write(*,*) 'Regular Diffusion solve ---failed to converge '
        write(*,*) '  maximum iterations = ', maxiter
        write(*,*) '  tolerance          = ', tol
        write(*,*) '  max residual       = ', maxval( abs(res) )
        write(*,*) '  stopping'
        
!!$        do i=imin,imax
!!$           do j=jmin,jmax
!!$              write(30,*) res(i,j)
!!$           end do
!!$        end do

        stop

     end if
     uold=u
  end do

!write(*,*)'SOR iterations',k
end subroutine solve
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! SUBROUTINE: RESIDUAL --- compute the resudual by
!                            res = rhs - Lu 
!                          where st is the stencil of L
!
subroutine residual(res,u,rhs,st,imin,imax,jmin,jmax,UVflag)
  real :: res( 1:nxb+1, 1:nyb+1)    ! residual
  real ::   u(-1:nxb+2,-1:nyb+2)    ! function
  real :: rhs( 1:nxb+1, 1:nyb+1)    ! right hand side'
  real ::  st( 1:nxb+1, 1:nyb+1, 5) ! stencil
  
  integer :: imin,imax,jmin,jmax    ! range of unknowns

  integer :: UVflag		    ! 1 for u 2 for v

  integer :: i,j                   ! loop counters

  ! fill ghost cells
  !
  if (UVflag.eq.1) then
     call fillgcu(u)
  else 
     call fillgcv(u)
  end if

  ! compute the residual
  !
  do i=imin,imax
     do j=jmin,jmax
        res(i,j) = rhs(i,j) - st(i,j,1)*u(i,j) - st(i,j,2)*u(i,j+1) - &
             st(i,j,3)*u(i+1,j) - st(i,j,4)*u(i,j-1) - st(i,j,5)*u(i-1,j)
     end do
  end do

end subroutine residual
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! SUBROUTINE: GSRELAX --- peform a single relaxation of u
!
subroutine gsrelax(u,rhs,st,imin,imax,jmin,jmax,UVflag)
  real ::   u(-1:nxb+2,-1:nyb+2   )  ! current solution
  real :: rhs( 1:nxb+1, 1:nyb+1   )  ! right hand side
  real ::  st( 1:nxb+1, 1:nyb+1, 5)  ! stencil of operator

  integer :: imin,imax,jmin,jmax     ! range of unknowns

  integer :: UVflag		     ! 1 for u 2 for v

  real :: uold(-1:nxb+2,-1:nyb+2)    ! store old solution

  integer :: i,j                     ! loop counters

  ! fill ghost cells
  !
  if (UVflag.eq.1) then
     call fillgcu(u)
  else
     call fillgcv(u)
  end if

  ! record old value
  !
  uold = u

  ! update by sor
  !
!write(*,*) maxval(st(:,:,1)), minval(st(:,:,1))

  do i=imin,imax
     do j=jmin,jmax
        u(i,j) = (1.0-omegaf)*uold(i,j) + omegaf*                       &
                 (rhs(i,j) - st(i,j,2)*u(i,j+1) - st(i,j,3)*u(i+1,j)  &
                           - st(i,j,4)*u(i,j-1) - st(i,j,5)*u(i-1,j)  &
                 )/st(i,j,1)
     end do
  end do

!write(*,*) 'rhs=', maxval(rhs), minval(rhs)
end subroutine gsrelax
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! SUBROUTINE: FILLGCU --- fill ghost cells for u 
!			  left boundary dirichlet inflow condition
!			  right boundary neumann outflow condtion 
!                         top/bottom are homogeneous dirichlet
!                        
!
!
subroutine fillgcu( v )
  real :: v(-1:nxb+2,-1:nyb+2)

  v(0, 1:nyb) = 2.0*hb*v(1,1:nyb) -v(2, 1:nyb)
  v(nxb+2 , 1:nyb) = v(nxb , 1:nyb)

  v(0:nxb+1,nyb+1) = -v(0:nxb+1,nyb)
  v(0:nxb+1,0) = -v(0:nxb+1,1)



end subroutine fillgcu
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! SUBROUTINE: FILLGCV --- fill ghost cells for v 
! 			  left boundary homogeneous dirichlet
!			  right boundary neumann outflow 
!                         top/bottom are homogeneous dirichlet
!                        
!
subroutine fillgcv( v )
  real :: v(-1:nxb+2,-1:nyb+2)

  v(  0  , 1:nyb) = -v(1, 1:nyb)
  v(nxb+1, 1:nyb) = v(nxb , 1:nyb)

  v(1:nxb,   0 ) = -v(1:nxb,  2)
  v(1:nxb, nyb+2) = -v(1:nxb, nyb)
end subroutine fillgcv
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


end module diffusion_mod


