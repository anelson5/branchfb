module mg_mod
  use bobconst_mod
  use grid_mod
  use gc_mod
  PRIVATE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! TYPE - MGMESH: define a type containging pointers to arrays
!  
  type mgmesh
     real, pointer :: u(:,:)
     real, pointer :: f(:,:)
     real, pointer :: r(:,:)
  end type mgmesh
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! TYPE - STENCIL:  define a 3x3 array as a type to make an array of stencils
!
  type stencil
     real :: S(-1:1,-1:1)
  end type stencil
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  
  ! multigrid module local variables
  !
  logical :: isnotInit = .true. ! state of initialization of module
  integer :: nu1=1              ! number of relaxations up and down Vcycles
  integer :: nu2=1       
  real    :: tol=1e-14          ! convergence tolerance
  integer :: maxiter=120        ! maximum iterations allowed

  type (mgmesh) :: grid(0:qx) ! data structure for coarser grids
  type (gcm)    :: gc(0:qx)   ! boundary data for all grids
  type (stencil) :: A(0:qx)   ! stencils for all grid levels
    
  ! list of public routines
  !
  public :: MG_solve, maxnorm, MG_settol

CONTAINS


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! SUBROUTINE - MG_SOLVE: solve the system (c1 I - c2 L)u=f, where
!                        I     - identity operator
!                        L     - 5 point, standard laplacian
!                        c1,c2 - scalars  
!
  subroutine MG_solve(u,f,c1,c2,gcu)
    real,target :: u(0:,0:)
    real :: f(0:,0:)
    real :: c1,c2
    type (gcm) :: gcu

    integer :: i
    integer :: n

    real :: r0,rold,rnew

    ! check for initialization
    ! 
    if( isnotInit ) call MG_init


    ! set up the stencils for the operators
    !
    do i=0,qx
       A(i)%S( 0, 1) = -c2
       A(i)%S(-1, 0) = -c2
       A(i)%S( 0, 0) = 4.0**(qx-i)*hsq*c1 + 4.0*c2
       A(i)%S( 1, 0) = -c2
       A(i)%S( 0,-1) = -c2
    end do

    
    ! initialize stuff
    !
    grid(qx)%u => u  
    gc(qx) = gcu

    ! multiply f through by h^2
    !
    grid(qx)%f(1:nxb,1:nyb) = hsq*f(1:nxb,1:nyb)


    ! initlize the ghostcells for the coarser grids
    !
    do i=0,qx-1
       call gc_init(gc(i), nxc*2**i, nyc*2**i, hb*2**(qx-i),   &
                    c1,c2,gc(qx)%type(LFT),                    &
                    gc(qx)%type(RGT), gc(qx)%type(BOT),         &
                    gc(qx)%type(TOP)  )
    end do
    
    n = 1
    call resid( grid(qx)%u, grid(qx)%f, grid(qx)%r, A(qx), gc(qx) )
    r0 = maxnorm( grid(qx)%r )
    do
       if( (r0 < tol).OR.(n>maxiter) ) exit

       call vcycle
       call resid( grid(qx)%u, grid(qx)%f, grid(qx)%r, A(qx), gc(qx) )
       r0 = maxnorm( grid(qx)%r )
       n = n+1
    end do
    
    !write(*,*) 'Multigrid V-cycle count ', n-1

    if( n > maxiter ) then
       write(*,*) 'Multigrid failed to converge!!  Terminating'
       write(*,*) '  Residual is  ', r0
       write(*,*) '  Tolerance is ', tol
       write(*,*) '  max Vcycles  ', maxiter
       stop
    end if

    
  end subroutine MG_solve
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! SUBROUTINE - MG_INIT: initialize the solver by allocating memory
!
  subroutine MG_init
    integer :: i

    ! we are now initialized
    !
    isnotInit = .false.    

    ! allocate memory for the multigrid solver to use
    !
    allocate(grid(qx)%f(0:nxb+1,0:nyb+1)) 
    allocate(grid(qx)%r(0:nxb+1,0:nyb+1))

    do i=0,(qx-1)
       allocate( grid(i)%u( 0:nxc*2**i+1, 0:nyc*2**i+1 ) )
       allocate( grid(i)%f( 0:nxc*2**i+1, 0:nyc*2**i+1 ) )
       allocate( grid(i)%r( 0:nxc*2**i+1, 0:nyc*2**i+1 ) )
       
       grid(i)%u = 0.0
       grid(i)%f = 0.0
       grid(i)%r = 0.0
    end do    

    ! initialize the stencils to zero
    !
    do i=0,qx
       A(i)%S = 0.0
    end do

  end subroutine MG_init
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! SUBROUTINE: MG_SETTOL --- set the tolerance of the multigrid solver
!
subroutine MG_settol( newtol )
  real :: newtol

  tol = newtol

end subroutine MG_settol
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! SUBROUTINE - RESTRICT: restrict from coarse to fine grid
!                        and multiply by 4
!
  subroutine restrict(uc,uf)
    real :: uc(0:,0:), uf(0:,0:)
    
    integer :: ic,jc
    integer :: if,jf

    ! loop over the coarse grid
    !
    do ic=1,ubound(uc,1)-1
       if = 2*ic
       do jc=1,ubound(uc,2)-1
          jf = 2*jc
          uc(ic,jc) = uf(if,jf) + uf(if-1,jf) + uf(if,jf-1) + uf(if-1,jf-1)
       end do
    end do
  end subroutine restrict
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! SUBROUTINE - CORRECT: takes the error on the coarse grid and adds
!                       it to the approximate solution on the fine 
!                       grid
!
  subroutine correct(ec,uf,gce)
    real :: ec(0:,0:), uf(0:,0:)
    type (gcm) :: gce

    integer :: if, jf
    integer :: ic, jc
    integer :: nx,ny
    real :: a(4) = (/ 9.0/16.0, 3.0/16.0, 3.0/16.0, 1.0/16.0 /)
    
    nx = ubound(uf,1)-1
    ny = ubound(uf,2)-1

    call gc_interp(ec,gce)

    do if=1,nx,2
       ic = (if+1)/2
       do jf=1,ny,2
          jc = (jf+1)/2
          uf(if,jf) = uf(if,jf) + a(1)*ec(ic,jc) + a(2)*ec(ic-1,jc) + &
                      a(3)*ec(ic,jc-1) + a(4)*ec(ic-1,jc-1)
       end do
       
       do jf=2,ny,2
          jc = jf/2
          uf(if,jf) = uf(if,jf) + a(1)*ec(ic,jc) + a(2)*ec(ic-1,jc) + &
                      a(3)*ec(ic,jc+1) + a(4)*ec(ic-1,jc+1)
       end do
    end do
    
    do if=2,nx,2
       ic = if/2
       do jf=1,ny,2
          jc = (jf+1)/2
          uf(if,jf) = uf(if,jf) + a(1)*ec(ic,jc) + a(2)*ec(ic+1,jc) + &
                      a(3)*ec(ic,jc-1) + a(4)*ec(ic+1,jc-1)
       end do
       
       do jf=2,ny,2
          jc = jf/2
          uf(if,jf) = uf(if,jf) + a(1)*ec(ic,jc) + a(2)*ec(ic+1,jc) + &
                      a(3)*ec(ic,jc+1) + a(4)*ec(ic+1,jc+1)
       end do
    end do
    
  end subroutine correct
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! SUBROUTINE - RELAX: relax on either red points or black point 
!                     color = 0 or 1 --- 0=even 1=odd
!
  subroutine relax(u,f,sten,par,gcu)
    real :: u(0:,0:), f(0:,0:)
    type (stencil) :: sten
    integer :: par
    type (gcm) :: gcu

    integer :: par2
    integer :: i,j
    integer :: nx,ny

    nx = ubound(u,1)-1
    ny = ubound(u,2)-1

    call gc_relax(u,f,gcu)

    do i=1+par,nx,2
       do j=1,ny,2
          u(i,j) = (f(i,j) &
                    - sten%S(-1,0)*u(i-1,j) - sten%S(1, 0)*u(i+1,j) &
                    - sten%S( 0,1)*u(i,j+1) - sten%S(0,-1)*u(i,j-1) &
                   ) / sten%S(0,0)
       end do
    end do
    
    par2 = mod(par+1,2)

    do i=1+par2,nx,2
       do j=2,ny,2
          u(i,j) = (f(i,j) &
                    - sten%S(-1,0)*u(i-1,j) - sten%S(1,0)*u(i+1,j) &
                    - sten%S(0,1)*u(i,j+1) - sten%S(0,-1)*u(i,j-1) &
                   ) / sten%S(0,0)
       end do
    end do

  end subroutine relax
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! SUBROUTINE - REDBLACK: perform redblack gauss sidel nu times 
!                     
!
subroutine redblack(u,f,sten,gcu,nu)
  real :: u(0:,0:), f(0:,0:)
  type (stencil) :: sten
  type (gcm) :: gcu
  integer :: nu

  integer :: i

  do i=1,nu
     call relax(u,f,sten,0,gcu)   ! relax even cells
     call relax(u,f,sten,1,gcu)   ! relax odd cells
  end do
end subroutine redblack

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! SUBROUTINE - RESID: compute the residual
!
  subroutine resid(u,f,r,sten,gcu)
    real :: u(0:,0:), f(0:,0:), r(0:,0:)
    type (stencil) :: sten
    type (gcm) :: gcu

    integer :: i,j
    integer :: nx,ny

    nx = ubound(u,1)-1
    ny = ubound(u,2)-1

    call gc_op(u,gcu)

    do i=1,nx
       do j=1,ny
          r(i,j) = f(i,j) - &
                   sten%S(-1,0)*u(i-1,j) - sten%S(1,0)*u(i+1,j) - &
                   sten%S(0,-1)*u(i,j-1) - sten%S(0,1)*u(i,j+1) - &
                   sten%S(0,0)*u(i,j)
       end do
    end do

  end subroutine resid
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! FUNCTION - MAXNORM: compute the maxnorm
!
  real function maxnorm( u )
    real :: u(0:,0:)
    integer :: nx, ny
    nx = ubound(u,1)-1
    ny = ubound(u,2)-1

    maxnorm = maxval( abs( u(1:nx,1:ny) ) )

  end function maxnorm
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! SUBROUTINE - VCYCLE: perform one vcycle
!
  subroutine vcycle
    integer :: L
   
    ! down the V
    !
    do L=qx,2,-1
       call redblack(grid(L)%u, grid(L)%f, A(L), gc(L), nu1)
       call resid(grid(L)%u, grid(L)%f, grid(L)%r, A(L), gc(L) )
       call restrict(grid(L-1)%f,grid(L)%r)
       grid(L-1)%u = 0
    end do
    
    call redblack(grid(1)%u, grid(1)%f, A(1), gc(1), 10)
    call resid(grid(1)%u, grid(1)%f, grid(1)%r, A(1), gc(1) )

    ! up the V
    !
    do L=2,qx
       call correct( grid(L-1)%u, grid(L)%u, gc(L-1) )
       call redblack( grid(L)%u, grid(L)%f, A(L), gc(L), nu2)
    end do


  end subroutine vcycle
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


end module mg_mod
