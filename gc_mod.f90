module gc_mod
  use bobconst_mod
  use grid_mod
  PRIVATE

  type gcm
     integer :: nx, ny    ! record the grid size
     real    :: h         ! record the grid spacing
     integer :: type(4)   ! give the type of boundary condition

     
     real, pointer :: lf_dat(:)  ! used for storing nonhomogeneous data
     real, pointer :: rt_dat(:)  !  these are only allocated for 
     real, pointer :: bt_dat(:)  !  dirichlet and neumann boundaries
     real, pointer :: tp_dat(:)

     real    :: lf_rel(4) ! weights for ghost cells for relaxing 
     real    :: rt_rel(4) !  1 - weights in +/- parallel direction
     real    :: bt_rel(4) !  2 - weight  in normal direction
     real    :: tp_rel(4) !  3 - weight of rhs
                          !  4 - weight of bndy data

     real    :: sw_rel(5) ! weights for corner ghost cells for relaxing
     real    :: se_rel(5) !  1 - weight normal to top and bottom
     real    :: nw_rel(5) !  2 - weight normal to left and right
     real    :: ne_rel(5) !  3 - weight of rhs
                          !  4 - weight of data on left and right
                          !  5 - weight of data on top and bottom

     integer :: lf_int    ! weights for setting interpolation 
     integer :: rt_int    !  currently this only applies for
     integer :: bt_int    !  homogeneous boundaries
     integer :: tp_int
     integer :: sw_int
     integer :: se_int
     integer :: nw_int
     integer :: ne_int

     logical :: isAllocated = .false.
  end type gcm

  ! overload the =
  !
  interface assignment(=)
     module procedure set_equal_op
  end interface

  ! List of public routines and types
  !
  public :: gc_init, gc_relax, gc_op, gc_opx,gc_opy, gc_interp
  public :: gc_edgex, gc_edgey
  public :: gcm



CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! SUBROUTINE - GC_INIT: sets the size, boundary type, problem type
!                       a and b refer to the problem aI-bL                    
!                       
!
  subroutine gc_init(gc,nx,ny,h,a,b,lf_type,rt_type,bt_type,tp_type)
    type (gcm) :: gc
    integer :: nx, ny
    real    :: h
    real    :: a,b
    integer :: lf_type,rt_type,bt_type,tp_type


    ! Make sure that if one side is periodic, 
    !  the other side is also
    !
    if( ((lf_type.eq.PER).and.(rt_type.ne.PER)).or. &
        ((lf_type.ne.PER).and.(rt_type.eq.PER)).or. & 
        ((bt_type.ne.PER).and.(tp_type.eq.PER)).or. &
        ((bt_type.eq.PER).and.(tp_type.ne.PER))     ) then
       write(*,*) "Invalid periodic boundary conditions!"
       stop
    end if

    ! set the boundary types
    !
    gc%type(LFT) = lf_type
    gc%type(RGT) = rt_type
    gc%type(BOT) = bt_type
    gc%type(TOP) = tp_type


    ! record the grid size
    !
    gc%nx = nx
    gc%ny = ny
    gc%h  = h

    ! allocate boundary data for dirichlet and neumann edges
    !    
    if( .not.gc%isAllocated ) then
       if( (lf_type.eq.DIR).or.(lf_type.eq.NMN) ) then
          allocate( gc%lf_dat(ny) )
          gc%lf_dat(:) = 0.0 
       end if
       
       if( (rt_type.eq.DIR).or.(rt_type.eq.NMN) ) then
          allocate( gc%rt_dat(ny) )
          gc%rt_dat(:) = 0.0 
       end if
       
       if( (bt_type.eq.DIR).or.(bt_type.eq.NMN) ) then
          allocate( gc%bt_dat(nx) )
          gc%bt_dat(:) = 0.0 
       end if
       
       if( (tp_type.eq.DIR).or.(tp_type.eq.NMN) ) then
          allocate( gc%tp_dat(nx) )
          gc%tp_dat(:) = 0.0 
       end if

       gc%isAllocated = .true.
    end if


    ! set weights for the left edge and the two of the corners
    !
    if( lf_type.eq.DIR ) then

       ! edge weights for dirichlet left edge
       !
       call set_dir(gc%lf_rel,a,b,h)
       gc%lf_int = -1

       ! southwest corner - for dirichlet on left
       !
       if( bt_type.eq.DIR ) then
          call set_ddc(gc%sw_rel,a,b,h)
          gc%sw_int = 1
       else if( bt_type.eq.NMN ) then
          call set_dnc(gc%sw_rel,a,b,h)
          gc%sw_int = -1
       end if
       
       ! northwest corner - for dirichlet on left
       !
       if( tp_type.eq.DIR ) then
          call set_ddc(gc%nw_rel,a,b,h)
          gc%nw_int = 1
       else if ( tp_type.eq.NMN ) then
          call set_dnc(gc%nw_rel,a,b,h)
          gc%nw_int = -1
       end if

    else if( lf_type.eq.NMN ) then

       ! edge weights for neumann left edge
       !
       call set_nmn(gc%lf_rel,a,b,h)
       gc%lf_int =  1

       ! southwest corner - for neumann on left
       !
       if( bt_type.eq.DIR ) then
          call set_ndc(gc%sw_rel,a,b,h)
          gc%sw_int = -1
       else if( bt_type.eq.NMN ) then
          call set_nnc(gc%sw_rel,a,b,h)
          gc%sw_int = 1
       end if
       
       ! northwest corner - for neumann on left
       !
       if( tp_type.eq.DIR ) then
          call set_ndc(gc%nw_rel,a,b,h)
          gc%nw_int = -1
       else if ( tp_type.eq.NMN ) then
          call set_nnc(gc%nw_rel,a,b,h)
          gc%nw_int = 1
       end if
    end if



    ! set weights for the right edge and two of the corners
    !
    if( rt_type.eq.DIR ) then

       ! edge weights for dirichlet right edge
       !
       call set_dir(gc%rt_rel,a,b,h)
       gc%rt_int = -1

       ! southeast corner - for dirichlet on right
       !
       if( bt_type.eq.DIR ) then
          call set_ddc(gc%se_rel,a,b,h)
          gc%se_int = 1
       else if( bt_type.eq.NMN ) then
          call set_dnc(gc%se_rel,a,b,h)
          gc%se_int = -1
       end if
       
       ! northeast corner - for dirichlet on right
       !
       if( tp_type.eq.DIR ) then
          call set_ddc(gc%ne_rel,a,b,h)
          gc%ne_int = 1
       else if ( tp_type.eq.NMN ) then
          call set_dnc(gc%ne_rel,a,b,h)
          gc%ne_int = -1
       end if

    else if ( rt_type.eq.NMN ) then

       ! edge weights for neumann right edge
       !
       call set_nmn(gc%rt_rel,a,b,h)
       gc%rt_int =  1

       ! southeast corner - for neumann on right
       !
       if( bt_type.eq.DIR ) then
          call set_ndc(gc%se_rel,a,b,h)
          gc%se_int = -1
       else if( bt_type.eq.NMN ) then
          call set_nnc(gc%se_rel,a,b,h)
          gc%se_int = 1
       end if
       
       ! northeast corner - for neumann on right
       !
       if( tp_type.eq.DIR ) then
          call set_ndc(gc%ne_rel,a,b,h)
          gc%ne_int = -1
       else if ( tp_type.eq.NMN ) then
          call set_nnc(gc%ne_rel,a,b,h)
          gc%ne_int = 1
       end if

    end if


    ! set weights for the bottom edge
    !
    if( bt_type.eq.DIR ) then
       call set_dir(gc%bt_rel,a,b,h)
       gc%bt_int = -1
    else if ( bt_type.eq.NMN ) then
       call set_nmn(gc%bt_rel,a,b,h)
       gc%bt_int =  1
    end if

    ! set weights for the top edge
    !
    if( tp_type.eq.DIR ) then
       call set_dir(gc%tp_rel,a,b,h)
       gc%tp_int = -1
    else if ( tp_type.eq.NMN ) then
       call set_nmn(gc%tp_rel,a,b,h)
       gc%tp_int =  1
    end if

  end subroutine gc_init
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! SUBROUTINE: GC_OP - set the ghost cells for applying an operator
!
subroutine gc_op(u,gc)
  real :: u(0:,0:)
  type (gcm) :: gc

  call gc_opx(u,gc)
  call gc_opy(u,gc)

end subroutine gc_op

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! SUBROUTINE: GC_OPX - set the ghost cells for applying an operator
!                      only in x-direction
!
subroutine gc_opx(u,gc)
  real :: u(0:,0:)
  type (gcm) :: gc

  if( gc%type(LFT).eq.PER ) then
     call set_x_per(u,gc)
  else
     call set( u(0,1:gc%ny), u(1,1:gc%ny), &
               u(2,1:gc%ny), u(3,1:gc%ny), &
               u(4,1:gc%ny),               &
               gc%lf_dat(:), gc%h, gc%type(LFT) )
     call set( u(gc%nx+1,1:gc%ny), u(gc%nx,1:gc%ny),   &
               u(gc%nx-1,1:gc%ny), u(gc%nx-2,1:gc%ny), &
               u(gc%nx-3,1:gc%ny),                     &
               gc%rt_dat(:), gc%h, gc%type(RGT) )
  end if

end subroutine gc_opx

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! SUBROUTINE: GC_OPY - set the ghost cells for applying an operator
!                      only in the y-direction
!
subroutine gc_opy(u,gc)
  real :: u(0:,0:)
  type (gcm) :: gc

  if( gc%type(BOT).eq.PER ) then
     call set_y_per(u,gc)
  else
     call set( u(1:gc%nx,0), u(1:gc%nx,1), &
               u(1:gc%nx,2), u(1:gc%nx,3), &
               u(1:gc%nx,4),               &
               gc%bt_dat, gc%h, gc%type(BOT) )
     call set( u(1:gc%nx,gc%ny+1), u(1:gc%nx,gc%ny),   &
               u(1:gc%nx,gc%ny-1), u(1:gc%nx,gc%ny-2), &
               u(1:gc%nx,gc%ny-3),                     &
               gc%tp_dat, gc%h, gc%type(TOP) )
  end if


end subroutine gc_opy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! SUBROUTINE - GC_RELAX
!
subroutine gc_relax(u,f,gc)
  real :: u(0:,0:), f(0:,0:)
  type (gcm) :: gc

  integer :: i,j
  integer :: nx, ny

  nx = gc%nx
  ny = gc%ny

  ! set periodic edges first
  !
  if( gc%type(LFT).eq.PER ) then
     call set_x_per(u,gc)
  end if

  if( gc%type(BOT).eq.PER ) then
     call set_y_per(u,gc)
  end if


  ! set the non-periodic edges
  !
  if( gc%type(LFT).ne.PER ) then
     u(0,1:ny) = gc%lf_rel(1)*u(1,2:ny+1) + gc%lf_rel(1)*u(1,0:ny-1) + &
                 gc%lf_rel(2)*u(2,1:ny  ) + gc%lf_rel(3)*f(1,1:ny  ) + &
                 gc%lf_rel(4)*gc%lf_dat(1:ny) 

     u(nx+1,1:ny) = gc%rt_rel(1)*u(nx  ,2:ny+1) + gc%rt_rel(1)*u(nx,0:ny-1) + &
                    gc%rt_rel(2)*u(nx-1,1:ny  ) + gc%rt_rel(3)*f(nx,1:ny  ) + &
                    gc%rt_rel(4)*gc%rt_dat(1:ny) 
  end if
  
  if( gc%type(BOT).ne.PER ) then
     u(1:nx,0) = gc%bt_rel(1)*u(0:nx-1,1) + gc%bt_rel(1)*u(2:nx+1,1) + &
                 gc%bt_rel(2)*u(1:nx  ,2) + gc%bt_rel(3)*f(1:nx  ,1) + &
                 gc%bt_rel(4)*gc%bt_dat(1:nx) 
     
     u(1:nx,ny+1) = gc%tp_rel(1)*u(0:nx-1,ny) + gc%tp_rel(1)*u(2:nx+1,ny) + &
                    gc%tp_rel(2)*u(1:nx  ,ny-1) + gc%tp_rel(3)*f(1:nx,ny) + &
                    gc%tp_rel(4)*gc%tp_dat(1:nx)
  end if


  ! if no edges are periodic, set the corners too
  !
  if( (gc%type(LFT).ne.PER).and.(gc%type(BOT).ne.PER) ) then
     u(0   ,  1) = gc%sw_rel(1)*u(1,2) + gc%sw_rel(2)*u(2,1)    + &
                   gc%sw_rel(3)*f(1,1) + gc%sw_rel(4)*gc%lf_dat(1) + &
                   gc%sw_rel(5)*gc%bt_dat(1) 

     u(nx+1,  1) = gc%se_rel(1)*u(nx,2) + gc%se_rel(2)*u(nx-1,1) + &
                   gc%se_rel(3)*f(nx,1) + gc%se_rel(4)*gc%rt_dat(1) + &
                   gc%se_rel(5)*gc%bt_dat(nx)

     u(0   , ny) = gc%nw_rel(1)*u(1,ny-1) + gc%nw_rel(2)*u(2,ny)    + &
                   gc%nw_rel(3)*f(1,ny  ) + gc%nw_rel(4)*gc%lf_dat(ny) + &
                   gc%nw_rel(5)*gc%tp_dat(1)

     u(nx+1, ny) = gc%ne_rel(1)*u(nx,ny-1) + gc%ne_rel(2)*u(nx-1,ny) + &
                   gc%ne_rel(3)*f(nx,ny  ) + gc%ne_rel(4)*gc%rt_dat(ny) + &
                   gc%ne_rel(5)*gc%tp_dat(nx)

     u(1 , 0   ) = 0.0 
     u(nx, 0   ) = 0.0 
     u(1 , ny+1) = 0.0 
     u(nx, ny+1) = 0.0      
  end if

end subroutine gc_relax
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! SUBROUTINE - GC_INTERP
!
subroutine gc_interp(u,gc)
  real :: u(0:,0:)
  type (gcm) :: gc

  integer :: i,j
  integer :: nx,ny

  nx = gc%nx
  ny = gc%ny


  ! set periodic edges first
  ! 
  if ( gc%type(LFT).eq.PER ) then
     call set_x_per(u,gc)
  end if
  if( gc%type(BOT).eq.PER ) then
     call set_y_per(u,gc)
  end if

  ! set the horizonal edges -- including the corners
  !
  if ( gc%type(LFT).ne.PER ) then
     u(0   ,0:ny+1) = gc%lf_int*u(1 ,0:ny+1)
     u(nx+1,0:ny+1) = gc%rt_int*u(nx,0:ny+1)
  end if

  ! set the vertical edges -- including the corners
  !
  if( gc%type(BOT).ne.PER ) then
     u(0:nx+1,0   ) = gc%bt_int*u(0:nx+1, 1)
     u(0:nx+1,ny+1) = gc%tp_int*u(0:nx+1,ny)
  end if

  ! set the corners if needed
  !
  if( (gc%type(LFT).eq.PER).and.(gc%type(BOT).eq.PER) ) then
     u(0   ,0   ) = u(nx,ny)
     u(0   ,ny+1) = u(nx, 1)
     u(nx+1,0   ) = u(1 ,ny)
     u(nx+1,ny+1) = u(1 , 1)
  else if( (gc%type(LFT).ne.PER).and.(gc%type(BOT).ne.PER) ) then
     u(0   ,0   ) = gc%sw_int*u(1 ,1 )
     u(0   ,ny+1) = gc%nw_int*u(1 ,ny)
     u(nx+1,0   ) = gc%se_int*u(nx,1 )
     u(nx+1,ny+1) = gc%ne_int*u(nx,ny)
  end if

end subroutine gc_interp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! SUBROUTINE: SET - used for gc_op
!
subroutine set(u0,u1,u2,u3,u4,dat,h,type)
  real :: u0(:), u1(:), u2(:), u3(:), u4(:),dat(:)
  real :: h
  integer :: type

  select case (type)
  case (DIR)
     u0 = -2.0*u1 + (u2 + 8.0*dat)/3.0 
  case (NMN)
     u0 = u1 - h*dat
  case (EXT)
     u0 = 2.0*(u1-u3)+u4
  end select

end subroutine set
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! SUBROUTINE: SET2 - used for computing edge velocities
!               the only difference from set is linear extrapolation
!               is used for Dirichlet edges
!
subroutine set2(u0,u1,dat,h,type)
  real :: u0(:), u1(:),dat(:)
  real :: h
  integer :: type

  select case (type)
  case (DIR)
     u0 = 2.0*dat - u1 
  case (NMN)
     u0 = u1 - h*dat
  case default
     write(*,*) 'Illegal velocity boundary condition'
     stop
  end select

end subroutine set2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! SUBROUTINE: GC_EDGEX
!
subroutine gc_edgex(u,gc)
  real :: u(0:,0:)
  type (gcm) :: gc

  if( gc%type(LFT).eq.PER ) then
     call set_x_per(u,gc)
  else
     call set2(u(0,1:gc%ny),u(1,1:gc%ny),  &
               gc%lf_dat(:), gc%h, gc%type(LFT))
     call set2(u(gc%nx+1,1:gc%ny),u(gc%nx,1:gc%ny), &
               gc%rt_dat(:), gc%h, gc%type(RGT))
  end if
end subroutine gc_edgex
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! SUBROUTINE: GC_EDGEY
!
subroutine gc_edgey(v,gc)
  real :: v(0:,0:)
  type (gcm) :: gc
  
  if( gc%type(BOT).eq.PER ) then
     call set_y_per(v,gc)
  else
     call set2(v(1:gc%nx,0),v(1:gc%nx,1),  &
               gc%bt_dat(:), gc%h, gc%type(BOT) )
     call set2(v(1:gc%nx,gc%ny+1), v(1:gc%nx,gc%ny),  &
               gc%tp_dat(:), gc%h, gc%type(TOP) )
  end if
end subroutine gc_edgey
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!*
!*
!* Routines for setting periodic edges
!*   all following routines are private
!*


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! SUBROUTINE - SET_X_PER: set the ghost cells for periodic domains in
!                         the x-direction
!
subroutine set_x_per(u,gc)
  real :: u(0:,0:)
  type (gcm) :: gc

  u(0      , 1:gc%ny) = u(gc%nx, 1:gc%ny)
  u(gc%nx+1, 1:gc%ny) = u(1    , 1:gc%ny)
end subroutine set_x_per
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! SUBROUTINE - SET_Y_PER: set the ghost cells for periodic domains in
!                           the y-direction
!
subroutine set_y_per(u,gc)
  real :: u(0:,0:)
  type (gcm) :: gc
  
  u(1:gc%nx, 0      ) = u(1:gc%nx, gc%ny)
  u(1:gc%nx, gc%ny+1) = u(1:gc%nx, 1    )

end subroutine set_y_per
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!*
!*
!* Routines to set the weights
!*
!*


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! SUBROUTINE - SET_DIR: sets the weights for a dirichlet edge
!
subroutine set_dir(wgt,a,b,h)
  real :: wgt(4)
  real :: a, b, h

  real :: ah2 
  ah2 = a*h*h

  wgt = (/-2.0*b, (ah2-2.0*b)/3.0, -2.0, (8.0*ah2+32.0*b)/3.0 /)
  wgt = wgt / (ah2+6.0*b)

end subroutine set_dir
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! SUBROUTINE - SET_NMN: sets the weighs for a neumann edge
!
subroutine set_nmn(wgt,a,b,h)
  real :: wgt(4)
  real :: a, b, h

  real :: ah2 
  ah2 = a*h*h

  wgt = (/b,b,1.0,(-ah2-4.0*b)*h/)
  wgt = wgt / (ah2+3.0*b)

end subroutine set_nmn
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! SUBROUTINE - SET_DNC: sets the weights for a dirichlet--neumann corner
!                       first type is always the vertical boundary
!
subroutine set_dnc(wgtc,a,b,h)
  real :: wgtc(5)
  real :: a,b,h

  real :: ah2 
  ah2 = a*h*h

  wgtc = (/-b, (ah2+b)/3.0, -1.0, (8.0*ah2+32.0*b)/3.0, -(ah2+4.0*b)*h/)
  wgtc = wgtc / (ah2 + 5.0*b)

end subroutine set_dnc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! SUBROUTINE - SET_NDC: sets the weights for a neumann--dirichlet corner
!                       first type is always the vertical boundary
!
subroutine set_ndc(wgtc,a,b,h)
  real :: wgtc(5)
  real :: a,b,h

  real :: ah2 
  ah2 = a*h*h


  wgtc = (/ (ah2+b)/3.0, -b, -1.0, -(ah2+4.0*b)*h, (8.0*ah2+32.0*b)/3.0 /)
  wgtc = wgtc / (ah2 + 5.0*b)

end subroutine set_ndc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! SUBROUTINE - SET_DDC: sets the weights for dirichlet--dirichlet corner
!
subroutine set_ddc(wgtc,a,b,h)
  real :: wgtc(5)
  real :: a,b,h

  real :: ah2 
  ah2 = a*h*h

  wgtc = (/ (ah2-8.0*b)/3.0, (ah2-8.0*b)/3.0, -4.0, (8.0*ah2+32.0*b)/3.0, &
            (8.0*ah2+32.0*b)/3.0 /)
  wgtc = wgtc / (ah2 + 8.0*b)

end subroutine set_ddc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! SUBROUTINE - SET_NNC: sets the weights for a neumann--neumann corner
!
subroutine set_nnc(wgtc,a,b,h)
  real :: wgtc(5)
  real :: a, b, h

  real :: ah2 
  ah2 = a*h*h

  wgtc = (/2.0*b, 2.0*b, 2.0, -(ah2+4.0*b)*h, -(ah2+4.0*b)*h/)
  wgtc = wgtc / (ah2 + 2.0*b)

end subroutine set_nnc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! SUBROUTINE: SET_EQUAL_OP - set opererators equal, used for overload
!
subroutine set_equal_op(g1,g2)
  type (gcm), intent(out) :: g1
  type (gcm), intent(in)  :: g2
  
  g1%nx = g2%nx
  g1%ny = g2%ny
  g1%h  = g2%h
  g1%type(:) = g2%type(:)

  g1%lf_dat = g2%lf_dat
  g1%rt_dat = g2%rt_dat
  g1%bt_dat = g2%bt_dat
  g1%tp_dat = g2%tp_dat
 
  g1%lf_rel(:) = g2%lf_rel(:)
  g1%rt_rel(:) = g2%rt_rel(:)
  g1%bt_rel(:) = g2%bt_rel(:)
  g1%tp_rel(:) = g2%tp_rel(:)

  g1%sw_rel(:) = g2%sw_rel(:)
  g1%se_rel(:) = g2%se_rel(:)
  g1%nw_rel(:) = g2%nw_rel(:)
  g1%ne_rel(:) = g2%ne_rel(:)

 
  g1%lf_int = g2%lf_int
  g1%rt_int = g2%rt_int
  g1%bt_int = g2%bt_int
  g1%tp_int = g2%tp_int
  g1%sw_int = g2%sw_int
  g1%se_int = g2%se_int
  g1%nw_int = g2%nw_int
  g1%ne_int = g2%ne_int
   
end subroutine set_equal_op

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


end module gc_mod
