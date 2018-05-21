module grid_mod

! information about the grid
! 
integer, parameter :: maspect = 4            ! define the shape of the
integer, parameter :: qx      = 5            !5  base grid -- 2^qx*maspact X 2^qx

integer, parameter :: nyb = 2**(qx)            ! Ny for base grid
integer, parameter :: nxb = maspect* 2**(qx)   ! Nx for base grid

integer, parameter :: nyc = 1                ! Ny for coarse grid
integer, parameter :: nxc = maspect          ! Nx for coarse grid

real, parameter :: ymin = 0.0                        ! y at bottom edge
real, parameter :: ymax = 1.0 !dheight/xchar !1.0 !0.5!60.0/3.0 ! y at top edge dheight/xchar

real, parameter :: xmin = 0.0                        ! x at left edge
real, parameter :: xmax = xmin + maspect*(ymax-ymin) ! x at right edge

real, parameter :: hb   = (ymax-ymin)/nyb            ! cell width
real, parameter :: hsq  = hb*hb

integer, parameter :: jjmin = -1                      ! min and max indices for 
integer, parameter :: jjmax = nxb+3                   ! Aaron's 2d arrays
integer, parameter :: llmin = -1
integer, parameter :: llmax = nyb+3

integer, parameter :: numacross = jjmax-jjmin+1  ! if jjmin > 0 drop +1
integer, parameter :: numup     = llmax-llmin+1  ! if llmin > 0 drop +1

integer, parameter :: nyp1 = nyb + 1
integer, parameter :: nyp3 = nyb + 3
integer, parameter :: nxp1 = nxb + 1
integer, parameter :: nxp2 = nxb + 2
integer, parameter :: nxp3 = nxb + 3

CONTAINS

subroutine outputgridinfo( un )
  integer :: un
  write(*,*) 'nyb=', nyb
  write(un,'(a)') 'DOMAIN INFO'
  write(un,'(a)') '-----------'
  write(un,'(a,f6.2)') 'domain  xmin ', xmin
  write(un,'(a,f6.2)') '        xmax ', xmax
  write(un,'(a,f6.2)') '        ymin ', ymin
  write(un,'(a,f6.2)') '        ymax ', ymax
  write(un,'(a,i6)')   'grid size Nx ', nxb
  write(un,'(a,i6)')   '          Ny ', nyb
end subroutine outputgridinfo
 

end module grid_mod 
