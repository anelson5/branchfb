module fluidsolve_mod
  use grid_mod
  use param_mod
  use bobconst_mod
  use gc_mod
  use mg_mod
  use diffusion_mod
!  use gel_mod
  implicit none

  ! parameter to choose the type of projection method
  !   1 or 3
  ! the difference is in the handeling of the pressure
  !
  integer :: pmMeth = 1


  ! fluid variables -- all of these now have two ghost cells
  !
  real :: u(-1:nxb+2,-1:nyb+2, 2)
  real :: p(-1:nxb+2,-1:nyb+2)
  real :: fns(-1:nxb+2,-1:nyb+2, 2)
  real :: px(-1:nxb+2,-1:nyb+2), py(-1:nxb+2,-1:nyb+2)
  !real :: asq(-1:nxb+2,-1:nyb+2, 2)

  ! differences in nonlinear terms
  !
  real :: H1old(-1:nxb+2, -1:nyb+2), H2old(-1:nxb+2, -1:nyb+2)
  real :: H1(   -1:nxb+2, -1:nyb+2), H2(   -1:nxb+2, -1:nyb+2)

  ! projection variables  -- note that only ui,vi have 2 gc
  !
  real :: ui(-1:nxb+2,-1:nyb+2), vi(-1:nxb+2,-1:nyb+2)
  real :: phidt( 0:nxb+1, 0:nyb+1)
  real :: divu(  0:nxb+1, 0:nyb+1)
  real :: phidtx(0:nxb+1, 0:nyb+1), phidty(0:nxb+1,0:nyb+1)

  real :: temp1(-1:nxb+2, -1:nyb+2)
  real :: temp2(-1:nxb+2, -1:nyb+2)


  ! ghost cell managers for variables with boundary conditions
  !  and those which have operators applied to them
  !
  type (gcm) :: gc_phidt

  ! flag for the projection step is singular
  !
  logical :: isSing

  ! flag for first step
  !
  logical :: isFirst = .true.


!!$  ! used for boundary condtions on u,v
!!$  !
!!$  real :: ulf(1:nyb), urt(1:nyb), vlf(1:nyb), vrt(1:nyb)
!!$  real :: utp(1:nxb), ubt(1:nxb), vtp(1:nxb), vbt(1:nyb)


CONTAINS


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! SUBROUTINE: MOMENTUM - solve the momentum equation for the temporary
!                        velocity
!
subroutine momentum
  integer :: i,j

  if( pmMeth.eq.3 ) then
     px = 0.0
     py = 0.0
  end if

  
  ! copy the velocity to the intermediate velocity
  !
  ui = u(:,:,1)
  vi = u(:,:,2)


  ! solve for ui
  ! removed the Re in front of the adv terms 4/25/08 KL
  temp1 = fns(:,:,1) - px + (1.5*H1 - 0.5*H1old)
  !write(*,*) 'dt=', dt, 'temp1=', minval(H1old)

  call step_diffusion(ui,temp1,1.0,dt,1)
  
  temp1 = fns(:,:,2) - py + (1.5*H2 - 0.5*H2old)
  call step_diffusion(vi,temp1,1.0,dt,2)
!stop
end subroutine momentum
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! SUBROUTINE: PROJECT - perform the projection and update the velocity
!                       and pressure
subroutine project

  integer :: i,j

  ! short circuit if no projection
  !
  if( pmMeth.eq.0 ) then
     u(:,:,1) = ui
     u(:,:,2) = vi
     return
  end if

  ! GET RID OF THIS FOR NOW...
  ! set outflow edges to be periodic
  ! 
  !ui(nxb+1, 1:nyb) = ui(1, 1:nyb)


  ! compute the divergece of the intermediate cell edge velocity
  !
  divu(1:nxb,1:nyb) = ( ui(2:nxb+1,1:nyb) - ui(1:nxb,1:nyb) + &
                        vi(1:nxb,2:nyb+1) - vi(1:nxb,1:nyb) ) / hb

  ! if problem is singular, ensure that is is solvable 
  !
  if( isSing ) then
     divu = divu - sum( divu(1:nxb,1:nyb) ) /(1.0*nxb*nyb)
  end if

  ! solve the poisson problem for phidt 
  !
  call MG_solve( phidt, divu, 0.0,-1.0, gc_phidt)

  ! compute the gradient of phidt at the cell edges
  !
  call gc_op(phidt, gc_phidt)
  phidtx(1:nxb+1,1:nyb) = (phidt(1:nxb+1,1:nyb)-phidt(0:nxb,1:nyb)) / hb
  phidty(1:nxb,1:nyb+1) = (phidt(1:nxb,1:nyb+1)-phidt(1:nxb,0:nyb)) / hb

  ! update the cell edge velocity
  ! only copy what needs to be copied u(2:nxb+1,1:nyb), v(1:nxb,2:nyb) 
  !
  u(2:nxb+1,1:nyb,1) = ui(2:nxb+1,1:nyb) - phidtx(2:nxb+1,1:nyb)
  u(1:nxb,2:nyb,2) = vi(1:nxb,2:nyb) - phidty(1:nxb,2:nyb)

  ! update the pressure gradient
  !
  px(1:nxb+1,1:nyb) = px(1:nxb+1,1:nyb) + phidtx(1:nxb+1,1:nyb)/dt
  py(1:nxb,1:nyb+1) = py(1:nxb,1:nyb+1) + phidty(1:nxb,1:nyb+1)/dt

  ! update the pressure
  !
  p(1:nxb,1:nyb) = p(1:nxb,1:nyb) + phidt(1:nxb,1:nyb)/dt

end subroutine project
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! SUBROUTINE: SETNEWH -- uses the current velocity to compute the 
!                        nonlinear terms
!
subroutine setnewH
  real :: uu(-1:nxb+2,-1:nyb+2), uv(-1:nxb+2,-1:nyb+2), vv(-1:nxb+2,-1:nyb+2)

!!$  ! no nonlinear terms
!!$  !
!!$  H1 = 0.0
!!$  H2 = 0.0
!!$  return

  ! fill ghost cells
  !
  call fillgc

  ! do the differencing (uu)x
  !
!write(*,*)'u-u=', u(2,:,1) - u(0,:,1)
  uu = u(:,:,1)*u(:,:,1)
  H1(1:nxb+1,1:nyb) = -0.5*( uu(2:nxb+2,1:nyb) - uu(0:nxb,1:nyb) )/hb
!  write(*,*) 'maxval(uu)=', abs(uu(1,1:nyb) - uu(0,1:nyb))
!  write(*,*) 'u(0)=', u(0,1:nyb,1)
!  write(*,*) 'u(1) =', u(1,1:nyb,1)
!  write(*,*) 'u(2)=', u(2,1:nyb,1)
  ! form uv at the vertical edges
  !
  uv(1:nxb+1,1:nyb) = u(1:nxb+1,1:nyb,1) *                   &
        0.25*( u(1:nxb+1, 1:nyb,2) + u(1:nxb+1 ,2:nyb+1,2) +   &
               u(0:nxb  , 1:nyb,2) + u(0:nxb   ,2:nyb+1,2)  )
!  write(*,*) 'max uv=', maxval(abs(uv))
  ! do the differencing (uv)y
  !
  H1(1:nxb+1,1:nyb) = H1(1:nxb+1,1:nyb) - &
       0.5*( uv(1:nxb+1,2:nyb+1) - uv(1:nxb+1,0:nyb-1) )/hb
  
  ! compute H2 = -(uv)x - (vv)y
  !
  vv = u(:,:,2)*u(:,:,2)
  H2(1:nxb,1:nyb) =-0.5*(vv(1:nxb,2:nyb+1) - vv(1:nxb,0:nyb-1) )/hb
!write(*,*) 'maxval(vv)=', maxval(abs(vv))
  ! form uv at the horizontal edges
  !
  uv(0:nxb+1,1:nyb) = u(0:nxb+1,1:nyb,2) *                   &
       0.25*(u(0:nxb+1 ,1:nyb,1) + u(0:nxb+1, 0:nyb-1,1)   + &
             u(1:nxb+2 ,1:nyb,1) + u(1:nxb+2, 0:nyb-1,1) )

  H2(1:nxb,1:nyb) = H2(1:nxb,1:nyb) - &
       0.5*(uv(2:nxb+1,1:nyb)-uv(0:nxb-1,1:nyb))/hb
  !write(*,*) 'H1,H2=', maxval(H1), maxval(H2), minval(H1), minval(H2)
  H2 = 0.0
  H1 = 0.0
end subroutine setnewH
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! SUBROUTINE: FLUIDSTEP - take a regular (not first) timestep
!                        boundary condtions and forces should be
!                        set externally before calling this routine
!
subroutine fluidstep
  !real:: asqin(-1:nxb+2,-1:nyb+2,2)
!  write(*,*) 'fluidstep', isFirst
  ! do something else if this is the first step
  !
!  asq = asqin
  if( isFirst ) then
     call firstfluidstep
     return
  end if
!write(*,*) 'fluidstep2'
  H1old = H1
  H2old = H2
  !call alpha_sq(t) !CZ added this step here
  call setnewH
  call momentum
  call project

end subroutine fluidstep
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! SUBROUTINE: FIRSTFLUIDSTEP - take the first time step given and initial 
!                         velocity.
!                         boundary condtions and forces should be
!                         set externally before calling this routine
!
subroutine firstfluidstep
  integer :: i
  real :: uinit(-1:nxb+2,-1:nyb+2), vinit(-1:nxb+2,-1:nyb+2) 
!write(*,*) 'first fluid step'
  ! decided if the projection is singular or not
  !
  isSing = all( gc_phidt%type.ne.DIR ) 
!write(*,*) 'isSing=', isSing

  ! set all variables but the velocity to zero
  !
  p      = 0.0
  px     = 0.0
  py     = 0.0
  H1old  = 0.0
  H2old  = 0.0
  H1     = 0.0
  H2     = 0.0
  phidt  = 0.0
  ui     = 0.0
  vi     = 0.0
  divu   = 0.0
  phidtx = 0.0
  phidty = 0.0
  temp1  = 0.0
  temp2  = 0.0

  ! save the initial velocity
  !
  uinit = u(:,:,1)
  vinit = u(:,:,2)
  !write(*,*) 'uinit(0)=', uinit(0,:)-uinit(2,:)
  !write(*,*) 'u(0)-u(2)=', u(2,:,1) - u(0,:,1)
  ! compute the advective terms from the initial velocity
  !
  call setnewH
  H1old = H1
  H2old = H2
!  write(*,*) 'Through first H'
  ! iterate -- find new velocity, compute pressure, compute new adv. terms
  !
  do i=1,5
 !    write(*,*) 'initial iteration ', i
     u(:,:,1) = uinit
     u(:,:,2) = vinit
 !    write(*,*) 'u-uint=', u(0,:,1)-u(2,:,1)
     !write(*,*) 'maxu=', maxval(uinit)
     !write(*,*) 'maxv=', maxval(vinit)
     call momentum
     !write(*,*) 'asq_mom=', maxval(asq)
     
  !   write(*,*) 'momentum'
     call project

     call setnewH
     H1 = (2.0*H1old + H1) / 3.0  ! makes an average of two H's for
     H2 = (2.0*H2old + H2) / 3.0  !   momentum solve
     H1 = 0.0
     H2 = 0.0
  end do

  ! set flag for other steps
  !
  isFirst = .false.

  omegaf=1.8
  !write(*,*) 'p=', p(:,1:4)
end subroutine firstfluidstep
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! SUBROUTINE: FILLGC -- set the ghost cells of the velocity
!
!    HOMOGENEOUS DIRICHLET TOP AND BOTTOM
!    DIRICHLET LEFT NEUMANN RIGHT
!
subroutine fillgc
  integer :: j
  real :: y
  real :: utemp(1,1:nyb)

  ! horizontal velocity (u) top and bottom
  !
!  write(*,*) 'Pre fillgc =', u(0,:,1)
!  write(*,*) '2=', u(2,:,1)
  u(1:nxb+1,   0  , 1) = -u(1:nxb+1,  1,  1)
  u(1:nxb+1, nyb+1, 1) = -u(1:nxb+1, nyb, 1)


  ! vertical velocity (v) top and bottom
  !
  u(1:nxb,   0  , 2) = -u(1:nxb,  2 , 2)
  u(1:nxb, nyb+2, 2) = -u(1:nxb, nyb, 2)



  ! set inflow and outflow condition
  !
  do j=1,nyb
    y = ymin + hb*(j-0.5)
    ! use nondimensional inflow profile from driver.f90
    !utemp(1,j)=4.0*(y-ymin)*(ymax-y)
    !utemp(1,j)=umax/uchar*0.01*(y-ymin)*(ymax-y)
    utemp(1,j)=umax/(uchar*ymax**2)*4.0*(y-ymin)*(ymax - y)
 end do

  ! u left and right
  u(0, 1:nyb, 1) = 2.0*utemp(1,1:nyb) -u(2, 1:nyb, 1)
!  write(*,*)  2.0*utemp(1,1:nyb) - u(2, 1:nyb, 1)
!  write(*,*) utemp(1,1:nyb)-u(2,1:nyb,1)
!  write(*,*) u(2,1:nyb,1)
  u(nxb+2,1:nyb,1) = u(nxb, 1:nyb, 1)

  ! v left and right
  u(0, 1:nyb+1, 2) = -u(2, 1:nyb+1, 2)
  u(nxb+2,1:nyb+1,2) = u(nxb, 1:nyb+1, 2)
!  write(*,*) 'Post fillgc =', u(0,:,1)
!  write(*,*) '2=', u(2,:,1)

end subroutine fillgc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


end module fluidsolve_mod
