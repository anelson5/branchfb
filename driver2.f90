
 module clot_mod
   use grid_mod
   use param_mod
   use bobconst_mod
   use gc_mod
   use diffusion_mod  
   use fluidsolve_mod

   implicit none

   ! constant background force on the fluid
   !
   real :: fbg(0:nxb+1,0:nyb+1,2)


   ! arrays to store the values of alpha^2 -- both sets of edges
   !
   real :: asq(1:nxb+1, 1:nyb+1, 2)


   ! tolerances
   !
   real :: mgtol   = 1e-10
   real :: difftol = 1e-7
  

 CONTAINS

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !
 ! SUBROUTINE: DRIVER --- main program driver
 !
 subroutine driver
   real    :: rk4tstep
   real    :: bcA(1:num_fp,1:nxb),bcB(1:num_fp,1:nxb)
   real    :: bcC(1:num_fp,1:nxb),bcCold(1:num_fp,1:nxb) ! for boundary conditions
   real    :: bcEta(1:nxb)
   real    :: pltsum(-1:nxb+2,-1:nyb+2)   ! for the density dependent advection
   real    :: pltsumh(-1:nxb+2,-1:nyb+2)  ! for the density dependent advection
   real    :: oldplt(1:nxb,1:nyb)         ! for coagreactions at each timestep
   real    :: y
   real    :: maxu                        ! for variable dt
   integer :: i,j,k                       ! for looping 
   integer :: outputflag,numsteps,numstepsm
   real    :: tnext,ttotal                ! for variable time
   real    :: dt_rxn                      ! time step for reactions
   integer :: mm,M                        ! for dt/M=dt_rxn
   real    :: thetap(1:nxb,1:nyb)         ! for hindered transport
   real    :: A_ones(-1:nxb+2,-1:nyb+2) ! coefficient for hindered ADP




   character(len=20) :: date
   character(len=20) :: time


   ! initialize everything
   !
   write(*,*) 'begin driver'
   call init
 !  call diffvc_init
 !  call diffvcplts_init

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! BEGIN SIMULATION 
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !  call date_and_time(date,time)
 !  write(*,'(a24,a2,a1,a2,a1,a4,a4,a2,a1,a2,a1,a2)')             &
 !       'Beginning simulation on ',date(5:6),'/',date(7:8),'/',date(1:4),' at ', &
 !       time(1:2),':',time(3:4),':',time(5:6)  


   ! set output flag to zero
   outputflag = 0

   ! set time (s) increments for output
   tnext = 10.0

   ! set total time passed (s) to zero
   ttotal = 0.0

   numsteps=0
   numstepsm=0

   M = 0

   ! set bcs for eta
   bcEta=0.0

   ! SET TIMESTEP FOR DIFFUSION AND REACTION
   dt_rxn  = 1.0e-3/tchar   ! 1 msec

   thetap(:,:)=0.0
   !for ADP
   A_ones(:,:)=1.0

   write(*,*)'time step for diff and rxn (sec) : ',dt_rxn*tchar

   ! set rk timestep to dt_rxn
   rk4tstep=dt_rxn

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! BEGIN TIME LOOP
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   write(*,*) 'begin time loop'
   do 

      ! maximum u
      maxu = maxval( abs(u) )

      ! figure out how many times to run through fluid and advection

      ! for testing:
      ! dt_rxn = 0.75*hb/maxu	
      ! M = 1

      M = ceiling(maxu*dt_rxn/(0.75*hb))
      if(numsteps.eq.0)then
         write(*,*)'M',M
      end if


      ! use 3/4 for constant CFL condition 
      dt = dt_rxn/M

      ! UPDATE TOTAL TIME 
      ttotal=ttotal+dt_rxn*tchar
      write(*,*) 'ttotal =', ttotal
     ! set new alpha in fluid solver
      call new_alpha(asq)

     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! ANYTHING NEEDED FOR FLUID AND ADVECTION DO HERE:
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      do mm=1,M
      ! set the background pressure gradient
      !
      fbg(:,:,1) = 0.0
      fbg(:,:,2) = 0.0

      ! set the background force - fns initialezed in fluidinit
      fns(1:nxb+1,1:nyb+1,:) = fbg(1:nxb+1,1:nyb+1,:)



      ! step the fluid solver
      call fluidstep


      tcount=tcount+dt;

      ! update sigmar and sigmaz
      end do
   if (ttotal.gt.tnext)then
      outputflag = 1
   end if
   if (ttotal.ge.tfinal) then
      write(*,*) 'Tfinal', tfinal
      write(*,*) 'Reached time', ttotal
      exit
   end if
end do  
end subroutine driver
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! SUBROUTINE: ALPHA_SQ --- evaluate alpha^2 function at a point
!
subroutine alpha_sq(p1,p2)
  real :: p1(-1:nxb+2,-1:nyb+2) ! bound platelets
  real :: p2(-1:nxb+2,-1:nyb+2) ! subendothelial platelets
  real :: volfrac(-1:nxb+2,-1:nyb+2)
  real :: Aalpha(1:nxb,1:nyb)   ! A = asqmax*volfrac.^6./(0.01+volfrac.^6)
  real :: Calpha(1:nxb-1,1:nyb-1) ! centers, see notes
  real :: ndp,nds
  integer :: i,j   

  ndp = Pmaxb*6.022e17
  nds = Pmaxse*6.022e17

  volfrac = 0.0

  ! set volume fraction and A  

  do i=1,nxb
     do j=1,nyb
        volfrac(i,j) = p1(i,j)+p2(i,j)
        Aalpha(i,j)  = asqmax*volfrac(i,j)**2/(0.5**2+volfrac(i,j)**2)
     end do
  end do

  ! set centers

  do i=1,nxb-1
     do j=1,nyb-1
         Calpha(i,j)=0.25*(Aalpha(i,j)+Aalpha(i+1,j)+Aalpha(i,j+1)+Aalpha(i+1,j+1));
     end do
  end do

  ! set alpha edges

  asq(1,1:nyb,1)=Aalpha(1,1:nyb);
  asq(nxb+1,1:nyb,1)=Aalpha(nxb,1:nyb);
  asq(2:nxb,1,1)=(Aalpha(1:nxb-1,1)+Aalpha(2:nxb,1))/2;
  asq(2:nxb,nyb,1)=(Aalpha(1:nxb-1,nyb)+Aalpha(2:nxb,nyb))/2;

  asq(1:nxb,1,1)=Aalpha(1:nxb,1);
  asq(1:nxb,nyb+1,2)=Aalpha(1:nxb,nyb);
  asq(1,2:nyb,2)=(Aalpha(1,1:nyb-1)+Aalpha(1,2:nyb))/2;
  asq(nxb,2:nyb,2)=(Aalpha(nxb,1:nyb-1)+Aalpha(nxb,2:nyb))/2;


  ! set alpha centers

  do i=2,nxb
     do j=2,nyb-1
            asq(i,j,1)=0.25*(Calpha(i-1,j-1) + Calpha(i-1,j) + &
                 Aalpha(i,j) + Aalpha(i-1,j));
     end do
  end do

  do i=2,nxb-1
     do j=2,nyb
            asq(i,j,2)=0.25*(Calpha(i-1,j-1) + Calpha(i,j-1) + &
                 Aalpha(i,j) + Aalpha(i,j-1));
     end do
  end do
  

end subroutine alpha_sq
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!I need to add a different type of alpha_sq file to reflect fibers !!!
!!!permeability information!!!!


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
write(*,*) 'Initialize1' 
  call loadparams
 
  ! put the rest of the data in the Data directory
  !
  runname = './Data/'//runname
  strlen = trimbt(runname)


  ! set the parameters
  !
write(*,*) 'Initialize2'
  call setparams

  ! initialize fluid variables
  !
write(*,*) 'Initialize3'
  call fluidinit

   
  ! initialize biochemical variable
  !
!  call biocheminit


  ! output all parameters and information about the run
  !

!  write(filename,'(2a)') runname(1:strlen),'.info'
!  open(30,file=filename)
!  write(30,'(a)') 'INFO FOR CLUMP SIMULATION'
!  write(30,'(a)') '---------------------------------'
!  call outputstepinfo(30)!
!  call outputparams(30)
  write(*,*)' init after outputparams'
!  write(30,*)
!  write(30,*)
!  call outputgridinfo(30)
!  close(30)


  ! output the initial data
  !
    
!  call output

end subroutine init
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  SUBROUTINE: FLUIDINIT --- initialize the fluid solver  
!
subroutine fluidinit
  real :: x,y
  integer :: i,j

  ! set alpha on vertical edges ( u )
  !
        
  asq(i,j,1) = 0.0

  !  set alpha on horizontal edges ( v )
  !
       
  asq(i,j,2) = 0.0

  ! initialize the diffusion solver
  !
  call diff_init( asq )

  
  ! initialize ghost cell managers for fluid variables
  ! note that gc_phidt refers to the ghost cells for
  ! the projection. these are related to the actual bcs 

  call gc_init(gc_phidt,nxb,nyb,hb,0.0,-1.0,NMN,DIR,NMN,NMN)

  ! set the tolerance of the projection and momentum solver
  !
  call MG_settol( mgtol )
  call diff_settol( difftol )


  ! initialize with a Poiseuille flow with max vel as umax
  ! zero outside the domain
 
 !ADDED THIS FOR TESTING CONVERGENCE
  u(:,:,1)=0.0

  do j=1,nyb
     y = ymin + hb*(j-0.5)
     u(:,j,1) = umax/(ymax**2)*4.0*(y-ymin)*(ymax - y) !this is dimensional!
     ! note using the 0.01 instead of 4 since the ymax = 20 and not 1
  end do


  ! set the vertical velocity and pressure to zero
  !
  u(:,:,2) = 0.0
  p(:,:) = 0.0


  ! nondimensionalize the velocities
  u(:,:,1) = u(:,:,1)/uchar
  

end subroutine fluidinit
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! SUBROUTINE: OUTPUT --- output data to files
!
subroutine output
  integer       :: strlen      ! length of runname string   
  integer       :: namelen     ! length of chem name string
  character(40) :: velfile     ! name of velocity file
  character(40) :: chemfile    ! name of chemical file
  character(40) :: filetype    ! binary or formatted
  integer       :: i,j,k       ! loop counters
  
  ! set the file type
  ! 
  if( isBinary ) then
     filetype = 'binary'
  else
     filetype = 'formatted'
  end if
  
  ! define file names
  !
  strlen = trimbt(runname)
  write(velfile,'(2a,i4.4,a,i4.4)') runname(1:strlen),'.',nyb,'.vel.',outcount

  ! open files for writing
  !
  open(20,file=velfile,form=filetype)

  ! for outputing the veloicties, we include the extra cell, which
  !  are not needed for all variables.  
  !

     do i=1,nxb+1
        do j=1,nyb+1
	  ! put dimensions of velocity back in
           write(20,'(5es14.5e3)') uchar*u(i,j,1),uchar*u(i,j,2),p(i,j),asq(i,j,1),asq(i,j,2)
        end do
     end do
  ! close the file
  !
  close(20)


 do i=1,num_fp
    
    namelen = len_trim(fp(i)%name)
    !write(*,*)'fp i',i,namelen,fp(i)%name
    write(chemfile,'(4a,i4.4,a,i4.4,a)') runname(1:strlen),'.' &
	,fp(i)%name(3:namelen),'.',nyb,'.',outcount,'.dat'
    open(50,file=chemfile,form=filetype)

          do j=1,nxb
             do k=1,nyb
                write(50,'(3es14.5e3)') fp(i)%value(j,k)*ndfp(i)
             end do
          end do
       close(50)   
 end do

do i=1,num_pl    
    namelen = len_trim(pl(i)%name)  
    !write(*,*)'pl i',i,namelen,pl(i)%name
    write(chemfile,'(4a,i4.4,a,i4.4,a)') runname(1:strlen),'.' &
	,pl(i)%name(3:namelen),'.',nyb,'.',outcount,'.dat'
    open(70,file=chemfile,form=filetype)
          do j=1,nxb
             do k=1,nyb
                write(70,'(3es14.5e3)') pl(i)%value(j,k)*ndpl(i)
             end do
          end do
       close(70)

 end do


 do i=1,num_pb
    namelen = len_trim(pb(i)%name)
    !write(*,*)'pb i',i,namelen,pb(i)%name
    write(chemfile,'(4a,i4.4,a,i4.4,a)') runname(1:strlen),'.' &
	,pb(i)%name(3:namelen),'.',nyb,'.',outcount,'.dat'
    open(61,file=chemfile,form=filetype)
          do j=1,nxb
             do k=1,nyb
              write(61,'(3es14.5e3)') pb(i)%value(j,k)*ndpb(i)   
             end do
          end do
    close(61)
 end do

 do i=1,num_se   
    namelen = len_trim(se(i)%name)
    !write(*,*)'se i',i,namelen,se(i)%name
    write(chemfile,'(4a,i4.4,a,i4.4,a)') runname(1:strlen),'.' &
	,se(i)%name(3:namelen),'.',nyb,'.',outcount,'.dat'
    open(60,file=chemfile,form=filetype) 


          do j=1,nxb
          !if too small, just make zero
          if(abs(se(i)%value(j)).lt.1.0e-98)then 
             write(60,'(3es14.5e3)') 0.0
          else
             write(60,'(3es14.5e3)') se(i)%value(j)*ndse(i)    
          end if

          end do
    close(60)
 end do

  do i=1,1
         
     namelen = len_trim(fp(i)%name)
     !write(*,*)'fp i',i,namelen,fp(i)%name
     write(chemfile,'(4a,i4.4,a,i4.4,a)') runname(1:strlen),'.' &
         ,fp(i)%name(3:namelen),'.',nyb,'.diff.',outcount,'.dat' 
     open(51,file=chemfile,form=filetype)
           
           do j=1,nxb
              do k=1,nyb
                 write(51,'(3es14.5e3)') D_h(j,k),A_h(j,k)
              end do
           end do
        close(51)
  end do







! DON'T FORGET TO ADD OUTPUT FOR PB TOO


  ! update the output counter
  !
  outcount = outcount + 1

end subroutine output
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


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




end module clot_mod


!--------------------------------------------------------------------------
!--------------------------------------------------------------------------



! FINALLY, THE MAIN PROGRAM
!
program clot
  use clot_mod
write(*,*) 'start driver'
  call driver
write(*,*) 'end driver'
end program clot

