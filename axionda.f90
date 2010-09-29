

program axionda

! Solve wave equation in effective 2D = axisymmetric 3D grid

  implicit none

  include 'mpif.h' 

  !MPI stuff
  integer :: ierr,request, rank, size, root ! error,request, process,total_processes, root_process

  integer :: status(MPI_STATUS_SIZE) ! status size message

  !internal
  integer :: i, j, k   !counters
  integer :: n_proc, part_r, part_z ! number of processes, grid divisions r,z

  integer :: Nr,Nz,ghost, Nt, fout
  integer :: l_Nr,l_Nz
  real*8  :: dr,dz,dt,courant

  integer :: allocst

  character*20 :: filename 
  character*20  :: chrank


  real*8, allocatable, dimension(:,:) :: r,z
  real*8, allocatable, dimension(:,:) :: phi, phi_p, phi_s 
  real*8, allocatable, dimension(:,:) :: xir, xir_p, xir_s 
  real*8, allocatable, dimension(:,:) :: xiz, xiz_p, xiz_s 
  real*8, allocatable, dimension(:,:) :: pi, pi_p, pi_s 

  real*8 :: vr=1.D0,vz=0.D0,v=1.D0,mass=0.D0
  real*8 :: r0,z0, sig, A

  root = 0

! initialize MPI

  call MPI_INIT( ierr )
  call MPI_COMM_RANK( MPI_COMM_WORLD, rank, ierr )
  call MPI_COMM_SIZE( MPI_COMM_WORLD, size, ierr )


  if(rank==root)then
! Params, should be read

     read*,part_r
     read*,part_z

     read*,Nr
     read*,Nz

     read*,ghost

     read*,dr
     read*,dz

     read*,v


     read*, courant

     read*,Nt 
     read*,fout 
     
     read*,     r0 
     read*,     z0 
     
     read*,     sig
     
     read*,     A 



!********************************

     n_proc = part_r*part_z
  
     dt = courant*MIN(dr/v,dz/v)

! Make shure Nr is divisible by parts, at most allocates one more row per
! processor
 
     Nr= Nr+part_r-MOD(Nr,part_r)
     Nz= Nz+part_z-MOD(Nz,part_z)
     
     if(n_proc>size)then
        print*, 'Not enough processors assigned, need ', n_proc
        call MPI_FINALIZE( ierr )
        stop
     end if

     print*, 'Program AXIONDA start on ',size,'procesors'
     print*, ''
     print*, 'iteration,time'
     print*, '0         0'

  end if

! Broadcast relevant info
  call MPI_BCAST(part_r,1,MPI_INTEGER,root,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(part_z,1,MPI_INTEGER,root,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(n_proc,1,MPI_INTEGER,root,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(Nr,1,MPI_INTEGER,root,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(Nz,1,MPI_INTEGER,root,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(ghost,1,MPI_INTEGER,root,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(Nt,1,MPI_INTEGER,root,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(fout,1,MPI_INTEGER,root,MPI_COMM_WORLD,ierr)

  call MPI_BCAST(dr,2,MPI_REAL,root,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(dz,2,MPI_REAL,root,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(v,2,MPI_REAL,root,MPI_COMM_WORLD,ierr)
!  call MPI_BCAST(v,2,MPI_REAL,root,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(courant,2,MPI_REAL,root,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(dt,2,MPI_REAL,root,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(r0,2,MPI_REAL,root,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(z0,2,MPI_REAL,root,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(A,2,MPI_REAL,root,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(sig,2,MPI_REAL,root,MPI_COMM_WORLD,ierr)

  n_proc = part_r*part_z

  if(rank>=n_proc)then
     print*, 'Process inactive, shuting down ', rank
     call MPI_FINALIZE( ierr )
     stop
  end if

  if(n_proc>size)then
     call MPI_FINALIZE( ierr )
     stop
  end if

! allocate and init data

  l_Nr=Nr/part_r
  l_Nz=Nz/part_z

  allocate(r(1-ghost:l_Nr+ghost,1-ghost:l_Nz+ghost), z(1-ghost:l_Nr+ghost,1-ghost:l_Nz+ghost))
  allocate(phi(1-ghost:l_Nr+ghost,1-ghost:l_Nz+ghost),phi_p(1-ghost:l_Nr+ghost,1-ghost:l_Nz+ghost))
  allocate(phi_s(1-ghost:l_Nr+ghost,1-ghost:l_Nz+ghost))
  allocate(xir(1-ghost:l_Nr+ghost,1-ghost:l_Nz+ghost),xir_p(1-ghost:l_Nr+ghost,1-ghost:l_Nz+ghost))
  allocate(xir_s(1-ghost:l_Nr+ghost,1-ghost:l_Nz+ghost))
  allocate(xiz(1-ghost:l_Nr+ghost,1-ghost:l_Nz+ghost),xiz_p(1-ghost:l_Nr+ghost,1-ghost:l_Nz+ghost))
  allocate(xiz_s(1-ghost:l_Nr+ghost,1-ghost:l_Nz+ghost))
  allocate(pi(1-ghost:l_Nr+ghost,1-ghost:l_Nz+ghost),pi_p(1-ghost:l_Nr+ghost,1-ghost:l_Nz+ghost))
  allocate(pi_s(1-ghost:l_Nr+ghost,1-ghost:l_Nz+ghost))



  write(chrank,*)rank
  chrank = adjustl(chrank)
  filename='phi'//trim(chrank)//'.rl'
  open(unit=66+rank,file=filename,status='replace')


!  print*, rank,dr,dy

! Populate gridpoints

  do j = 1-ghost, l_Nz+ghost 
     do i = 1-ghost, l_Nr+ghost
        r(i,j)= MOD(rank,part_r)*l_Nr*dr + dr*i-1 +dr*.5D0
        z(i,j)= (rank/part_r)*l_Nz*dz + dz*j-1
     end do
  end do

  phi = A*exp( -  ( (r - r0)**2+(z- z0)**2)/sig**2)
  phi_p = 0.0
  phi_s = 0.0

  xir= -A*exp( -  ( (r - r0)**2+(z- z0)**2)/sig**2)*2.D0*(r-r0)/sig**2
  xir_p = 0.0
  xir_s = 0.0

  xiz= -A*exp( -  ( (r - r0)**2+(z- z0)**2)/sig**2)*2.D0*(z-z0)/sig**2
  xiz_p = 0.0
  xiz_s = 0.0

  pi=0.D0
  pi_p=0.D0
  pi_s=0.D0
  

! write initial data  
  do i = 1, l_Nz
     do j= 1, l_Nr
        write(66+rank,*) r(j,i),z(j,i), phi(j,i)
     end do
  end do
  write(66+rank,*) ''
  write(66+rank,*) ''

! Loop
  do i = 1, Nt
     phi_p= phi
     xir_p= xir
     xiz_p= xiz
     pi_p = pi
     ! Inner
! use ICN scheme
!
     !first step
     do k = 2-ghost, l_Nz+ghost-1
        do j = 2-ghost, l_Nr+ghost-1

!           phi_s(j,k) = -v*.5D0*(phi(j+1,k)-phi(j-1,k))/dr
           phi_s(j,k) = pi(j,k)
           xir_s(j,k) = .5D0*(pi(j+1,k)-pi(j-1,k))/dr
           xiz_s(j,k) = .5D0*(pi(j,k+1)-pi(j,k-1))/dz
           pi_s(j,k)  = v**2*.5D0* ( (xir(j+1,k)-xir(j-1,k))/dr &
                +(xiz(j,k+1)-xiz(j,k-1))/dz)+v**2*xir(j,k)/r(j,k)

           phi(j,k) = phi_p(j,k) + dt/2.0*phi_s(j,k)
           xir(j,k) = xir_p(j,k) + dt/2.0*xir_s(j,k)
           xiz(j,k) = xiz_p(j,k) + dt/2.0*xiz_s(j,k)
           pi(j,k)  = pi_p(j,k)  + dt/2.0*pi_s(j,k)

        end do
     end do

     !Second step
     do k = 3-ghost, l_Nz+ghost-2
        do j = 3-ghost, l_Nr+ghost-2
!           phi_s(j,k) = -v*.5D0*(phi(j+1,k)-phi(j-1,k))/dr
           phi_s(j,k) = pi(j,k)
           xir_s(j,k) = .5D0*(pi(j+1,k)-pi(j-1,k))/dr
           xiz_s(j,k) = .5D0*(pi(j,k+1)-pi(j,k-1))/dz
           pi_s(j,k)  = v**2*.5D0* ( (xir(j+1,k)-xir(j-1,k))/dr &
                +(xiz(j,k+1)-xiz(j,k-1))/dz)+v**2*xir(j,k)/r(j,k)

           phi(j,k) = phi_p(j,k) + dt/2.0*phi_s(j,k)
           xir(j,k) = xir_p(j,k) + dt/2.0*xir_s(j,k)
           xiz(j,k) = xiz_p(j,k) + dt/2.0*xiz_s(j,k)
           pi(j,k)  = pi_p(j,k)  + dt/2.0*pi_s(j,k)
        end do
     end do

     !Final step
     do k = 4-ghost, l_Nz+ghost-3
        do j = 4-ghost, l_Nr+ghost-3
!!$           phi_s(j,k) = -v*.5D0*(phi(j+1,k)-phi(j-1,k))/dr
           phi_s(j,k) = pi(j,k)
           xir_s(j,k) = .5D0*(pi(j+1,k)-pi(j-1,k))/dr
           xiz_s(j,k) = .5D0*(pi(j,k+1)-pi(j,k-1))/dz
           pi_s(j,k)  = v**2*.5D0* ( (xir(j+1,k)-xir(j-1,k))/dr &
                +(xiz(j,k+1)-xiz(j,k-1))/dz)+v**2*xir(j,k)/r(j,k)

           phi(j,k) = phi_p(j,k) + dt*phi_s(j,k)
           xir(j,k) = xir_p(j,k) + dt*xir_s(j,k)
           xiz(j,k) = xiz_p(j,k) + dt*xiz_s(j,k)
           pi(j,k)  = pi_p(j,k)  + dt*pi_s(j,k)
        end do
     end do


     ! Boundaries

     !Axis
     !we enforce regular behaviour
     if(Mod(rank,part_r)==0)then
        do j=1,ghost
           phi(1-j,:)=phi(j,:)
           xir(1-j,:)=-xir(j,:)
           xiz(1-j,:)=xiz(j,:)
           pi(1-j,:)=pi(j,:)
        end do
     end if

     !Domain Boundaries
     ! FIXED!!!
     if(Mod(rank,part_r)==partr-1)then
        phi(l_Nr,:)=phi_p(l_Nr,:)
        xir(l_Nr,:)=xir_p(l_Nr,:)
        xiz(l_Nr,:)=xiz_p(l_Nr,:)
         pi(l_Nr,:)= pi_p(l_Nr,:)
     end if
     if(rank/part_r==0)then
        phi(:,1)=phi_p(:,1)
        xir(:,1)=xir_p(:,1)
        xiz(:,1)=xiz_p(:,1)
         pi(:,1)= pi_p(:,1)
     end if
     if(rank/part_r==part_z-1)then
        phi(:,l_Nz)=phi_p(:,l_Nz)
        xir(:,l_Nz)=xir_p(:,l_Nz)
        xiz(:,l_Nz)=xiz_p(:,l_Nz)
         pi(:,l_Nz)= pi_p(:,l_Nz)
     end if
     !Sync ghostzones

     do k=1-ghost, l_Nz+ghost
        do j=1, ghost
           
           !ask for left ghostzone
           !left recievers
           if(Mod(rank,part_r)/=0)then
              call MPI_IRECV(phi(1-j,k),2,MPI_REAL,&
                   rank-1,1,MPI_COMM_WORLD,request,ierr) 
           end if
           !right senders
           if(Mod(rank,part_r)/=part_r-1)then
              call MPI_SEND(phi(l_Nr+1-j,k),2,MPI_REAL,&
                   rank+1,1,MPI_COMM_WORLD,ierr) 
           end if
           !Hang
           if(Mod(rank,part_r)/=0)then
              call MPI_WAIT(request,status,ierr) 
           end if
           
           !ask for right ghostzone
           !right recievers
           if(Mod(rank,part_r)/=part_r-1)then
              call MPI_IRECV(phi(l_Nr+j,k),2,MPI_REAL,&
                   rank+1,1,MPI_COMM_WORLD,request,ierr) 
           end if
           !left senders
           if(Mod(rank,part_r)/=0)then
              call MPI_SEND(phi(j,k),2,MPI_REAL,&
                   rank-1,1,MPI_COMM_WORLD,ierr) 
           end if
           !Hang
           if(Mod(rank,part_r)/=part_r-1)then
              call MPI_WAIT(request,status,ierr) 
           end if

           !XIR

           if(Mod(rank,part_r)/=0)then
              call MPI_IRECV(xir(1-j,k),2,MPI_REAL,&
                   rank-1,1,MPI_COMM_WORLD,request,ierr) 
           end if
           !right senders
           if(Mod(rank,part_r)/=part_r-1)then
              call MPI_SEND(xir(l_Nr+1-j,k),2,MPI_REAL,&
                   rank+1,1,MPI_COMM_WORLD,ierr) 
           end if
           !Hang
           if(Mod(rank,part_r)/=0)then
              call MPI_WAIT(request,status,ierr) 
           end if
           
           !ask for right ghostzone
           !right recievers
           if(Mod(rank,part_r)/=part_r-1)then
              call MPI_IRECV(xir(l_Nr+j,k),2,MPI_REAL,&
                   rank+1,1,MPI_COMM_WORLD,request,ierr) 
           end if
           !left senders
           if(Mod(rank,part_r)/=0)then
              call MPI_SEND(xir(j,k),2,MPI_REAL,&
                   rank-1,1,MPI_COMM_WORLD,ierr) 
           end if
           !Hang
           if(Mod(rank,part_r)/=part_r-1)then
              call MPI_WAIT(request,status,ierr) 
           end if

           !XIY
           if(Mod(rank,part_r)/=0)then
              call MPI_IRECV(xiz(1-j,k),2,MPI_REAL,&
                   rank-1,1,MPI_COMM_WORLD,request,ierr) 
           end if
           !right senders
           if(Mod(rank,part_r)/=part_r-1)then
              call MPI_SEND(xiz(l_Nr+1-j,k),2,MPI_REAL,&
                   rank+1,1,MPI_COMM_WORLD,ierr) 
           end if
           !Hang
           if(Mod(rank,part_r)/=0)then
              call MPI_WAIT(request,status,ierr) 
           end if
           
           !ask for right ghostzone
           !right recievers
           if(Mod(rank,part_r)/=part_r-1)then
              call MPI_IRECV(xiz(l_Nr+j,k),2,MPI_REAL,&
                   rank+1,1,MPI_COMM_WORLD,request,ierr) 
           end if
           !left senders
           if(Mod(rank,part_r)/=0)then
              call MPI_SEND(xiz(j,k),2,MPI_REAL,&
                   rank-1,1,MPI_COMM_WORLD,ierr) 
           end if
           !Hang
           if(Mod(rank,part_r)/=part_r-1)then
              call MPI_WAIT(request,status,ierr) 
           end if
           
           !Pi
           if(Mod(rank,part_r)/=0)then
              call MPI_IRECV(pi(1-j,k),2,MPI_REAL,&
                   rank-1,1,MPI_COMM_WORLD,request,ierr) 
           end if
           !right senders
           if(Mod(rank,part_r)/=part_r-1)then
              call MPI_SEND(pi(l_Nr+1-j,k),2,MPI_REAL,&
                   rank+1,1,MPI_COMM_WORLD,ierr) 
           end if
           !Hang
           if(Mod(rank,part_r)/=0)then
              call MPI_WAIT(request,status,ierr) 
           end if
           
           !ask for right ghostzone
           !right recievers
           if(Mod(rank,part_r)/=part_r-1)then
              call MPI_IRECV(pi(l_Nr+j,k),2,MPI_REAL,&
                   rank+1,1,MPI_COMM_WORLD,request,ierr) 
           end if
           !left senders
           if(Mod(rank,part_r)/=0)then
              call MPI_SEND(pi(j,k),2,MPI_REAL,&
                   rank-1,1,MPI_COMM_WORLD,ierr) 
           end if
           !Hang
           if(Mod(rank,part_r)/=part_r-1)then
              call MPI_WAIT(request,status,ierr) 
           end if

        end do
     end do
     
     do k=1, ghost
        do j=1-ghost, l_Nr+ghost
           !ask for lower ghostzone
           !low recievers
           if( rank/part_r/=0)then
              call MPI_IRECV(phi(j,1-k),2,MPI_REAL,&
                   rank-part_r,1,MPI_COMM_WORLD,request,ierr) 
           end if
           !up senders
           if( rank/part_r/=part_z-1)then
              call MPI_SEND(phi(j,l_Nz+1-k),2,MPI_REAL,&
                   rank+part_r,1,MPI_COMM_WORLD,ierr) 
           end if
           !Hang
           if( rank/part_r/=0)then
              call MPI_WAIT(request,status,ierr) 
           end if
           
           !ask for upper ghostzone
           !up recievers
           if( rank/part_r/=part_z-1)then
              call MPI_IRECV(phi(j,l_Nz+k),2,MPI_REAL,&
                   rank+part_r,1,MPI_COMM_WORLD,request,ierr) 
           end if
           !low senders
           if( rank/part_r/=0)then
              call MPI_SEND(phi(j,k),2,MPI_REAL,&
                   rank-part_r,1,MPI_COMM_WORLD,ierr) 
           end if
           !Hang
           if( rank/part_r/=part_z-1)then
              call MPI_WAIT(request,status,ierr) 
           end if

           !Rir
           !ask for lower ghostzone
           !low recievers
           if( rank/part_r/=0)then
              call MPI_IRECV(xir(j,1-k),2,MPI_REAL,&
                   rank-part_r,1,MPI_COMM_WORLD,request,ierr) 
           end if
           !up senders
           if( rank/part_r/=part_z-1)then
              call MPI_SEND(xir(j,l_Nz+1-k),2,MPI_REAL,&
                   rank+part_r,1,MPI_COMM_WORLD,ierr) 
           end if
           !Hang
           if( rank/part_r/=0)then
              call MPI_WAIT(request,status,ierr) 
           end if
           
           !ask for upper ghostzone
           !up recievers
           if( rank/part_r/=part_z-1)then
              call MPI_IRECV(xir(j,l_Nz+k),2,MPI_REAL,&
                   rank+part_r,1,MPI_COMM_WORLD,request,ierr) 
           end if
           !low senders
           if( rank/part_r/=0)then
              call MPI_SEND(xir(j,k),2,MPI_REAL,&
                   rank-part_r,1,MPI_COMM_WORLD,ierr) 
           end if
           !Hang
           if( rank/part_r/=part_z-1)then
              call MPI_WAIT(request,status,ierr) 
           end if
           ! Xiz
           !ask for lower ghostzone
           !low recievers
           if( rank/part_r/=0)then
              call MPI_IRECV(xiz(j,1-k),2,MPI_REAL,&
                   rank-part_r,1,MPI_COMM_WORLD,request,ierr) 
           end if
           !up senders
           if( rank/part_r/=part_z-1)then
              call MPI_SEND(xiz(j,l_Nz+1-k),2,MPI_REAL,&
                   rank+part_r,1,MPI_COMM_WORLD,ierr) 
           end if
           !Hang
           if( rank/part_r/=0)then
              call MPI_WAIT(request,status,ierr) 
           end if
           
           !ask for upper ghostzone
           !up recievers
           if( rank/part_r/=part_z-1)then
              call MPI_IRECV(xiz(j,l_Nz+k),2,MPI_REAL,&
                   rank+part_r,1,MPI_COMM_WORLD,request,ierr) 
           end if
           !low senders
           if( rank/part_r/=0)then
              call MPI_SEND(xiz(j,k),2,MPI_REAL,&
                   rank-part_r,1,MPI_COMM_WORLD,ierr) 
           end if
           !Hang
           if( rank/part_r/=part_z-1)then
              call MPI_WAIT(request,status,ierr) 
           end if
           ! Pi
           !ask for lower ghostzone
           !low recievers
           if( rank/part_r/=0)then
              call MPI_IRECV(pi(j,1-k),2,MPI_REAL,&
                   rank-part_r,1,MPI_COMM_WORLD,request,ierr) 
           end if
           !up senders
           if( rank/part_r/=part_z-1)then
              call MPI_SEND(pi(j,l_Nz+1-k),2,MPI_REAL,&
                   rank+part_r,1,MPI_COMM_WORLD,ierr) 
           end if
           !Hang
           if( rank/part_r/=0)then
              call MPI_WAIT(request,status,ierr) 
           end if
           
           !ask for upper ghostzone
           !up recievers
           if( rank/part_r/=part_z-1)then
              call MPI_IRECV(pi(j,l_Nz+k),2,MPI_REAL,&
                   rank+part_r,1,MPI_COMM_WORLD,request,ierr) 
           end if
           !low senders
           if( rank/part_r/=0)then
              call MPI_SEND(pi(j,k),2,MPI_REAL,&
                   rank-part_r,1,MPI_COMM_WORLD,ierr) 
           end if
           !Hang
           if( rank/part_r/=part_z-1)then
              call MPI_WAIT(request,status,ierr) 
           end if
        end do
     end do
     ! Output
     
     if(Mod(i,fout).eq.0) then
        do k = 1, l_Nz
           do j = 1, l_Nr
              write(66+rank,*) r(j,k),z(j,k), phi(j,k)
           end do
        end do
        write(66+rank,*) ''     
        write(66+rank,*) ''     
        if(rank==root)then
           print*, i, i*dt
        end if
     end if

  end do


  close(66+rank)

  deallocate( r,z)
  deallocate( phi,phi_p,phi_s)
  deallocate( xiz,xiz_p,xiz_s)
  deallocate( xir,xir_p,xir_s)
  deallocate( pi,pi_p,pi_s)


  if(rank==root)then
     print*, '----------------------'
     print*, 'Program ONDA finalized'
     print*, ''
  end if
  call MPI_FINALIZE( ierr )

end program axionda
