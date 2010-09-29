

program onda

! Solve wave equation in 2D, cartesian grid

  implicit none

  include 'mpif.h' 



  integer :: ierr,request, rank, size, root, i, j, k,status(MPI_STATUS_SIZE)

  integer :: n_proc, part_x, part_y

  integer :: Nx,Ny,ghost, Nt, fout
  integer :: l_Nx,l_Ny
  real*8  :: dx,dy,dt,courant

  integer :: allocst

  character*20 :: filename 
  character*20  :: chrank


  real*8, allocatable, dimension(:,:) :: x,y
  real*8, allocatable, dimension(:,:) :: phi, phi_p, phi_s 
  real*8, allocatable, dimension(:,:) :: xix, xix_p, xix_s 
  real*8, allocatable, dimension(:,:) :: xiy, xiy_p, xiy_s 
  real*8, allocatable, dimension(:,:) :: pi, pi_p, pi_s 

  real*8 :: vx=1.D0,vy=0.D0,v=1.D0,mass=0.D0
  real*8 :: x0,y0, sig, A

  root = 0

! initialize MPI

  call MPI_INIT( ierr )
  call MPI_COMM_RANK( MPI_COMM_WORLD, rank, ierr )
  call MPI_COMM_SIZE( MPI_COMM_WORLD, size, ierr )


  if(rank==root)then
! Params, should be read

     read*,part_x
     read*,part_y

     read*,Nx
     read*,Ny

     read*,ghost

     read*,dx
     read*,dy

     read*,v


     read*, courant

     read*,Nt 
     read*,fout 
     
     read*,     x0 
     read*,     y0 
     
     read*,     sig
     
     read*,     A 



!********************************

     n_proc = part_x*part_y
  
     dt = courant*MIN(dx/v,dy/v)

! Make shure Nx is divisible by parts, at most allocates one more row per
! processor
 
     Nx= Nx+part_x-MOD(Nx,part_x)
     Ny= Ny+part_y-MOD(Ny,part_y)
     
     if(n_proc>size)then
        print*, 'Not enough processors assigned, need ', n_proc
        call MPI_FINALIZE( ierr )
        stop
     end if

     print*, 'Program ONDA Start on ',size,'procesors'
     print*, ''
     print*, 'iteration,time'
     print*, '0         0'

  end if

! Broadcast relevant info
  call MPI_BCAST(part_x,1,MPI_INTEGER,root,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(part_y,1,MPI_INTEGER,root,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(n_proc,1,MPI_INTEGER,root,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(Nx,1,MPI_INTEGER,root,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(Ny,1,MPI_INTEGER,root,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(ghost,1,MPI_INTEGER,root,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(Nt,1,MPI_INTEGER,root,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(fout,1,MPI_INTEGER,root,MPI_COMM_WORLD,ierr)

  call MPI_BCAST(dx,2,MPI_REAL,root,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(dy,2,MPI_REAL,root,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(v,2,MPI_REAL,root,MPI_COMM_WORLD,ierr)
!  call MPI_BCAST(v,2,MPI_REAL,root,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(courant,2,MPI_REAL,root,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(dt,2,MPI_REAL,root,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(x0,2,MPI_REAL,root,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(y0,2,MPI_REAL,root,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(A,2,MPI_REAL,root,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(sig,2,MPI_REAL,root,MPI_COMM_WORLD,ierr)

  n_proc = part_x*part_y

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

  l_Nx=Nx/part_x
  l_Ny=Ny/part_y

  allocate(x(1-ghost:l_Nx+ghost,1-ghost:l_Ny+ghost), y(1-ghost:l_Nx+ghost,1-ghost:l_Ny+ghost))
  allocate(phi(1-ghost:l_Nx+ghost,1-ghost:l_Ny+ghost),phi_p(1-ghost:l_Nx+ghost,1-ghost:l_Ny+ghost))
  allocate(phi_s(1-ghost:l_Nx+ghost,1-ghost:l_Ny+ghost))
  allocate(xix(1-ghost:l_Nx+ghost,1-ghost:l_Ny+ghost),xix_p(1-ghost:l_Nx+ghost,1-ghost:l_Ny+ghost))
  allocate(xix_s(1-ghost:l_Nx+ghost,1-ghost:l_Ny+ghost))
  allocate(xiy(1-ghost:l_Nx+ghost,1-ghost:l_Ny+ghost),xiy_p(1-ghost:l_Nx+ghost,1-ghost:l_Ny+ghost))
  allocate(xiy_s(1-ghost:l_Nx+ghost,1-ghost:l_Ny+ghost))
  allocate(pi(1-ghost:l_Nx+ghost,1-ghost:l_Ny+ghost),pi_p(1-ghost:l_Nx+ghost,1-ghost:l_Ny+ghost))
  allocate(pi_s(1-ghost:l_Nx+ghost,1-ghost:l_Ny+ghost))



  write(chrank,*)rank
  chrank = adjustl(chrank)
  filename='phi'//trim(chrank)//'.rl'
  open(unit=66+rank,file=filename,status='replace')


!  print*, rank,dx,dy

! Populate gridpoints

  do j = 1-ghost, l_Ny+ghost 
     do i = 1-ghost, l_Nx+ghost
        x(i,j)= MOD(rank,part_x)*l_Nx*dx + dx*i-1
        y(i,j)= (rank/part_x)*l_Ny*dy + dy*j-1
     end do
  end do

  phi = A*exp( -  ( (x - x0)**2+(y- y0)**2)/sig**2)
  phi_p = 0.0
  phi_s = 0.0

  xix= -A*exp( -  ( (x - x0)**2+(y- y0)**2)/sig**2)*2.D0*(x-x0)/sig**2
  xix_p = 0.0
  xix_s = 0.0

  xiy= -A*exp( -  ( (x - x0)**2+(y- y0)**2)/sig**2)*2.D0*(y-y0)/sig**2
  xiy_p = 0.0
  xiy_s = 0.0

  pi=0.D0
  pi_p=0.D0
  pi_s=0.D0
  

! write initial data  
  do i = 1, l_Ny
     do j= 1, l_Nx
        write(66+rank,*) x(j,i),y(j,i), phi(j,i)
     end do
  end do
  write(66+rank,*) ''
  write(66+rank,*) ''

! Loop
  do i = 1, Nt
     phi_p= phi
     xix_p= xix
     xiy_p= xiy
     pi_p = pi
     ! Inner
! use ICN scheme
!
     !first step
     do k = 2-ghost, l_Ny+ghost-1
        do j = 2-ghost, l_Nx+ghost-1

!           phi_s(j,k) = -v*.5D0*(phi(j+1,k)-phi(j-1,k))/dx
           phi_s(j,k) = pi(j,k)
           xix_s(j,k) = .5D0*(pi(j+1,k)-pi(j-1,k))/dx
           xiy_s(j,k) = .5D0*(pi(j,k+1)-pi(j,k-1))/dy
           pi_s(j,k)  = v**2*.5D0* ( (xix(j+1,k)-xix(j-1,k))/dx &
                +(xiy(j,k+1)-xiy(j,k-1))/dy)

           phi(j,k) = phi_p(j,k) + dt/2.0*phi_s(j,k)
           xix(j,k) = xix_p(j,k) + dt/2.0*xix_s(j,k)
           xiy(j,k) = xiy_p(j,k) + dt/2.0*xiy_s(j,k)
           pi(j,k)  = pi_p(j,k) + dt/2.0*pi_s(j,k)

        end do
     end do

     !Second step
     do k = 3-ghost, l_Ny+ghost-2
        do j = 3-ghost, l_Nx+ghost-2
!           phi_s(j,k) = -v*.5D0*(phi(j+1,k)-phi(j-1,k))/dx
           phi_s(j,k) = pi(j,k)
           xix_s(j,k) = .5D0*(pi(j+1,k)-pi(j-1,k))/dx
           xiy_s(j,k) = .5D0*(pi(j,k+1)-pi(j,k-1))/dy
           pi_s(j,k)  = v**2*.5D0* ( (xix(j+1,k)-xix(j-1,k))/dx &
                +(xiy(j,k+1)-xiy(j,k-1))/dy)

           phi(j,k) = phi_p(j,k) + dt/2.0*phi_s(j,k)
           xix(j,k) = xix_p(j,k) + dt/2.0*xix_s(j,k)
           xiy(j,k) = xiy_p(j,k) + dt/2.0*xiy_s(j,k)
           pi(j,k)  = pi_p(j,k) + dt/2.0*pi_s(j,k)
        end do
     end do

     !Final step
     do k = 4-ghost, l_Ny+ghost-3
        do j = 4-ghost, l_Nx+ghost-3
!!$           phi_s(j,k) = -v*.5D0*(phi(j+1,k)-phi(j-1,k))/dx
           phi_s(j,k) = pi(j,k)
           xix_s(j,k) = .5D0*(pi(j+1,k)-pi(j-1,k))/dx
           xiy_s(j,k) = .5D0*(pi(j,k+1)-pi(j,k-1))/dy
           pi_s(j,k)  = v**2*.5D0* ( (xix(j+1,k)-xix(j-1,k))/dx &
                +(xiy(j,k+1)-xiy(j,k-1))/dy)

           phi(j,k) = phi_p(j,k) + dt*phi_s(j,k)
           xix(j,k) = xix_p(j,k) + dt*xix_s(j,k)
           xiy(j,k) = xiy_p(j,k) + dt*xiy_s(j,k)
           pi(j,k)  = pi_p(j,k) + dt*pi_s(j,k)
        end do
     end do


     ! Nothing on boundaries

     !Sync ghostzones

     do k=1-ghost, l_Ny+ghost
        do j=1, ghost
           
           !ask for left ghostzone
           !left recievers
           if(Mod(rank,part_x)/=0)then
              call MPI_IRECV(phi(1-j,k),2,MPI_REAL,&
                   rank-1,1,MPI_COMM_WORLD,request,ierr) 
           end if
           !right senders
           if(Mod(rank,part_x)/=part_x-1)then
              call MPI_SEND(phi(l_Nx+1-j,k),2,MPI_REAL,&
                   rank+1,1,MPI_COMM_WORLD,ierr) 
           end if
           !Hang
           if(Mod(rank,part_x)/=0)then
              call MPI_WAIT(request,status,ierr) 
           end if
           
           !ask for right ghostzone
           !right recievers
           if(Mod(rank,part_x)/=part_x-1)then
              call MPI_IRECV(phi(l_Nx+j,k),2,MPI_REAL,&
                   rank+1,1,MPI_COMM_WORLD,request,ierr) 
           end if
           !left senders
           if(Mod(rank,part_x)/=0)then
              call MPI_SEND(phi(j,k),2,MPI_REAL,&
                   rank-1,1,MPI_COMM_WORLD,ierr) 
           end if
           !Hang
           if(Mod(rank,part_x)/=part_x-1)then
              call MPI_WAIT(request,status,ierr) 
           end if

           !XIX

           if(Mod(rank,part_x)/=0)then
              call MPI_IRECV(xix(1-j,k),2,MPI_REAL,&
                   rank-1,1,MPI_COMM_WORLD,request,ierr) 
           end if
           !right senders
           if(Mod(rank,part_x)/=part_x-1)then
              call MPI_SEND(xix(l_Nx+1-j,k),2,MPI_REAL,&
                   rank+1,1,MPI_COMM_WORLD,ierr) 
           end if
           !Hang
           if(Mod(rank,part_x)/=0)then
              call MPI_WAIT(request,status,ierr) 
           end if
           
           !ask for right ghostzone
           !right recievers
           if(Mod(rank,part_x)/=part_x-1)then
              call MPI_IRECV(xix(l_Nx+j,k),2,MPI_REAL,&
                   rank+1,1,MPI_COMM_WORLD,request,ierr) 
           end if
           !left senders
           if(Mod(rank,part_x)/=0)then
              call MPI_SEND(xix(j,k),2,MPI_REAL,&
                   rank-1,1,MPI_COMM_WORLD,ierr) 
           end if
           !Hang
           if(Mod(rank,part_x)/=part_x-1)then
              call MPI_WAIT(request,status,ierr) 
           end if

           !XIY
           if(Mod(rank,part_x)/=0)then
              call MPI_IRECV(xiy(1-j,k),2,MPI_REAL,&
                   rank-1,1,MPI_COMM_WORLD,request,ierr) 
           end if
           !right senders
           if(Mod(rank,part_x)/=part_x-1)then
              call MPI_SEND(xiy(l_Nx+1-j,k),2,MPI_REAL,&
                   rank+1,1,MPI_COMM_WORLD,ierr) 
           end if
           !Hang
           if(Mod(rank,part_x)/=0)then
              call MPI_WAIT(request,status,ierr) 
           end if
           
           !ask for right ghostzone
           !right recievers
           if(Mod(rank,part_x)/=part_x-1)then
              call MPI_IRECV(xiy(l_Nx+j,k),2,MPI_REAL,&
                   rank+1,1,MPI_COMM_WORLD,request,ierr) 
           end if
           !left senders
           if(Mod(rank,part_x)/=0)then
              call MPI_SEND(xiy(j,k),2,MPI_REAL,&
                   rank-1,1,MPI_COMM_WORLD,ierr) 
           end if
           !Hang
           if(Mod(rank,part_x)/=part_x-1)then
              call MPI_WAIT(request,status,ierr) 
           end if
           
           !Pi
           if(Mod(rank,part_x)/=0)then
              call MPI_IRECV(pi(1-j,k),2,MPI_REAL,&
                   rank-1,1,MPI_COMM_WORLD,request,ierr) 
           end if
           !right senders
           if(Mod(rank,part_x)/=part_x-1)then
              call MPI_SEND(pi(l_Nx+1-j,k),2,MPI_REAL,&
                   rank+1,1,MPI_COMM_WORLD,ierr) 
           end if
           !Hang
           if(Mod(rank,part_x)/=0)then
              call MPI_WAIT(request,status,ierr) 
           end if
           
           !ask for right ghostzone
           !right recievers
           if(Mod(rank,part_x)/=part_x-1)then
              call MPI_IRECV(pi(l_Nx+j,k),2,MPI_REAL,&
                   rank+1,1,MPI_COMM_WORLD,request,ierr) 
           end if
           !left senders
           if(Mod(rank,part_x)/=0)then
              call MPI_SEND(pi(j,k),2,MPI_REAL,&
                   rank-1,1,MPI_COMM_WORLD,ierr) 
           end if
           !Hang
           if(Mod(rank,part_x)/=part_x-1)then
              call MPI_WAIT(request,status,ierr) 
           end if

        end do
     end do
     
     do k=1, ghost
        do j=1-ghost, l_Nx+ghost
           !ask for lower ghostzone
           !low recievers
           if( rank/part_x/=0)then
              call MPI_IRECV(phi(j,1-k),2,MPI_REAL,&
                   rank-part_x,1,MPI_COMM_WORLD,request,ierr) 
           end if
           !up senders
           if( rank/part_x/=part_y-1)then
              call MPI_SEND(phi(j,l_Ny+1-k),2,MPI_REAL,&
                   rank+part_x,1,MPI_COMM_WORLD,ierr) 
           end if
           !Hang
           if( rank/part_x/=0)then
              call MPI_WAIT(request,status,ierr) 
           end if
           
           !ask for upper ghostzone
           !up recievers
           if( rank/part_x/=part_y-1)then
              call MPI_IRECV(phi(j,l_Ny+k),2,MPI_REAL,&
                   rank+part_x,1,MPI_COMM_WORLD,request,ierr) 
           end if
           !low senders
           if( rank/part_x/=0)then
              call MPI_SEND(phi(j,k),2,MPI_REAL,&
                   rank-part_x,1,MPI_COMM_WORLD,ierr) 
           end if
           !Hang
           if( rank/part_x/=part_y-1)then
              call MPI_WAIT(request,status,ierr) 
           end if

           !Xix
           !ask for lower ghostzone
           !low recievers
           if( rank/part_x/=0)then
              call MPI_IRECV(xix(j,1-k),2,MPI_REAL,&
                   rank-part_x,1,MPI_COMM_WORLD,request,ierr) 
           end if
           !up senders
           if( rank/part_x/=part_y-1)then
              call MPI_SEND(xix(j,l_Ny+1-k),2,MPI_REAL,&
                   rank+part_x,1,MPI_COMM_WORLD,ierr) 
           end if
           !Hang
           if( rank/part_x/=0)then
              call MPI_WAIT(request,status,ierr) 
           end if
           
           !ask for upper ghostzone
           !up recievers
           if( rank/part_x/=part_y-1)then
              call MPI_IRECV(xix(j,l_Ny+k),2,MPI_REAL,&
                   rank+part_x,1,MPI_COMM_WORLD,request,ierr) 
           end if
           !low senders
           if( rank/part_x/=0)then
              call MPI_SEND(xix(j,k),2,MPI_REAL,&
                   rank-part_x,1,MPI_COMM_WORLD,ierr) 
           end if
           !Hang
           if( rank/part_x/=part_y-1)then
              call MPI_WAIT(request,status,ierr) 
           end if
           ! Xiy
           !ask for lower ghostzone
           !low recievers
           if( rank/part_x/=0)then
              call MPI_IRECV(xiy(j,1-k),2,MPI_REAL,&
                   rank-part_x,1,MPI_COMM_WORLD,request,ierr) 
           end if
           !up senders
           if( rank/part_x/=part_y-1)then
              call MPI_SEND(xiy(j,l_Ny+1-k),2,MPI_REAL,&
                   rank+part_x,1,MPI_COMM_WORLD,ierr) 
           end if
           !Hang
           if( rank/part_x/=0)then
              call MPI_WAIT(request,status,ierr) 
           end if
           
           !ask for upper ghostzone
           !up recievers
           if( rank/part_x/=part_y-1)then
              call MPI_IRECV(xiy(j,l_Ny+k),2,MPI_REAL,&
                   rank+part_x,1,MPI_COMM_WORLD,request,ierr) 
           end if
           !low senders
           if( rank/part_x/=0)then
              call MPI_SEND(xiy(j,k),2,MPI_REAL,&
                   rank-part_x,1,MPI_COMM_WORLD,ierr) 
           end if
           !Hang
           if( rank/part_x/=part_y-1)then
              call MPI_WAIT(request,status,ierr) 
           end if
           ! Pi
           !ask for lower ghostzone
           !low recievers
           if( rank/part_x/=0)then
              call MPI_IRECV(pi(j,1-k),2,MPI_REAL,&
                   rank-part_x,1,MPI_COMM_WORLD,request,ierr) 
           end if
           !up senders
           if( rank/part_x/=part_y-1)then
              call MPI_SEND(pi(j,l_Ny+1-k),2,MPI_REAL,&
                   rank+part_x,1,MPI_COMM_WORLD,ierr) 
           end if
           !Hang
           if( rank/part_x/=0)then
              call MPI_WAIT(request,status,ierr) 
           end if
           
           !ask for upper ghostzone
           !up recievers
           if( rank/part_x/=part_y-1)then
              call MPI_IRECV(pi(j,l_Ny+k),2,MPI_REAL,&
                   rank+part_x,1,MPI_COMM_WORLD,request,ierr) 
           end if
           !low senders
           if( rank/part_x/=0)then
              call MPI_SEND(pi(j,k),2,MPI_REAL,&
                   rank-part_x,1,MPI_COMM_WORLD,ierr) 
           end if
           !Hang
           if( rank/part_x/=part_y-1)then
              call MPI_WAIT(request,status,ierr) 
           end if
        end do
     end do
     ! Output
     
     if(Mod(i,fout).eq.0) then
        do k = 1, l_Ny
           do j = 1, l_Nx
              write(66+rank,*) x(j,k),y(j,k), phi(j,k)
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

  deallocate( x,y)
  deallocate( phi,phi_p,phi_s)
  deallocate( xiy,xiy_p,xiy_s)
  deallocate( xix,xix_p,xix_s)
  deallocate( pi,pi_p,pi_s)


  if(rank==root)then
     print*, '----------------------'
     print*, 'Program ONDA finalized'
     print*, ''
  end if
  call MPI_FINALIZE( ierr )

end program onda
