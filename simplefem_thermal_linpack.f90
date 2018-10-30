!==============================================!
!                                              !
! C. thieulot ; March 2016                     !
!                                              !
!==============================================!
                                               !
program fcubed                                 !
                                               !
implicit none                                  !
                                               !
integer, parameter :: m=4                      ! number of nodes which constitute an element
integer, parameter :: ndofV=2                  ! number of velocity dofs per node
integer, parameter :: ndofT=1                  ! number of temperature dofs per node
integer nnx                                    ! number of grid points in the x direction
integer nny                                    ! number of grid points in the y direction
integer np                                     ! number of grid points
integer nelx                                   ! number of elements in the x direction
integer nely                                   ! number of elements in the y direction
integer nel                                    ! number of elements
integer NfemV                                  ! size of the FEM Stokes matrix 
integer NfemT                                  ! size of the FEM temperature matrix 
integer nstep                                  ! number of timesteps 
integer output_freq                            ! how often vtu output is triggered
integer, dimension(:,:), allocatable :: icon   ! connectivity array
integer, dimension(:,:), allocatable :: icon_inv! inverse connectivity array

integer, dimension(:), allocatable :: ipvt     ! work array needed by the solver 
                                               !
integer i1,i2,i,j,k,iel,counter,iq,jq,istep    ! integer variables needed 
integer ik,jk,ikk,jkk,m1,m2,k1,k2,job          ! by the code
                                               !  
real(8) Lx,Ly                                  ! size of the numerical domain
real(8) viscosity                              ! dynamic viscosity $\mu$ of the material
real(8) rho0                                   ! reference mass density $\rho_0$ of the material
real(8) gx,gy                                  ! gravity acceleration
real(8), dimension(:),   allocatable :: x,y    ! node coordinates arrays
real(8), dimension(:),   allocatable :: u,v    ! node velocity arrays
real(8), dimension(:),   allocatable :: press_el !pressure , elemental
real(8), dimension(:),   allocatable :: press_node !also pressure, but nodal 
real(8), dimension(:),   allocatable :: BV     ! right hand side of Stokes system
real(8), dimension(:,:), allocatable :: AV     ! FEM matrix of Stokes system
real(8), dimension(:),   allocatable :: B_T    ! right hand side of temperature system
real(8), dimension(:,:), allocatable :: A_T    ! FEM matrix of temperature system
real(8), dimension(:),   allocatable :: work   ! work array needed by the solver
real(8), dimension(:),   allocatable :: bc_valV! array containing bc values for velocity
real(8), dimension(:),   allocatable :: bc_valT! array containing bc values for temperature
real(8), dimension(:),   allocatable :: bc_valP
real(8), dimension(:),   allocatable :: T      ! node temperature array
real(8), dimension(:),   allocatable :: density! node density array
real(8), dimension(:),   allocatable :: qx,qy  ! heat flux vector 
real(8), dimension(:),   allocatable :: Tavrg  !
                                               !
real(8) rq,sq,wq                               ! local coordinate and weight of qpoint
real(8) xq,yq                                  ! global coordinate of qpoint
real(8) uq,vq                                  ! velocity at qpoint
real(8) Tq                                     ! temperature at qpoint
real(8) exxq,eyyq,exyq                         ! strain-rate components at qpoint  
real(8) AelV(m*ndofV,m*ndofV)                  ! elemental Stokes FEM matrix
real(8) BelV(m*ndofV)                          ! elemental Stokes right hand side
real(8) AelT(m*ndofT,m*ndofT)                  ! elemental temperature FEM matrix
real(8) BelT(m*ndofT)                          ! elemental temperature right hand side
real(8) N(m),dNdx(m),dNdy(m),dNdr(m),dNds(m)   ! shape fcts and derivatives
real(8) jcob                                   ! determinant of jacobian matrix
real(8) jcb(2,2)                               ! jacobian matrix
real(8) jcbi(2,2)                              ! inverse of jacobian matrix
real(8) BmatV(3,ndofV*m)                       ! B matrix
real(8), dimension(3,3) :: Kmat                ! K matrix 
real(8), dimension(3,3) :: Cmat                ! C matrix
real(8) Aref                                   !
real(8) hcapa                                  ! heat capacity
real(8) hcond                                  ! heat conductivity
real(8) alpha                                  ! thermal expansion 
real(8) time                                   ! real time
real(8) dt                                     ! timestep
real(8) courant_nb                             ! Courant number for CFL criterion
real(8) Ra                                     ! Rayleigh number 
real(8) chi_T
real(8), parameter :: theta = 0.5              ! mid-point timestepping parameter
real(8), parameter :: eps=1.d-10               !
real(8), parameter :: pi = 3.14159265359d0     !
real(8), parameter :: tol=1.d-9                ! convergence tolerance
real(8) M_T(4,4),Ka(4,4),Kc(4,4),KK(4,4),Kb(4,4) ! various FEM arrays
real(8) Nvect(1,4),NvectT(4,1),NvectTstar(4,1) !
real(8) velnorm,vel2D(1,2),temp(4)             !
real(8) BmatT(2,4),BmatTT(4,2)                 !
real(8) hx,hy                                  ! grid spacing
real(8) umax,vmax,rho                          !
real(8) dTdxq,dTdyq                            ! temperature gradient at quad point
real(8) vrms,vrms_prev,Nu,Nu_prev              ! v_rms and Nu number
real(8) chi_Nu,chi_vrms                        ! convergence indicator
real(8) rcond                                  ! parameter needed for solver
                                               !
logical, dimension(:), allocatable :: bc_fixV  ! prescribed b.c. array for velocities
logical, dimension(:), allocatable :: bc_fixT  ! prescribed b.c. array for temperature
logical, dimension(:), allocatable :: bc_fixP
                                               !

real(8) Kel(m*ndofV,m*ndofV)                   ! elemental FEM matrix
real(8) fel(m*ndofV)                           ! elemental right hand side Vel
real(8) hel(1)                                 ! elemental right hand side Press
real(8) Gel(m*ndofV,1)                         !
real(8) rho2,drhody2,drhodx2    
real(8), dimension(:,:), allocatable :: KKKmem ! FEM matrix
real(8), dimension(:,:), allocatable :: KKK    ! FEM matrix

real(8) Bmat(3,ndofV*m)                        ! B matrix
real(8), dimension(:),   allocatable :: rhs_f  ! right hand side
real(8), dimension(:),   allocatable :: rhs_h  ! right hand side
real(8), dimension(:,:), allocatable :: G      ! gradient operator matrix

integer ipicard,npicard,solver
real(8) P_diff,V_diff
real(8), dimension(:),   allocatable :: Vsol
real(8), dimension(:),   allocatable :: Vsolmem
real(8), dimension(:),   allocatable :: Psolmem  
real(8), dimension(:),   allocatable :: Rv
real(8), dimension(:),   allocatable :: Rp
real(8) mass,ke,gpe,ke_density,gpe_density
real(8) thermal_energy, thermal_energy_density
real(8) ee, ee_density
real(8) phi_q_el,phi_q_top,phi_q_bottom
real(8) element_boundary_value
logical node_is_boundary
real(8) pe, pe_density,pe2

real(8) work_energy, work_el, normal_y, normal_x

real(8) T0, Di
real(8) mean_surface_pressure
integer number_of_surface_nodes

real(8),  parameter :: D = 3.d0
real(8),  parameter :: Dinv1 = -1.d0/D
real(8),  parameter :: Dinv2 = 1 + Dinv1

real(8) phi,rho_el,v_el

logical compressible, ALA, TALA, BA, EBA,buoyancy
integer approximation

approximation = 1
!Chose which approx to use
compressible=.false.
ALA = .false.
TALA = .false.
BA = .false.
EBA = .false.
buoyancy = .false.
select case(approximation)
case(1)
  !ALA
  ALA = .true.
  compressible = .true.
  buoyancy = .true.
case(2)
  !TALA
  TALA = .true.
  compressible = .true.
  buoyancy = .true.
case(3)
  !BA 
  BA = .true.
  buoyancy = .false.
case(4)
  !EBA 
  EBA = .true.
  buoyancy = .true.
end select



!==============================================!
!=====[setup]==================================!
!==============================================!

! Lx=1.d0
! Ly=1.d0

! nnx=21
! nny=21

! nelx=nnx-1
! nely=nny-1

! viscosity=1.d-3
! rho0=1.d0
! hcapa=1.d0
! hcond=1.d0
! alpha=1.d-2
 chi_T=1.d-5
 T0=1.d0

! gx=0.d0
! gy=-10

Lx=1.d0
Ly=1.d0

nnx=21
nny=21

nelx=nnx-1
nely=nny-1

viscosity=1.d-2
rho0=1.d0
hcapa=1.d0 !Massive performance loss associated with parameter
hcond=1.d0
alpha=1.d-2

Ra=1.d2

gx=0
gy=-Ra/alpha

courant_nb=0.5

nstep=1000

output_freq=10

npicard=1000

solver=2

Ra = -alpha*gy*T0*lx**3 / (hcond*viscosity)
Di = -alpha*gy*Lx / hcapa
write(*,*) "Ra = ", Ra, "Di = ", Di
!==============================================!
!==============================================!

time=0.d0

nel=nelx*nely ! total number of elements

np=nnx*nny ! total number of nodes

NfemV=np*ndofV ! size of Stokes matrix
NfemT=np*ndofT ! size of temperature matrix

hx=Lx/(nnx-1) ! grid spacing in x direction
hy=Ly/(nny-1) ! grid spacing in y direction

!==============================================!

Kmat(1,1)=1.d0 ; Kmat(1,2)=1.d0 ; Kmat(1,3)=0.d0  
Kmat(2,1)=1.d0 ; Kmat(2,2)=1.d0 ; Kmat(2,3)=0.d0  
Kmat(3,1)=0.d0 ; Kmat(3,2)=0.d0 ; Kmat(3,3)=0.d0  

Cmat(1,1)=2.d0 ; Cmat(1,2)=0.d0 ; Cmat(1,3)=0.d0  
Cmat(2,1)=0.d0 ; Cmat(2,2)=2.d0 ; Cmat(2,3)=0.d0  
Cmat(3,1)=0.d0 ; Cmat(3,2)=0.d0 ; Cmat(3,3)=1.d0

!==============================================!

print *,'Lx=',Lx
print *,'Ly=',Ly
print *,'hx=',hx
print *,'hy=',hy
print *,'nelx',nelx
print *,'nely',nely
print *,'nel',nel
print *,'nnx',nnx
print *,'nny',nny
print *,'np',np
print *,'gx=',gx
print *,'gy=',gy
print *,'NfemV=',NfemV
print *,'NfemT=',NfemT
print *,'courant_nb=',courant_nb
print *,'Ra=',Ra

!==============================================!
!===[open files]===============================!
!==============================================!

open(unit=1000,file='OUT/velocity_stats.dat')
open(unit=1001,file='OUT/pressure_stats.dat')
open(unit=1002,file='OUT/temperature_stats.dat')
open(unit=1003,file='OUT/vrms.dat')
open(unit=1004,file='OUT/Nu.dat')
open(unit=1005,file='OUT/dt.dat')
open(unit=1006,file='OUT/conv_Nu.dat')
open(unit=1007,file='OUT/conv_vrms.dat')
open(unit=1008,file='OUT/energies.dat',status='replace')

!==============================================!
!===[allocate memory]==========================!
!==============================================!

allocate(x(np))
allocate(y(np))
allocate(u(np))
allocate(v(np))
allocate(icon(m,nel))
allocate(icon_inv(m,np))
allocate(AV(NfemV,NfemV))
allocate(BV(NfemV))
allocate(A_T(NfemT,NfemT))
allocate(B_T(NfemT))
allocate(bc_fixV(NfemV))
allocate(bc_fixT(NfemT))
allocate(bc_fixP(Nel))
allocate(bc_valV(NfemV))
allocate(bc_valT(NfemT))
allocate(bc_valP(nel))
allocate(press_el(nel))
allocate(T(np))
allocate(density(np))
allocate(qx(nel),qy(nel))
allocate(Tavrg(nny))

allocate(KKK(NfemV,NfemV))
allocate(KKKmem(NfemV,NfemV))
allocate(rhs_f(NfemV))
allocate(rhs_h(nel))
allocate(G(NfemV,nel))
allocate(Vsol(NfemV))
allocate(Vsolmem(NfemV))
allocate(press_node(np))
allocate(Psolmem(Nel))
allocate(Rv(NfemV))
allocate(Rp(nel))


!==============================================!
!===[grid points setup]========================!
!==============================================!

counter=0
do j=0,nely
   do i=0,nelx
      counter=counter+1
      x(counter)=dble(i)*hx
      y(counter)=dble(j)*hy
   end do
end do

!open(unit=123,file='OUT/gridnodes.dat',status='replace')
!write(123,'(a)') '#     xpos      ypos    node '
!do i=1,np
!   write(123,'(2f10.5,i8)') x(i),y(i),i
!end do
!close(123)

!==============================================!
!===[connectivity]=============================!
!==============================================!

counter=0
do j=1,nely
   do i=1,nelx
      counter=counter+1
      icon(1,counter)=i+(j-1)*(nelx+1)
      icon(2,counter)=i+1+(j-1)*(nelx+1)
      icon(3,counter)=i+1+j*(nelx+1)
      icon(4,counter)=i+j*(nelx+1)
   end do
end do

call inverse_icon(icon,icon_inv,nel,np,m)

!open(unit=123,file='OUT/icon.dat',status='replace')
!do iel=1,nel
!   write(123,'(a)') '----------------------------'
!   write(123,'(a,i4,a)') '---element #',iel,' -----------'
!   write(123,'(a)') '----------------------------'
!   write(123,'(a,i8,a,2f20.10)') '  node 1 ', icon(1,iel),' at pos. ',x(icon(1,iel)),y(icon(1,iel))
!   write(123,'(a,i8,a,2f20.10)') '  node 2 ', icon(2,iel),' at pos. ',x(icon(2,iel)),y(icon(2,iel))
!   write(123,'(a,i8,a,2f20.10)') '  node 3 ', icon(3,iel),' at pos. ',x(icon(3,iel)),y(icon(3,iel))
!   write(123,'(a,i8,a,2f20.10)') '  node 4 ', icon(4,iel),' at pos. ',x(icon(4,iel)),y(icon(4,iel))
!end do
!close(123)

!==============================================!
!=====[define bc for velocity]=================!
!==============================================!

bc_fixV=.false.

do i=1,np
   if (x(i).lt.eps) then
      bc_fixV((i-1)*ndofV+1)=.true. ; bc_valV((i-1)*ndofV+1)=0.d0
      !bc_fixV((i-1)*ndofV+2)=.true. ; bc_valV((i-1)*ndofV+2)=0.d0
   endif
   if (x(i).gt.(Lx-eps)) then
      bc_fixV((i-1)*ndofV+1)=.true. ; bc_valV((i-1)*ndofV+1)=0.d0
      !bc_fixV((i-1)*ndofV+2)=.true. ; bc_valV((i-1)*ndofV+2)=0.d0
   endif
   if (y(i).lt.eps) then
      !bc_fixV((i-1)*ndofV+1)=.true. ; bc_valV((i-1)*ndofV+1)=0.d0
      bc_fixV((i-1)*ndofV+2)=.true. ; bc_valV((i-1)*ndofV+2)=0.d0
   endif
   if (y(i).gt.(Ly-eps) ) then
      !bc_fixV((i-1)*ndofV+1)=.true. ; bc_valV((i-1)*ndofV+1)=0.d0
      bc_fixV((i-1)*ndofV+2)=.true. ; bc_valV((i-1)*ndofV+2)=0.d0
   endif
end do

!open(unit=123,file='OUT/bc_u.dat',status='replace')
!open(unit=234,file='OUT/bc_v.dat',status='replace')
!do i=1,np
!   if (bc_fixV((i-1)*ndofV+1)) write(123,'(3f20.10)') x(i),y(i),bc_valV((i-1)*ndofV+1) 
!   if (bc_fixV((i-1)*ndofV+2)) write(234,'(3f20.10)') x(i),y(i),bc_valV((i-1)*ndofV+2) 
!end do
!close(123)
!close(234)

!==============================================!
!=====[define bc for temperature]==============!
!==============================================!

bc_fixT=.false.

do i=1,np
   if (y(i).lt.eps) then 
      bc_fixT(i)=.true. ; bc_valT(i)=T0
   endif
   if (y(i).gt.(Ly-eps) ) then 
      bc_fixT(i)=.true. ; bc_valT(i)=0.d0
      !write(*,*) "temp boundary set ", i
   endif
end do

!==============================================!
!=====[define bc for pressure]=================!
!==============================================!
bc_fixP=.false.
do i=1,nel
   do k=1,m
      if(y(icon(k,i)) .gt. (Ly-eps)) then
         bc_fixP(i)=.true. ; bc_valP = 0.d0
         !write(*,*) "Pressure boundary set ", i
      end if 
   end do 
end do 


!==============================================!
!=====[Initial temperature field]==============!
!==============================================!

do i=1,np
   T(i)=(1.d0-y(i)) - 0.01d0*cos(pi*x(i)/Lx)*sin(pi*y(i)/Ly)
end do

!==============================================!
!===========[Initial Pressure Field]===========!
!==============================================!

do i=1,iel
   press_el(i) = rho0*gy*(Lx-y(icon(1,i))+hx/2.d0)
   !this is a bit of a hack
   !hydrostatic pressure to initialise
end do



!********************************************************************************************
!********************************************************************************************
! T I M E S T E P P I N G 
!********************************************************************************************
!********************************************************************************************

do istep=1,nstep

print *,'***********************************************************'
print *,'*** istep= ',istep,'*** time=',time
print *,'***********************************************************'


!*********************************************************!

open(unit=666,file='conv_picard.dat',status='replace')

do ipicard=1,npicard ! picard iterations 

!=========Set up density ======================!
!open(unit=123,file='density.dat',status='replace')

!Get nodal pressure
call elemental_to_nodal(press_el,press_node,icon_inv,np,nel,m)

do i=1,np
  if ((ALA) .or. (BA) .or. (EBA)) then
   density(i)=rho0*(1.d0-alpha*T(i)+chi_T*press_node(i))
  else if(TALA) then
   density(i)=rho0*(1.d0-alpha*T(i))
  end if  

  !write(123,*) density(i)
end do
!close(123)


!==============================================!
!=====[build FE matrix]========================!
!==============================================!



KKK=0.d0
G=0.d0
rhs_f=0.d0
rhs_h=0.d0

do iel=1,nel

   Kel=0.d0
   fel=0.d0
   Gel=0.d0
   hel=0.d0

   do iq=-1,1,2
   do jq=-1,1,2

      rq=iq/sqrt(3.d0)
      sq=jq/sqrt(3.d0)
      wq=1.d0*1.d0

      N(1)=0.25d0*(1.d0-rq)*(1.d0-sq)
      N(2)=0.25d0*(1.d0+rq)*(1.d0-sq)
      N(3)=0.25d0*(1.d0+rq)*(1.d0+sq)
      N(4)=0.25d0*(1.d0-rq)*(1.d0+sq)

      dNdr(1)= - 0.25d0*(1.d0-sq)   ;   dNds(1)= - 0.25d0*(1.d0-rq)
      dNdr(2)= + 0.25d0*(1.d0-sq)   ;   dNds(2)= - 0.25d0*(1.d0+rq)
      dNdr(3)= + 0.25d0*(1.d0+sq)   ;   dNds(3)= + 0.25d0*(1.d0+rq)
      dNdr(4)= - 0.25d0*(1.d0+sq)   ;   dNds(4)= + 0.25d0*(1.d0-rq)

      jcb=0.d0
      do k=1,m
         jcb(1,1)=jcb(1,1)+dNdr(k)*x(icon(k,iel))
         jcb(1,2)=jcb(1,2)+dNdr(k)*y(icon(k,iel))
         jcb(2,1)=jcb(2,1)+dNds(k)*x(icon(k,iel))
         jcb(2,2)=jcb(2,2)+dNds(k)*y(icon(k,iel))
      enddo

      jcob=jcb(1,1)*jcb(2,2)-jcb(2,1)*jcb(1,2)

      jcbi(1,1)=    jcb(2,2) /jcob
      jcbi(1,2)=  - jcb(1,2) /jcob
      jcbi(2,1)=  - jcb(2,1) /jcob
      jcbi(2,2)=    jcb(1,1) /jcob

      xq=0.d0
      yq=0.d0
      uq=0.d0
      vq=0.d0
      exxq=0.d0
      eyyq=0.d0
      exyq=0.d0

      drhodx2=0.d0
      drhody2=0.d0
      rho2=0.d0
      do k=1,m
         xq=xq+N(k)*x(icon(k,iel))
         yq=yq+N(k)*y(icon(k,iel))
         uq=uq+N(k)*u(icon(k,iel))
         vq=vq+N(k)*v(icon(k,iel))
         dNdx(k)=jcbi(1,1)*dNdr(k)+jcbi(1,2)*dNds(k)
         dNdy(k)=jcbi(2,1)*dNdr(k)+jcbi(2,2)*dNds(k)
         exxq=exxq+ dNdx(k)*u(icon(k,iel))
         eyyq=eyyq+ dNdy(k)*v(icon(k,iel))
         exyq=exyq+ dNdx(k)*v(icon(k,iel)) *0.5d0 &
                  + dNdy(k)*u(icon(k,iel)) *0.5d0
         drhodx2=drhodx2 + dNdx(k)*density(icon(k,iel))
         drhody2=drhody2 + dNdy(k)*density(icon(k,iel))
         rho2=rho2 + N(k)*density(icon(k,iel))
      end do

      do i=1,m
         i1=2*i-1
         i2=2*i
         Bmat(1,i1)=Dinv2*dNdx(i) ; Bmat(1,i2)=Dinv1*dNdy(i)
         Bmat(2,i1)=Dinv1*dNdx(i) ; Bmat(2,i2)=Dinv2*dNdy(i)
         Bmat(3,i1)=dNdy(i)       ; Bmat(3,i2)=dNdx(i)
         ! Bmat(1,i1)=+0.5*dNdx(i) ; Bmat(1,i2)=-0.5*dNdy(i)
         ! Bmat(2,i1)=-0.5*dNdx(i) ; Bmat(2,i2)=+0.5*dNdy(i)
         ! Bmat(3,i1)=dNdy(i)      ; Bmat(3,i2)=dNdx(i)
      end do

      Kel=Kel + matmul(transpose(Bmat),matmul(viscosity*Cmat,Bmat))*wq*jcob

      do i=1,m
      i1=2*i-1
      i2=2*i

      fel(i1)=fel(i1)+N(i)*jcob*wq*rho2*gx
      fel(i2)=fel(i2)+N(i)*jcob*wq*rho2*gy
      Gel(i1,1)=Gel(i1,1)-dNdx(i)*jcob*wq
      Gel(i2,1)=Gel(i2,1)-dNdy(i)*jcob*wq
      end do
      !hel(1)=hel(1)+1.d0/rho(xq,yq,ibench)*(uq* drhodx(xq,yq,ibench) + vq* drhody(xq,yq,ibench) )*wq*jcob ! approach 1
      if (compressible) then
      hel(1)=hel(1)+1.d0/rho2*(uq*drhodx2+vq*drhody2)*wq*jcob !approach 2
      end if 

   end do
   end do

   !=====================
   !=====[assemble]======
   !=====================

   do k1=1,m
      ik=icon(k1,iel)
      do i1=1,ndofV
         ikk=ndofV*(k1-1)+i1
         m1=ndofV*(ik-1)+i1
         do k2=1,m
            jk=icon(k2,iel)
            do i2=1,ndofV
               jkk=ndofV*(k2-1)+i2
               m2=ndofV*(jk-1)+i2
               KKK(m1,m2)=KKK(m1,m2)+Kel(ikk,jkk)
            end do
         end do
         rhs_f(m1)=rhs_f(m1)+fel(ikk)
         G(m1,iel)=G(m1,iel)+Gel(ikk,1)
      end do
   end do

   rhs_h(iel)=hel(1)

end do

!==============================================!
!=====[impose b.c.]============================!
!==============================================!

do i=1,NfemV
    if (bc_fixV(i)) then 
      Aref=KKK(i,i)
      do j=1,NfemV
         rhs_f(j)=rhs_f(j)-KKK(i,j)*bc_valV(i)
         KKK(i,j)=0.d0
         KKK(j,i)=0.d0
      enddo
      KKK(i,i)=Aref
      rhs_f(i)=Aref*bc_valV(i)

      !G(i,:)=0

      do j=1,nel
         rhs_h(j)=rhs_h(j)-G(i,j)*bc_valV(i)
         G(i,j)=0.d0
      end do

   endif
enddo





! print *,minval(KKK),maxval(KKK)
! print *,minval(G),maxval(G)
! print *,minval(rhs_f),maxval(rhs_f)
! print *,minval(rhs_h),maxval(rhs_h)

!==============================================!
!=====[solve system]===========================!
!==============================================!

KKKmem=KKK

!write(*,*) ipicard, maxval(abs(KKKmem)),maxval(abs(Vsol)),maxval(abs(G)),maxval(abs(Psol)),maxval(abs(rhs_f))


select case(solver)
case(1)
call solve_uzawa1(KKK,G,NfemV,nel,rhs_f,rhs_h,Vsol,press_el)
case(2)
call solve_uzawa2(KKK,G,NfemV,nel,rhs_f,rhs_h,Vsol,press_el)
case(3)
stop 'uzawa3 does not work'
call solve_uzawa3(KKK,G,NfemV,nel,rhs_f,rhs_h,Vsol,press_el)
case(4)
call solve_full(KKK,G,NfemV,nel,rhs_f,rhs_h,Vsol,press_el)
end select

do i=1,np
   u(i)=Vsol((i-1)*ndofV+1)
   v(i)=Vsol((i-1)*ndofV+2)
end do
Rv=matmul(KKKmem,Vsol)+matmul(G,press_el)-rhs_f 
Rp=matmul(transpose(G),Vsol)-rhs_h          

!write(*,*) ipicard, maxval(abs(Rv))/maxval(abs(rhs_f)), maxval(abs(Rp))/maxval(abs(rhs_h)),maxval(abs(KKKmem)),maxval(abs(Vsol)),maxval(abs(G)),maxval(abs(Psol)),maxval(abs(rhs_f))

!==============================================!
!=====[iterations converged?]==================!
!==============================================!

V_diff=maxval(abs(Vsol-Vsolmem))/maxval(abs(Vsol))

P_diff=maxval(abs(press_el-Psolmem))/maxval(abs(press_el))

Vsolmem=Vsol
Psolmem=press_el


write(666,*) ipicard,V_diff,P_diff

write(*,'(a,i3,a,3es15.6,a,2es15.6)') &
'Picard it. ',ipicard,&
' | <|u|>,<|v|>,<|p|> ',sum(abs(u))/np,sum(abs(v))/np,sum(abs(press_el))/nel,&
' | V_diff, P_diff ',V_diff,P_diff

if (V_diff<tol .and. P_diff<tol) then
print *,'-> Picard iterations have converged'
exit
end if


end do 

close(666)



!==============================================!
!========[Pressure Boundary Conditions]========!
!==============================================!

mean_surface_pressure = 0.d0
number_of_surface_nodes = 0
do i=1,iel
  if (bc_fixP(i)) then
    !write(*,*) i
    mean_surface_pressure = mean_surface_pressure + press_el(i)
    number_of_surface_nodes = number_of_surface_nodes + 1
  end if 
end do 
!write(*,*) number_of_surface_nodes
mean_surface_pressure = mean_surface_pressure/real(number_of_surface_nodes)
!write(*,*) mean_surface_pressure
do i=1,iel
  press_el(i) =  press_el(i) - mean_surface_pressure
end do 

!=====[Compute vrms]===========================!
!==============================================!

vrms=0.d0

iel=0
do j=1,nely
do i=1,nelx
   iel=iel+1

   do iq=-1,1,2
   do jq=-1,1,2

      rq=iq/sqrt(3.d0)
      sq=jq/sqrt(3.d0)
      wq=1.d0*1.d0

      N(1)=0.25d0*(1.d0-rq)*(1.d0-sq)
      N(2)=0.25d0*(1.d0+rq)*(1.d0-sq)
      N(3)=0.25d0*(1.d0+rq)*(1.d0+sq)
      N(4)=0.25d0*(1.d0-rq)*(1.d0+sq)

      dNdr(1)= - 0.25d0*(1.d0-sq)   ;   dNds(1)= - 0.25d0*(1.d0-rq)
      dNdr(2)= + 0.25d0*(1.d0-sq)   ;   dNds(2)= - 0.25d0*(1.d0+rq)
      dNdr(3)= + 0.25d0*(1.d0+sq)   ;   dNds(3)= + 0.25d0*(1.d0+rq)
      dNdr(4)= - 0.25d0*(1.d0+sq)   ;   dNds(4)= + 0.25d0*(1.d0-rq)

      jcb=0.d0
      do k=1,m
         jcb(1,1)=jcb(1,1)+dNdr(k)*x(icon(k,iel))
         jcb(1,2)=jcb(1,2)+dNdr(k)*y(icon(k,iel))
         jcb(2,1)=jcb(2,1)+dNds(k)*x(icon(k,iel))
         jcb(2,2)=jcb(2,2)+dNds(k)*y(icon(k,iel))
      enddo

      jcob=jcb(1,1)*jcb(2,2)-jcb(2,1)*jcb(1,2)

      uq=0.d0
      vq=0.d0
      do k=1,m
         uq=uq+N(k)*u(icon(k,iel))
         vq=vq+N(k)*v(icon(k,iel))
      end do

      vrms=vrms+(uq**2+vq**2)*jcob*wq

   end do
   end do

end do
end do

vrms=sqrt(vrms/Lx/Ly)

write(1003,*) time,vrms ; call flush(1003)

!==============================================!
!=====[Compute timestep]=======================!
!==============================================!

umax=maxval(abs(u))
vmax=maxval(abs(v))

dt=courant_nb*min(hx,hy)/max(umax,vmax)

time=time+dt

write(1005,*) istep,dt ; call flush(1005)

!==============================================!
!=====[Build temperature matrix]===============!
!==============================================!



call elemental_to_nodal(press_el,press_node,icon_inv,np,nel,m)

do i=1,np
  if ((ALA) .or. (BA) .or. (EBA)) then
   density(i)=rho0*(1.d0-alpha*T(i)+chi_T*press_node(i))
  else if(TALA) then
   density(i)=rho0*(1.d0-alpha*T(i))
  end if  
end do

A_T=0.d0
B_T=0.d0

do iel=1,nel

   AelT=0.d0
   BelT=0.d0

   temp(1)=T(icon(1,iel))
   temp(2)=T(icon(2,iel))
   temp(3)=T(icon(3,iel))
   temp(4)=T(icon(4,iel))

   do iq=-1,1,2
   do jq=-1,1,2

      rq=iq/sqrt(3.d0)
      sq=jq/sqrt(3.d0)
      wq=1.d0*1.d0

      N(1)=0.25d0*(1.d0-rq)*(1.d0-sq)
      N(2)=0.25d0*(1.d0+rq)*(1.d0-sq)
      N(3)=0.25d0*(1.d0+rq)*(1.d0+sq)
      N(4)=0.25d0*(1.d0-rq)*(1.d0+sq)

      dNdr(1)= - 0.25d0*(1.d0-sq)   ;   dNds(1)= - 0.25d0*(1.d0-rq)
      dNdr(2)= + 0.25d0*(1.d0-sq)   ;   dNds(2)= - 0.25d0*(1.d0+rq)
      dNdr(3)= + 0.25d0*(1.d0+sq)   ;   dNds(3)= + 0.25d0*(1.d0+rq)
      dNdr(4)= - 0.25d0*(1.d0+sq)   ;   dNds(4)= + 0.25d0*(1.d0-rq)

      jcb=0.d0
      do k=1,m
         jcb(1,1)=jcb(1,1)+dNdr(k)*x(icon(k,iel))
         jcb(1,2)=jcb(1,2)+dNdr(k)*y(icon(k,iel))
         jcb(2,1)=jcb(2,1)+dNds(k)*x(icon(k,iel))
         jcb(2,2)=jcb(2,2)+dNds(k)*y(icon(k,iel))
      enddo

      jcob=jcb(1,1)*jcb(2,2)-jcb(2,1)*jcb(1,2)

      jcbi(1,1)=    jcb(2,2) /jcob
      jcbi(1,2)=  - jcb(1,2) /jcob
      jcbi(2,1)=  - jcb(2,1) /jcob
      jcbi(2,2)=    jcb(1,1) /jcob

      xq=0.d0
      yq=0.d0
      uq=0.d0
      vq=0.d0
      phi=0.d0
      rho_el=0.d0
      v_el=0.d0
      do k=1,m
         xq=xq+N(k)*x(icon(k,iel))
         yq=yq+N(k)*y(icon(k,iel))
         uq=uq+N(k)*u(icon(k,iel))
         vq=vq+N(k)*v(icon(k,iel))
         dNdx(k)=jcbi(1,1)*dNdr(k)+jcbi(1,2)*dNds(k)
         dNdy(k)=jcbi(2,1)*dNdr(k)+jcbi(2,2)*dNds(k)

         phi = phi + viscosity*(2.d0*(dNdx(k)*u(icon(k,iel)))**2 - 2/D*(dNdx(k)*u(icon(k,iel)) + dNdx(k)*u(icon(k,iel)))**2)
         phi = phi + viscosity*(2.d0*(dNdy(k)*y(icon(k,iel)))**2 - 2/D*(dNdx(k)*u(icon(k,iel)) + dNdx(k)*u(icon(k,iel)))**2)
         phi = phi + viscosity*((dNdx(k)*v(icon(k,iel)))**2 + dNdx(k)*u(icon(k,iel))*dNdy(k)*v(icon(k,iel)))
         phi = phi + viscosity*((dNdy(k)*u(icon(k,iel)))**2 + dNdx(k)*u(icon(k,iel))*dNdy(k)*v(icon(k,iel)))

         rho_el = rho_el + density(icon(k,iel))*N(k)
         v_el = v_el + v(icon(k,iel))*N(k)
      end do
      phi = phi*wq*jcob
      velnorm=sqrt(uq**2+vq**2)

      vel2D(1,1)=uq
      vel2D(1,2)=vq

      Nvect(1,:)=N(:)
      NvectT(:,1)=N(:)

      BmatT(1,:)=dNdx
      BmatT(2,:)=dNdy

      BmatTT=transpose(BmatT)

      if (buoyancy) then
        Kb=-matmul(NvectT,Nvect)*wq*jcob*alpha*gy*v_el*rho_el*0.d0
      else
        Kb = 0.d0
      end if 

      Ka=matmul(NvectT,matmul(vel2D,BmatT))*wq*jcob*rho0*hcapa

      Kc=matmul(BmatTT,BmatT)*wq*jcob*hcond

      KK=(Ka+Kc+Kb)

      M_T=matmul(NvectT,Nvect)*wq*jcob*rho0*hcapa

      AelT=AelT+(M_T+KK*theta*dt)

      BelT=BelT+matmul(M_T-KK*(1.d0-theta)*dt,temp(1:4))

      do k=1,m
        !BelT(k) = BelT(k) + phi
      end do 

   end do
   end do

   ! impose boundary conditions

   do i=1,m
      ik=icon(i,iel)
      if (bc_fixT(ik)) then
         Aref=AelT(i,i)
         do j=1,m
         BelT(j)=BelT(j)-AelT(j,i)*bc_valT(ik)
         AelT(i,j)=0.d0
         AelT(j,i)=0.d0
         enddo
         AelT(i,i)=Aref
         BelT(i)=Aref*bc_valT(ik)
      endif
   enddo

   !=====================
   !=====[assemble]======
   !=====================

   do k1=1,m
      ik=icon(k1,iel)
      do i1=1,ndofT
         ikk=ndofT*(k1-1)+i1
         m1=ndofT*(ik-1)+i1
         do k2=1,m
            jk=icon(k2,iel)
            do i2=1,ndofT
               jkk=ndofT*(k2-1)+i2
               m2=ndofT*(jk-1)+i2
               A_T(m1,m2)=A_T(m1,m2)+AelT(ikk,jkk)
            end do
         end do
         B_T(m1)=B_T(m1)+BelT(ikk)
      end do
   end do

end do ! end of loop over cells

!print *,minval(B_T),maxval(B_T)

!==============================================!
!=====[solve system]===========================!
!==============================================!

job=0
allocate(work(NfemT))
allocate(ipvt(NfemT))
call DGECO (A_T, NfemT, NfemT, ipvt, rcond, work)
call DGESL (A_T, NfemT, NfemT, ipvt, B_T, job)
deallocate(ipvt)
deallocate(work)

do i=1,np
   T(i)=B_T(i)
end do

write(*,*) 'min/max T',minval(T),maxval(T),sum(T)

open(unit=123,file='OUT/solution_T.dat',status='replace')
do i=1,np
   write(123,'(5f20.10)') x(i),y(i),u(i),T(i),T(i)-(1-y(i))
end do
close(123)

write(1002,*) time,minval(T),maxval(T)

!==============================================!
!=====[Compute Nusselt number]=================!
!==============================================!

qx=0
qy=0
Nu=0.d0

iel=0
do j=1,nely
do i=1,nelx
   iel=iel+1

   rq=0.d0
   sq=0.d0
      
   N(1)=0.25d0*(1.d0-rq)*(1.d0-sq)
   N(2)=0.25d0*(1.d0+rq)*(1.d0-sq)
   N(3)=0.25d0*(1.d0+rq)*(1.d0+sq)
   N(4)=0.25d0*(1.d0-rq)*(1.d0+sq)

   dNdr(1)= - 0.25d0*(1.d0-sq)   ;   dNds(1)= - 0.25d0*(1.d0-rq)
   dNdr(2)= + 0.25d0*(1.d0-sq)   ;   dNds(2)= - 0.25d0*(1.d0+rq)
   dNdr(3)= + 0.25d0*(1.d0+sq)   ;   dNds(3)= + 0.25d0*(1.d0+rq)
   dNdr(4)= - 0.25d0*(1.d0+sq)   ;   dNds(4)= + 0.25d0*(1.d0-rq)

   jcb=0.d0
   do k=1,m
      jcb(1,1)=jcb(1,1)+dNdr(k)*x(icon(k,iel))
      jcb(1,2)=jcb(1,2)+dNdr(k)*y(icon(k,iel))
      jcb(2,1)=jcb(2,1)+dNds(k)*x(icon(k,iel))
      jcb(2,2)=jcb(2,2)+dNds(k)*y(icon(k,iel))
   enddo

   jcob=jcb(1,1)*jcb(2,2)-jcb(2,1)*jcb(1,2)

   jcbi(1,1)=    jcb(2,2) /jcob
   jcbi(1,2)=  - jcb(1,2) /jcob
   jcbi(2,1)=  - jcb(2,1) /jcob
   jcbi(2,2)=    jcb(1,1) /jcob

   do k=1,m
      dNdx(k)=jcbi(1,1)*dNdr(k)+jcbi(1,2)*dNds(k)
      dNdy(k)=jcbi(2,1)*dNdr(k)+jcbi(2,2)*dNds(k)
   end do

   dTdxq=0.d0
   dTdyq=0.d0
   do k=1,m
      dTdxq=dTdxq+ dNdx(k)*T(icon(k,iel))
      dTdyq=dTdyq+ dNdy(k)*T(icon(k,iel))
   end do

   if (j==nely) then
      Nu=Nu+abs(dTdyq)*hx
   end if

   qx(iel)=-hcond*dTdxq
   qy(iel)=-hcond*dTdyq

end do
end do

write(1004,*) time,Nu ; call flush(1004)

write(*,*) 'Ra=',Ra,'Nu=',Nu


!==============================================!
!=====[compute nodal density for plotting]=====! 
!==============================================!

do i=1,np
   density(i)=rho0*(1.d0-alpha*T(i))
end do


call elemental_to_nodal(press_el,press_node,icon_inv,np,nel,m)

! ==============================================!
! ====[compute the total mass of the system]====!
! ==============================================!
mass=0.d0
gpe=0.d0
ke=0.d0
thermal_energy=0.d0
ee=0.d0
pe=0.d0
pe2 = 0.d0
do iel=1,nel

   do iq=-1,1,2
   do jq=-1,1,2

      rq=iq/sqrt(3.d0)
      sq=jq/sqrt(3.d0)
      wq=1.d0*1.d0

      N(1)=0.25d0*(1.d0-rq)*(1.d0-sq)
      N(2)=0.25d0*(1.d0+rq)*(1.d0-sq)
      N(3)=0.25d0*(1.d0+rq)*(1.d0+sq)
      N(4)=0.25d0*(1.d0-rq)*(1.d0+sq)

      dNdr(1)= - 0.25d0*(1.d0-sq)   ;   dNds(1)= - 0.25d0*(1.d0-rq)
      dNdr(2)= + 0.25d0*(1.d0-sq)   ;   dNds(2)= - 0.25d0*(1.d0+rq)
      dNdr(3)= + 0.25d0*(1.d0+sq)   ;   dNds(3)= + 0.25d0*(1.d0+rq)
      dNdr(4)= - 0.25d0*(1.d0+sq)   ;   dNds(4)= + 0.25d0*(1.d0-rq)

      jcb=0.d0
      do k=1,m
         jcb(1,1)=jcb(1,1)+dNdr(k)*x(icon(k,iel))
         jcb(1,2)=jcb(1,2)+dNdr(k)*y(icon(k,iel))
         jcb(2,1)=jcb(2,1)+dNds(k)*x(icon(k,iel))
         jcb(2,2)=jcb(2,2)+dNds(k)*y(icon(k,iel))
      enddo

      jcob=jcb(1,1)*jcb(2,2)-jcb(2,1)*jcb(1,2)

      jcbi(1,1)=    jcb(2,2) /jcob
      jcbi(1,2)=  - jcb(1,2) /jcob
      jcbi(2,1)=  - jcb(2,1) /jcob
      jcbi(2,2)=    jcb(1,1) /jcob

      do k=1,m
         dNdx(k)=jcbi(1,1)*dNdr(k)+jcbi(1,2)*dNds(k)
         dNdy(k)=jcbi(2,1)*dNdr(k)+jcbi(2,2)*dNds(k)
      end do

      rho2=0.d0
      ke_density=0.d0
      gpe_density=0.d0
      thermal_energy_density=0.d0
      ee_density = 0.d0
      pe_density=0.d0
      do k=1,m
         rho2=rho2 + N(k)*density(icon(k,iel))
         ke_density = ke_density + 0.5*N(k)*density(icon(k,iel))*(u(icon(k,iel))**2 + v(icon(k,iel))**2)
         gpe_density=gpe_density + N(k)*density(icon(k,iel))*(Ly-y(icon(k,iel)))*gy
         thermal_energy_density = thermal_energy_density + hcapa*T(icon(k,iel))*N(k)

         !ee_density = ee_density + 2.d0*viscosity*((dNdx(k)*u(icon(k,iel)))**2 +(dNdy(k)*v(icon(k,iel)))**2 )
         !ee_density = ee_density + viscosity*(dNdx(k)*v(icon(k,iel)) + dNdy(k)*u(icon(k,iel)) )**2 

         !This is the version according to Rock Mechanics Book
         ee_density = ee_density + viscosity*(2.d0*dNdx(k)*u(icon(k,iel)) -2.d0/D*(dNdx(k)*u(icon(k,iel)) + dNdy(k)*v(icon(k,iel))))*dNdx(k)*u(icon(k,iel))
         ee_density = ee_density + viscosity*(2.d0*dNdy(k)*v(icon(k,iel)) -2.d0/D*(dNdx(k)*u(icon(k,iel)) + dNdy(k)*v(icon(k,iel))))*dNdy(k)*v(icon(k,iel))
         ee_density = ee_density + 2.d0*viscosity*(dNdx(k)*v(icon(k,iel)) + dNdy(k)*u(icon(k,iel)))*0.5*(dNdx(k)*v(icon(k,iel)) + dNdy(k)*u(icon(k,iel)))

         pe_density = pe_density + press_node(icon(k,iel))*N(k)
      end do
      mass = mass+rho2*wq*jcob
      ke = ke + ke_density*wq*jcob
      gpe=gpe+gpe_density*wq*jcob
      thermal_energy = thermal_energy + thermal_energy_density*wq*jcob
      ee = ee + ee_density*wq*jcob
      pe = pe + pe_density*wq*jcob
      pe2 = pe2 + press_el(iel)*wq*jcob
   end do
   end do

end do

! write(*,*) 'Mass =', mass
! write(*,*) 'GPE  =', gpe
! write(*,*) 'Ke   =', ke
!write(1008,*) mass, ' ', gpe, ' ',ke, ' ', thermal_energy, ' ', ee; call flush(1008)



!==============================================!
!========[Compute Heat Flux]===================!
!==============================================!


phi_q_top=0.d0
phi_q_bottom=0.d0

do iel=1,nel

node_is_boundary = .false.
  do k=1,m
   if (bc_fixT(icon(k,iel))) then
      node_is_boundary = .true.
      element_boundary_value = bc_valT(icon(k,iel))
   end if 
  end do 

if(node_is_boundary) then
    iq=0
    jq=0

      rq=iq/sqrt(3.d0)
      sq=jq/sqrt(3.d0)
      wq=1.d0*1.d0

      N(1)=0.25d0*(1.d0-rq)*(1.d0-sq)
      N(2)=0.25d0*(1.d0+rq)*(1.d0-sq)
      N(3)=0.25d0*(1.d0+rq)*(1.d0+sq)
      N(4)=0.25d0*(1.d0-rq)*(1.d0+sq)

      dNdr(1)= - 0.25d0*(1.d0-sq)   ;   dNds(1)= - 0.25d0*(1.d0-rq)
      dNdr(2)= + 0.25d0*(1.d0-sq)   ;   dNds(2)= - 0.25d0*(1.d0+rq)
      dNdr(3)= + 0.25d0*(1.d0+sq)   ;   dNds(3)= + 0.25d0*(1.d0+rq)
      dNdr(4)= - 0.25d0*(1.d0+sq)   ;   dNds(4)= + 0.25d0*(1.d0-rq)

      jcb=0.d0
      do k=1,m
         jcb(1,1)=jcb(1,1)+dNdr(k)*x(icon(k,iel))
         jcb(1,2)=jcb(1,2)+dNdr(k)*y(icon(k,iel))
         jcb(2,1)=jcb(2,1)+dNds(k)*x(icon(k,iel))
         jcb(2,2)=jcb(2,2)+dNds(k)*y(icon(k,iel))
      enddo

      jcob=jcb(1,1)*jcb(2,2)-jcb(2,1)*jcb(1,2)

      jcbi(1,1)=    jcb(2,2) /jcob
      jcbi(1,2)=  - jcb(1,2) /jcob
      jcbi(2,1)=  - jcb(2,1) /jcob
      jcbi(2,2)=    jcb(1,1) /jcob
      phi_q_el = 0.d0
      do k=1,m
         dNdx(k)=jcbi(1,1)*dNdr(k)+jcbi(1,2)*dNds(k)
         dNdy(k)=jcbi(2,1)*dNdr(k)+jcbi(2,2)*dNds(k)
         !phi_q_el = phi_q_el - hcond*T(icon(k,iel))*(dNdr(k) + dNds(k))*wq*jcob
      end do

      do k=1,m
        phi_q_el = phi_q_el + dNdy(k)*T(icon(k,iel))*hcond*hcapa*hx
         
      end do 



      if(element_boundary_value .eq. 0.d0) then
         phi_q_top = phi_q_top + phi_q_el
      else if(element_boundary_value .eq. T0) then
         phi_q_bottom = phi_q_bottom + phi_q_el
      else
         write(*,*) "This should not be printed. Check boundaries in flux calculations."
      end if




end if 
end do

!=================================================!
!====[Compute the work done by the boundaries]====!
!=================================================!

work_energy = 0.d0
do iel=1,nel
  node_is_boundary = .false.
  normal_x = 0.d0
  normal_y = 0.d0

  !Calculate if an element is adjacent to a boundary
  !and if so, the normal to the boundary
  if (iel .le. nnx) then
    node_is_boundary = .true.
    normal_y = -1.d0
  else if (iel .gt. np-nnx) then
    node_is_boundary = .true.
    normal_y = 1.d0
  end if 
  !Implement as two if statements for corner elements
  if (modulo(iel,nnx) .eq. 0)  then
    node_is_boundary = .true.
    normal_x = 1.d0
  else if (modulo(iel-1,nnx) .eq. 0) then
    node_is_boundary = .true.
    normal_x = -1.d0
  end if 

  !normalise normals
  if ((normal_x .ne. 0.d0) .and. (normal_y .ne. 0.d0)) then
    normal_x = normal_x/sqrt(2.d0)
    normal_y = normal_y/sqrt(2.d0)
  end if 

  if(node_is_boundary) then
    iq=0
    jq=0

      rq=iq/sqrt(3.d0)
      sq=jq/sqrt(3.d0)
      wq=1.d0*1.d0

      N(1)=0.25d0*(1.d0-rq)*(1.d0-sq)
      N(2)=0.25d0*(1.d0+rq)*(1.d0-sq)
      N(3)=0.25d0*(1.d0+rq)*(1.d0+sq)
      N(4)=0.25d0*(1.d0-rq)*(1.d0+sq)

      dNdr(1)= - 0.25d0*(1.d0-sq)   ;   dNds(1)= - 0.25d0*(1.d0-rq)
      dNdr(2)= + 0.25d0*(1.d0-sq)   ;   dNds(2)= - 0.25d0*(1.d0+rq)
      dNdr(3)= + 0.25d0*(1.d0+sq)   ;   dNds(3)= + 0.25d0*(1.d0+rq)
      dNdr(4)= - 0.25d0*(1.d0+sq)   ;   dNds(4)= + 0.25d0*(1.d0-rq)

      jcb=0.d0
      do k=1,m
         jcb(1,1)=jcb(1,1)+dNdr(k)*x(icon(k,iel))
         jcb(1,2)=jcb(1,2)+dNdr(k)*y(icon(k,iel))
         jcb(2,1)=jcb(2,1)+dNds(k)*x(icon(k,iel))
         jcb(2,2)=jcb(2,2)+dNds(k)*y(icon(k,iel))
      enddo

      jcob=jcb(1,1)*jcb(2,2)-jcb(2,1)*jcb(1,2)

      jcbi(1,1)=    jcb(2,2) /jcob
      jcbi(1,2)=  - jcb(1,2) /jcob
      jcbi(2,1)=  - jcb(2,1) /jcob
      jcbi(2,2)=    jcb(1,1) /jcob
      do k=1,m
         dNdx(k)=jcbi(1,1)*dNdr(k)+jcbi(1,2)*dNds(k)
         dNdy(k)=jcbi(2,1)*dNdr(k)+jcbi(2,2)*dNds(k)
      end do
      work_el = 0.d0
      do k=1,m
        !term xx
        work_el = work_el + (viscosity*(2*dNdx(k)*u(icon(k,iel)) -2/D*(dNdx(k)*u(icon(k,iel))+dNdy(k)*v(icon(k,iel)))) + press_node(icon(k,iel)))*normal_x*u(icon(k,iel))
        !term xy
        work_el = work_el + (viscosity*(dNdx(k)*v(icon(k,iel)) + dNdy(k)*u(icon(k,iel)))  + press_node(icon(k,iel)))*normal_y*v(icon(k,iel))
        !term yx
        work_el = work_el + (viscosity*(dNdx(k)*v(icon(k,iel)) + dNdy(k)*u(icon(k,iel)))  + press_node(icon(k,iel)))*normal_x*u(icon(k,iel))
          !term yy
        work_el = work_el + (viscosity*(2*dNdy(k)*v(icon(k,iel)) -2/D*(dNdx(k)*u(icon(k,iel))+dNdy(k)*v(icon(k,iel)))) + press_node(icon(k,iel)))*normal_y*v(icon(k,iel))
      end do 
      work_energy = work_energy + work_el*hx

  end if 
end do 





write(1008,*) mass, ' ', gpe, ' ',ke, ' ', thermal_energy, ' ', ee, ' ', pe,' ', pe2,' ', phi_q_bottom, ' ', phi_q_top, ' ', work_energy; call flush(1008)


!==============================================!
!=====[output data in vtu format]==============!
!==============================================!

if (istep==1 .or. mod(istep,output_freq)==0) &
call output_for_paraview (np,nel,x,y,u,v,press_el,T,density,qx,qy,icon,istep)

!==============================================!
!====[is steady state reached ?]===============!
!==============================================!

chi_Nu=abs(Nu-Nu_prev)/Nu
chi_vrms=abs(vrms-vrms_prev)/vrms 

write(1006,*) time,chi_Nu   ; call flush(1006)
write(1007,*) time,chi_vrms ; call flush(1007)

write(*,*) 'chi_Nu  =',chi_Nu,'tol=',tol
write(*,*) 'chi_vrms=',chi_vrms,'tol=',tol

if (chi_Nu<tol .and. chi_vrms<tol) then

   ! output avrg temperature profile
   counter=0
   do j=1,nny
   do i=1,nnx
      counter=counter+1
      Tavrg(j)=Tavrg(j)+T(counter)
   end do
   end do
   Tavrg=Tavrg/nnx
   open (unit=888,file='OUT/Tavrg.dat')
   do j=1,nny
   write(888,*) y(j*nnx),Tavrg(j)
   end do
   close(888)

   print *,'CONVERGED'
   exit
end if

Nu_prev=Nu
vrms_prev=vrms

end do ! timestepping

!********************************************************************************************
!********************************************************************************************

print *,'***********************************************************'

end program

