!==============================================!
!                                              !
! C. thieulot ; March 2016                     !
!                                              !
!==============================================!
                                               !
!Command line argument documentation           !
!Argument 1 : nnx                              !
!Argument 2 : approximation - ALA,TALA,BA,EBA  !
!Argument 3 : power (of Ra)                    !
!Argument 4 : Di                               !
!Argument 5 : MUMPS_T                          !
!Argument 5 : MUMPS_VP                         !


program fcubed                                 !
                                               !
implicit none                                  !

include 'mpif.h'                               !   
include 'dmumps_struc.h'      
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

integer ierr,iii,inode,iproc,nproc,idof        !
integer LELTVAR,NA_ELT,counter_mumps,ii,ij     !

                                               !  
real(8) Lx,Ly                                  ! size of the numerical domain
real(8) viscosity                              ! dynamic viscosity $\mu$ of the material
real(8) rho0                                   ! reference mass density $\rho_0$ of the material
real(8) gx,gy                                  ! gravity acceleration
real(8), dimension(:),   allocatable :: x,y    ! node coordinates arrays
real(8), dimension(:),   allocatable :: u,v    ! node velocity arrays
real(8), dimension(:),   allocatable :: u_prev,v_prev! node velocity arrays - previous timestep
real(8), dimension(:),   allocatable :: press_el !pressure , elemental
real(8), dimension(:),   allocatable :: press_node !also pressure, but nodal 
real(8), dimension(:),   allocatable :: BV     ! right hand side of Stokes system
real(8), dimension(:,:), allocatable :: AV     ! FEM matrix of Stokes system
real(8), dimension(:),   allocatable :: B_T    ! right hand side of temperature system
real(8), dimension(:,:), allocatable :: A_T    ! FEM matrix of temperature system
real(8), dimension(:),   allocatable :: work   ! work array needed by the solver
real(8), dimension(:),   allocatable :: bc_valV! array containing bc values for velocity
real(8), dimension(:),   allocatable :: bc_valT! array containing bc values for temperature
real(8), dimension(:),   allocatable :: bc_valP! array containing bc values for pressure
real(8), dimension(:),   allocatable :: T      ! node temperature array
real(8), dimension(:),   allocatable :: T_prev ! previous node temperature array
real(8), dimension(:),   allocatable :: density! node density array
real(8), dimension(:),   allocatable :: qx,qy  ! heat flux vector 
real(8), dimension(:),   allocatable :: Tavrg  ! average temperature
                                               !
real(8) rq,sq,wq                               ! local coordinate and weight of qpoint
real(8) xq,yq                                  ! global coordinate of qpoint
real(8) uq,vq                                  ! velocity at qpoint
real(8) exxq,eyyq,exyq                         ! strain-rate components at qpoint  
real(8) AelT(m*ndofT,m*ndofT)                  ! elemental temperature FEM matrix
real(8) BelT(m*ndofT)                          ! elemental temperature right hand side
real(8) N(m),dNdx(m),dNdy(m),dNdr(m),dNds(m)   ! shape fcts and derivatives
real(8) jcob                                   ! determinant of jacobian matrix
real(8) jcb(2,2)                               ! jacobian matrix
real(8) jcbi(2,2)                              ! inverse of jacobian matrix
real(8), dimension(3,3) :: Kmat                ! K matrix 
real(8), dimension(3,3) :: Cmat                ! C matrix
real(8) Aref                                   !
real(8) hcapa                                  ! heat capacity
real(8) hcond                                  ! heat conductivity
real(8) alpha                                  ! thermal expansion 
real(8) gamma                                  ! Gruneissen Parameter
real(8) time                                   ! real time
real(8) dt                                     ! timestep
real(8) courant_nb                             ! Courant number for CFL criterion
real(8) Ra                                     ! Rayleigh number 
real(8) chi_T
real(8), parameter :: theta = 0.5              ! mid-point timestepping parameter
real(8), parameter :: eps=1.d-10               !
real(8), parameter :: pi = 3.14159265359d0     !
real(8), parameter :: tol=1.d-6                ! convergence tolerance
real(8), parameter :: tol_nu = 1.d-7
real(8) M_T(4,4),Ka(4,4),Kc(4,4),KK(4,4),Kb(4,4) ! various FEM arrays
real(8) Nvect(1,4),NvectT(4,1)                 !
real(8) velnorm,vel2D(1,2),temp(4),temp_2(4)   !
real(8) BmatT(2,4),BmatTT(4,2)                 !
real(8) hx,hy                                  ! grid spacing
real(8) umax,vmax                              !
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

integer ipicard,npicard                        ! Loop variables for picard iteration
integer solver                                 ! Switch for which solver to use 
real(8) P_diff,V_diff,V_diff_mem,V_diff_mem_mem! Variables and arrays needed
real(8), dimension(:),   allocatable :: Vsol   ! by the solver(s)
real(8), dimension(:),   allocatable :: Vsolmem!
real(8), dimension(:),   allocatable :: Psolmem! 
real(8), dimension(:),   allocatable :: Rv     !
real(8), dimension(:),   allocatable :: Rp     !

integer, dimension(:), allocatable :: element_boundary_value ! 
logical, dimension(:), allocatable :: element_is_boundary    !
real(8), dimension(:), allocatable :: dxu_elemental,dxv_elemental,dyu_elemental,dyv_elemental
real(8), dimension(:), allocatable :: dxu_nodal,dxv_nodal,dyu_nodal,dyv_nodal
real(8), dimension(:), allocatable :: dxu_deviatoric,dxv_deviatoric,dyu_deviatoric,dyv_deviatoric, trace_strain_nodal
real(8), dimension(:), allocatable :: phi_nodal,phi_elemental,adiabatic_heating
real(8) str_inv2,e_kk

real(8) mean_surface_temp                      ! Mean surface temperature
real(8), dimension(:),   allocatable :: dTdx_nodal,dTdy_nodal
real(8), dimension(:),   allocatable :: dTdx_elemental, dTdy_elemental
real(8) T0, Di                                 ! Basal temperature, dissipation number
real(8) mean_surface_pressure                  !
integer number_of_surface_nodes                !
real(8),  parameter :: D = 3.d0                ! parameter determining dimensionality of simulation
real(8),  parameter :: Dinv1 = -1.d0/D         ! allowing switching between 2D/3D
real(8) phi,rho_el,v_el                        ! 
logical compressible, ALA, TALA, BA, EBA,buoyancy! Booleans for approximations
integer approximation                          ! Switch for which approximation to use
real(8) time_start_vel,time_end_vel,time_start_temp,time_end_temp
type(dmumps_struc) idV                         ! Structure for MUMPS
type(dmumps_struc) idT                         ! Structure for MUMPS
logical MUMPS_T,MUMPS_VP                       ! Switches for MUMPS
real(8) fixt                                   ! USed in applying boundary conditions

real(8) F(4)                                   ! 
character(len=4) :: arg                        ! Used to read command line arguments
integer power,MUMPS_T_int,MUMPS_VP_int         ! Further command line arguments
real(8) Vsurf                                  ! Mean velocity at the surface

real(8) vr,er                                  ! Scaling factors
real(8), external :: rho_eos                   ! Equation of State function
real(8) rho_eos_el,drhodx_eos,drhody_eos       ! Equation of State values
logical use_eos                                ! Switch for using equation of state
real(8) temp_surf                              ! Surface temperature
real(8), dimension(4) :: temp_ref              ! Elemenetal reference temperature

!=====End declarations=====!


!Read from command line
!Whether to use MUMPS or LINPACK
!As solver
call get_command_argument(5,arg)
read(arg, '(i1)') MUMPS_T_int

call get_command_argument(6,arg)
read(arg, '(i1)') MUMPS_VP_int

MUMPS_T = MUMPS_T_int
MUMPS_VP = MUMPS_VP_int


!==============================================!
                                               !

if (MUMPS_VP .or. MUMPS_T) then
CALL mpi_init(ierr)                            !  
call mpi_comm_size (mpi_comm_world,nproc,ierr) !
call mpi_comm_rank (mpi_comm_world,iproc,ierr) !
                                               !
end if 
if(MUMPS_T) then
idT%COMM = MPI_COMM_WORLD                      ! Define a communicator for the package 
idT%SYM = 0                                    ! Ask for symmetric matrix storage 
idT%par=1                                      ! Host working 
idT%JOB = -1                                   ! Initialize an instance of the package 
call DMUMPS(idT)                               ! MUMPS initialisation
end if 

if(MUMPS_VP) then
idV%COMM = MPI_COMM_WORLD                      ! Define a communicator for the package 
idV%SYM = 1                                    ! Ask for symmetric matrix storage 
idV%par=1                                      ! Host working 
idV%JOB = -1                                   ! Initialize an instance of the package 
call DMUMPS(idV)                               ! MUMPS initialisation
end if


!Read approximation from command line
call get_command_argument(2,arg)
read(arg,'(i2)') approximation

!Chose which approximation to use
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
  buoyancy = .false.
case(3)
  !BA 
  BA = .true.
  buoyancy = .false.
case(4)
  !EBA 
  EBA = .true.
  buoyancy = .true.
end select

use_eos = .true.
!==============================================!
!=====[setup]==================================!
!==============================================!

chi_T=1.d-1
T0=1.d0

gamma = 1 

Lx=1.d0
Ly=1.d0


call get_command_argument(1,arg)

read(arg,'(i2)') nnx

nny=nnx

nelx=nnx-1
nely=nny-1

rho0=1.d0
hcapa=1.d0 
hcond=1.d0

call get_command_argument(3,arg)
read(arg, '(i1)') power

Ra = 10.0**power

gx=0


call get_command_argument(4,arg)
read(arg,'(f4.2)') Di
write(*,*) "Di = ", Di
gy=-10
alpha=1.d-1*Di
viscosity=Di/Ra

courant_nb=0.5

nstep=5000

output_freq=10

if (compressible) then
npicard=1000
else
npicard=1
end if
solver=2

if (compressible) then
temp_surf = 0.091*T0
else
temp_surf=0.d0
end if 

Ra = -alpha*gy*T0*hcapa*lx**3*rho0**2 / (hcond*viscosity)
Di = -alpha*gy*Lx / hcapa

chi_T = (alpha + alpha**2*T0*gamma)/(rho0*hcapa*gamma)

write(*,*) "Ra = ", Ra, "Di = ", Di
write(*,*) "Lx = ", Lx, "gy = ", gy
write(*,*) "Rho0 = ",Rho0, "T0 = ",T0
write(*,*) "Hcapa = ", hcapa, "Hcond = ", hcond
write(*,*) "alpha = ",alpha, "viscosity = ", viscosity
write(*,*) "chi_T = ",chi_T

vr = hcond/(rho0*hcapa*Lx)
er = viscosity*(hcond/(rho0*hcapa*Lx**2))**2*Lx**3


!==============================================!
!==============================================!

time=0.d0

nel=nelx*nely  ! total number of elements

np=nnx*nny     ! total number of nodes

NfemV=np*ndofV ! size of Stokes matrix
NfemT=np*ndofT ! size of temperature matrix

hx=Lx/(nnx-1)  ! grid spacing in x direction
hy=Ly/(nny-1)  ! grid spacing in y direction

!==============================================!
if(compressible) then
Cmat(1,1)=4.d0/3.d0  ; Cmat(1,2)=-2.d0/3.d0 ; Cmat(1,3)=0.d0  
Cmat(2,1)=-2.d0/3.d0 ; Cmat(2,2)=4.d0/3.d0  ; Cmat(2,3)=0.d0  
Cmat(3,1)=0.d0       ; Cmat(3,2)=0.d0       ; Cmat(3,3)=1.d0  
else
Cmat(1,1)=2.d0 ; Cmat(1,2)=0.d0 ; Cmat(1,3)=0.d0  
Cmat(2,1)=0.d0 ; Cmat(2,2)=2.d0 ; Cmat(2,3)=0.d0  
Cmat(3,1)=0.d0 ; Cmat(3,2)=0.d0 ; Cmat(3,3)=1.d0
end if 

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
open(unit=1009,file='OUT/temp.dat',status='replace')

open(unit=2000,file='OUT/times.dat',status='replace')
open(unit=1066,file='OUT/temp_eq.dat')
open(unit=1067,file='OUT/dTdx_nodal.dat')
open(unit=1068,file='OUT/dTdy_nodal.dat')
open(unit=1069,file='OUT/momentum.dat')

!==============================================!
!===[allocate memory]==========================!
!==============================================!

allocate(x(np))
allocate(y(np))
allocate(u(np))
allocate(v(np))
allocate(u_prev(np))
allocate(v_prev(np))
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
allocate(T_prev(np))
allocate(density(np))
allocate(qx(nel),qy(nel))
allocate(Tavrg(nny))
if(.not. MUMPS_VP) then
allocate(KKK(NfemV,NfemV))
allocate(KKKmem(NfemV,NfemV))
end if
allocate(rhs_f(NfemV))
allocate(rhs_h(nel))
allocate(G(NfemV,nel))
allocate(Vsol(NfemV))
allocate(Vsolmem(NfemV))
allocate(press_node(np))
allocate(Psolmem(Nel))
allocate(Rv(NfemV))
allocate(Rp(nel))

allocate(dTdx_nodal(np))
allocate(dTdy_nodal(np))
allocate(dTdx_elemental(nel))
allocate(dTdy_elemental(nel))

allocate(element_is_boundary(nel))
allocate(element_boundary_value(nel))

allocate(dxu_elemental(nel))
allocate(dxv_elemental(nel))
allocate(dyu_elemental(nel))
allocate(dyv_elemental(nel))
allocate(dxu_nodal(np))
allocate(dxv_nodal(np))
allocate(dyu_nodal(np))
allocate(dyv_nodal(np))
allocate(dxu_deviatoric(np))
allocate(dxv_deviatoric(np))
allocate(dyu_deviatoric(np))
allocate(dyv_deviatoric(np))
allocate(trace_strain_nodal(np))

allocate(phi_nodal(np))
allocate(phi_elemental(nel))
allocate(adiabatic_heating(nel))

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

!==============================================!
! MUMPS arrays for idT
!==============================================!

idT%N=NfemT                                    ! total number of degrees of freedom, size of FEM matrix
idT%NELT=nel                                   ! number of elements
LELTVAR=nel*m                                  ! nb of elts X size of elemental matrix
NA_ELT=nel*m*m                                 ! nb of elts X nb of nbs in elemental matrix !NEW

allocate(idT%A_ELT (NA_ELT))
allocate(idT%RHS   (idT%N))

if (iproc==0) then

allocate(idT%ELTPTR(idT%NELT+1))
allocate(idT%ELTVAR(LELTVAR))

do i=1,nel                                     !
   idT%ELTPTR(i)=1+(i-1)*m                     ! building ELTPTR array
end do                                         !
idT%ELTPTR(nel+1)=1+nel*m                      !

counter=0                                      !
do iel=1,nel                                   !
   do k=1,m                                    !
      inode=icon(k,iel)                        !
      counter=counter+1                        !
      idT%ELTVAR(counter)=inode                ! building ELTVAR
   end do                                      !
end do                                         !

end if ! iproc

idT%ICNTL(3) = 1020                            !
idT%ICNTL(4) = 2                               !




!==============================================!
! MUMPS arrays for idV
!==============================================!

idV%N=NfemV                                    ! total number of degrees of freedom, size of FEM matrix
idV%NELT=nel                                   ! number of elements
LELTVAR=nel*(m*ndofV)                          ! nb of elts X size of elemental matrix
NA_ELT=nel*(m*ndofV)*(m*ndofV+1)/2             ! nb of elts X nb of nbs in elemental matrix !NEW

allocate(idV%A_ELT (NA_ELT))
allocate(idV%RHS   (idV%N))

if (iproc==0) then

allocate(idV%ELTPTR(idV%NELT+1))
allocate(idV%ELTVAR(LELTVAR))

do i=1,nel                                     !
   idV%ELTPTR(i)=1+(i-1)*(ndofV*m)             ! building ELTPTR array
end do                                         !
idV%ELTPTR(i)=1+nel*(ndofV*m)                  !

counter=0                                      !
do iel=1,nel                                   !
   do k=1,m                                    !
      inode=icon(k,iel)                        !
      do idof=1,ndofV                          !
         iii=(inode-1)*ndofV+idof              !
         counter=counter+1                     !
         idV%ELTVAR(counter)=iii               ! building ELTVAR
      end do                                   !
   end do                                      !
end do                                         !

end if ! iproc

idV%ICNTL(3) = 6                               !
idV%ICNTL(4) = 1                               !

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
      end if 
   end do 
end do 

!==============================================!
!=====[Initial temperature field]==============!
!==============================================!

do i=1,np
   T(i)=T0*((1.d0-y(i)/Ly) - 0.01d0*cos(pi*x(i)/Lx)*sin(pi*y(i)/Ly))
end do

!==============================================!
!===========[Initial Pressure Field]===========!
!==============================================!

do i=1,nel
   press_el(i) = rho0*gy*(Lx-y(icon(1,i))+hx/2.d0) 
   !hydrostatic pressure to initialise
   !On the offchance anyone ever upgrades this to other geometries
   !Note that this does not generalise, and would need updating
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
V_diff_mem = 0.d0
V_diff_mem_mem = 0.d0

open(unit=666,file='conv_picard.dat',status='replace')

call cpu_time(time_start_vel)

do ipicard=1,npicard ! picard iterations 

!=========Set up density ======================!

!Get nodal pressure
call elemental_to_nodal(press_el,press_node,icon_inv,np,nel,m)

do i=1,np
  if (ALA .and. .not. use_eos) then
   density(i)=rho0*(1.d0-alpha*T(i)+chi_T*press_node(i))
  else if((BA) .or. (EBA) .or. (TALA) .and. .not. use_eos) then
   density(i)=rho0*(1.d0-alpha*T(i))
  else if (ALA .and.  use_eos) then
   density(i)=rho_eos(y(i),rho0,chi_T,gy,Ly)*(1.d0-alpha*T(i)+chi_T*press_node(i))
  else if((BA) .or. (EBA) .or. (TALA) .and.  use_eos) then
   density(i)=rho_eos(y(i),rho0,chi_T,gy,Ly)*(1.d0-alpha*T(i))
  end if  
end do

!==============================================!
!=====[build FE matrix]========================!
!==============================================!

u_prev=u
v_prev=v
if (.not. MUMPS_VP) then
KKK=0.d0
end if

G=0.d0
rhs_f=0.d0
rhs_h=0.d0

idT%RHS=0.d0
idT%A_ELT=0.d0
counter_mumps = 0

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

      rho_eos_el=0.d0
      drhodx_eos=0.d0
      drhody_eos=0.d0
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
         rho_eos_el = rho_eos_el + N(k)*rho_eos(y(icon(k,iel)),rho0,chi_T,gy,Ly)
         drhodx_eos=drhodx_eos + dNdx(k)*rho_eos(y(icon(k,iel)),rho0,chi_T,gy,Ly)
         drhody_eos=drhody_eos + dNdy(k)*rho_eos(y(icon(k,iel)),rho0,chi_T,gy,Ly)
      end do

      do i=1,m
         i1=2*i-1
         i2=2*i

         Bmat(1,i1)=dNdx(i)   ; Bmat(1,i2)=0.d0
         Bmat(2,i1)=0.d0      ; Bmat(2,i2)=dNdy(i)
         Bmat(3,i1)=dNdy(i)   ; Bmat(3,i2)=dNdx(i)     

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
      if (compressible) then
      !if (use_eos) then
      if (.true.) then
      hel(1)=hel(1)+1.d0/rho_eos_el*(uq*drhodx_eos+vq*drhody_eos)*wq*jcob
      else
      hel(1)=hel(1)+1.d0/rho2*(uq*drhodx2+vq*drhody2)*wq*jcob !approach 2
      end if 
      end if 

   end do
   end do

   if (MUMPS_VP) then
    !Impose BC
      do ii=1,m
        inode=icon(ii,iel)
         do k=1,ndofV
            ij=(inode-1)*ndofV+k
            if (bc_fixV(ij)) then
            fixt=bc_valV(ij)
            i=(ii-1)*ndofV+k
            Aref=Kel(i,i)
            do j=1,m*ndofV
               fel(j)=fel(j)-Kel(j,i)*fixt
               Kel(i,j)=0.d0
               Kel(j,i)=0.d0
            enddo
            Kel(i,i)=Aref
            fel(i)=Aref*fixt
            endif
         enddo
      enddo

    end if 




   !=====================
   !=====[assemble]======
   !=====================

   do k1=1,m
      ik=icon(k1,iel)
      do i1=1,ndofV
         ikk=ndofV*(k1-1)+i1
         m1= ndofV*(ik-1)+i1
         do k2=1,m
            jk=icon(k2,iel)
            do i2=1,ndofV
               jkk=ndofV*(k2-1)+i2
               if(MUMPS_VP) then
                  if (jkk>=ikk) then 
                     counter_mumps=counter_mumps+1   
                     idV%A_ELT(counter_mumps)=Kel(ikk,jkk) 
                  end if  
                else
               m2=ndofV*(jk-1)+i2
               KKK(m1,m2)=KKK(m1,m2)+Kel(ikk,jkk)
               end if
            end do
         end do
         if (MUMPS_VP) then
          idV%RHS(m1)=idV%RHS(m1)+fel(ikk)  
        else
         rhs_f(m1)=rhs_f(m1)+fel(ikk)
        end if 
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
      if (.not. MUMPS_VP) then
      Aref=KKK(i,i)
      do j=1,NfemV
         rhs_f(j)=rhs_f(j)-KKK(i,j)*bc_valV(i)
         KKK(i,j)=0.d0
         KKK(j,i)=0.d0
      enddo
      KKK(i,i)=Aref
      rhs_f(i)=Aref*bc_valV(i)
      end if
      !G(i,:)=0

      do j=1,nel
         rhs_h(j)=rhs_h(j)-G(i,j)*bc_valV(i)
         G(i,j)=0.d0
      end do

   endif
enddo

!==============================================!
!=====[solve system]===========================!
!==============================================!
if (.not.MUMPS_VP) then
KKKmem=KKK
end if

if (MUMPS_VP) then

call solve_uzawa_MUMPS(idV,G,NfemV,nel,rhs_h,Vsol,press_el)

else

call solve_uzawa2(KKK,G,NfemV,nel,rhs_f,rhs_h,Vsol,press_el)

end if 

if (MUMPS_VP) then
  do i=1,np
   u(i)=Vsol((i-1)*ndofV+1)
   v(i)=Vsol((i-1)*ndofV+2)
end do
else
do i=1,np
   u(i)=Vsol((i-1)*ndofV+1)
   v(i)=Vsol((i-1)*ndofV+2)
end do
end if

!==============================================!
!=====[iterations converged?]==================!
!==============================================!

V_diff=maxval(abs(Vsol-Vsolmem))/maxval(abs(Vsol))

P_diff=maxval(abs(press_el-Psolmem))/maxval(abs(press_el))

Vsolmem=Vsol
Psolmem=press_el

V_diff_mem_mem = V_diff_mem
V_diff_mem = V_diff


write(666,*) ipicard,V_diff,P_diff

write(*,'(a,i3,a,3es15.6,a,2es15.6)') &
'Picard it. ',ipicard,&
' | <|u|>,<|v|>,<|p|> ',sum(abs(u))/np,sum(abs(v))/np,sum(abs(press_el))/nel,&
' | V_diff, P_diff ',V_diff,P_diff

if (V_diff<tol .and. P_diff<tol) then
! print *,'-> Picard iterations have converged'
exit
else if (abs(V_diff_mem_mem - V_diff) .lt. eps) then
! print *,'-> Picard iterations have stalled'
exit
end if

end do 

close(666)

call cpu_time(time_end_vel)

!==============================================!
!========[Pressure Boundary Conditions]========!
!==============================================!

press_el=press_el-sum(press_el)/nel

!=========Set up density ======================!

!Get nodal pressure
call elemental_to_nodal(press_el,press_node,icon_inv,np,nel,m)

do i=1,np
  if (ALA .and. .not. use_eos) then
   density(i)=rho0*(1.d0-alpha*T(i)+chi_T*press_node(i))
  else if((BA) .or. (EBA) .or. (TALA) .and. .not. use_eos) then
   density(i)=rho0*(1.d0-alpha*T(i))
  else if (ALA .and.  use_eos) then
   density(i)=rho_eos(y(i),rho0,chi_T,gy,Ly)*(1.d0-alpha*T(i)+chi_T*press_node(i))
  else if((BA) .or. (EBA) .or. (TALA) .and.  use_eos) then
   density(i)=rho_eos(y(i),rho0,chi_T,gy,Ly)*(1.d0-alpha*T(i))
  end if  
end do


!==============================================!
!=======[Compute vrms and strainrates]=========!
!==============================================!

vrms=0.d0
dxu_elemental=0.d0
dxv_elemental=0.d0
dyu_elemental=0.d0
dyv_elemental=0.d0


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
         dNdx(k)=jcbi(1,1)*dNdr(k)+jcbi(1,2)*dNds(k)
         dNdy(k)=jcbi(2,1)*dNdr(k)+jcbi(2,2)*dNds(k)
         uq=uq+N(k)*u(icon(k,iel))
         vq=vq+N(k)*v(icon(k,iel))
         dxu_elemental(iel) = dxu_elemental(iel) + dNdx(k)*u(icon(k,iel))*wq*jcob/(hx*hy)
         dxv_elemental(iel) = dxv_elemental(iel) + dNdx(k)*v(icon(k,iel))*wq*jcob/(hx*hy)
         dyu_elemental(iel) = dyu_elemental(iel) + dNdy(k)*u(icon(k,iel))*wq*jcob/(hx*hy)
         dyv_elemental(iel) = dyv_elemental(iel) + dNdy(k)*v(icon(k,iel))*wq*jcob/(hx*hy)
      end do

      vrms=vrms+(uq**2+vq**2)*jcob*wq

   end do
   end do

end do
end do

open(unit=1070,file='OUT/dxu.dat')
open(unit=1071,file='OUT/dxv.dat')
open(unit=1072,file='OUT/dyu.dat')
open(unit=1073,file='OUT/dyv.dat')

iel=0
do j=1,nely
do i=1,nelx
   iel=iel+1
   write(1070,*) iel,dxu_elemental(iel)
   write(1071,*) iel,dxv_elemental(iel)
   write(1072,*) iel,dyu_elemental(iel)
   write(1073,*) iel,dyv_elemental(iel)
end do 
end do 

close(1070)
close(1071)
close(1072)
close(1073)

call elemental_to_nodal(dxu_elemental,dxu_nodal,icon_inv,np,nel,m)
call elemental_to_nodal(dxv_elemental,dxv_nodal,icon_inv,np,nel,m)
call elemental_to_nodal(dyu_elemental,dyu_nodal,icon_inv,np,nel,m)
call elemental_to_nodal(dyv_elemental,dyv_nodal,icon_inv,np,nel,m)

trace_strain_nodal = 1.d0/2.d0*(dxu_nodal + dyv_nodal)

dxu_deviatoric = dxu_nodal - trace_strain_nodal
dxv_deviatoric = dxv_nodal - trace_strain_nodal
dyu_deviatoric = dyu_nodal - trace_strain_nodal
dyv_deviatoric = dyv_nodal - trace_strain_nodal

vrms=sqrt(vrms/Lx/Ly)

write(1003,*) time,vrms ; call flush(1003)



open(unit=1074,file='OUT/velocities.dat')
do i=1,np
  write(1074,*) i, u(i), v(i)
end do

!===============================================!
!======[Calculate nodal viscous dissipation]====!
!===============================================!

phi_nodal=0.d0
if (compressible) then 
do i=1,np
   phi_nodal(i) =                4.d0/3.d0*viscosity*dxu_nodal(i)**2 - 2.d0/3.d0*dxu_nodal(i)*dyv_nodal(i)*viscosity
   phi_nodal(i) = phi_nodal(i) + 4.d0/3.d0*viscosity*dyv_nodal(i)**2 - 2.d0/3.d0*dxu_nodal(i)*dyv_nodal(i)*viscosity
   phi_nodal(i) = phi_nodal(i) + viscosity*(dyu_nodal(i) + dxv_nodal(i))**2
end do
else
do i=1,np
  phi_nodal(i) = 2.d0*viscosity*(dxu_nodal(i)**2 + dyv_nodal(i)**2) + viscosity*(dyu_nodal(i) + dxv_nodal(i))**2
end do
end if

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

call cpu_time(time_start_temp)

call elemental_to_nodal(press_el,press_node,icon_inv,np,nel,m)

do i=1,np
  if (ALA) then
   density(i)=rho0*(1.d0-alpha*T(i)+chi_T*press_node(i))
  else if((BA) .or. (EBA) .or. (TALA)) then
   density(i)=rho0*(1.d0-alpha*T(i))
  end if  
  !write(123,*) density(i)
end do
A_T=0.d0
B_T=0.d0

idT%RHS=0.d0
idT%A_ELT=0.d0
counter_mumps=0

do iel=1,nel

   AelT=0.d0
   BelT=0.d0

   temp(1)=T(icon(1,iel))
   temp(2)=T(icon(2,iel))
   temp(3)=T(icon(3,iel))
   temp(4)=T(icon(4,iel))

   temp_2(1)=T_prev(icon(1,iel))
   temp_2(2)=T_prev(icon(2,iel))
   temp_2(3)=T_prev(icon(3,iel))
   temp_2(4)=T_prev(icon(4,iel))

   temp_ref(1) = temp_surf*exp(-alpha*gy*(Ly-y(icon(1,iel)))/hcapa)
   temp_ref(2) = temp_surf*exp(-alpha*gy*(Ly-y(icon(2,iel)))/hcapa)
   temp_ref(3) = temp_surf*exp(-alpha*gy*(Ly-y(icon(3,iel)))/hcapa)
   temp_ref(4) = temp_surf*exp(-alpha*gy*(Ly-y(icon(4,iel)))/hcapa)

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
      str_inv2=0.d0
      e_kk = 0.d0
      rho2 = 0.d0
      do k=1,m
         xq=xq+N(k)*x(icon(k,iel))
         yq=yq+N(k)*y(icon(k,iel))
         uq=uq+N(k)*u(icon(k,iel))
         vq=vq+N(k)*v(icon(k,iel))
         dNdx(k)=jcbi(1,1)*dNdr(k)+jcbi(1,2)*dNds(k)
         dNdy(k)=jcbi(2,1)*dNdr(k)+jcbi(2,2)*dNds(k)

         phi = phi + phi_nodal(icon(k,iel))*N(k)*wq*jcob!/(hx*hy)

         rho_el = rho_el + density(icon(k,iel))*N(k)
         v_el = v_el + v(icon(k,iel))*N(k)
         !if (use_eos) then
         if (.false.) then
         rho2 = rho2 + N(k)*rho_eos(y(icon(k,iel)),rho0,chi_T,gy,Ly)
         else
         rho2 = rho2 + N(k)*density(icon(k,iel))
         end if
      end do
      if (BA) then
        F=0.d0
      else
        F= N*phi
      end if 
      velnorm=sqrt(uq**2+vq**2)

      vel2D(1,1)=uq
      vel2D(1,2)=vq

      Nvect(1,:)=N(:)
      NvectT(:,1)=N(:)

      BmatT(1,:)=dNdx
      BmatT(2,:)=dNdy

      BmatTT=transpose(BmatT)

      if (BA) then 
        Kb=0.d0
      else
        Kb=-matmul(NvectT,Nvect)*wq*jcob*alpha*gy*v_el*rho2
      end if 

      Ka=matmul(NvectT,matmul(vel2D,BmatT))*wq*jcob*rho2*hcapa

      Kc=matmul(BmatTT,BmatT)*wq*jcob*hcond

      KK=(Ka+Kc+Kb)

      M_T=matmul(NvectT,Nvect)*wq*jcob*rho2*hcapa

      AelT=AelT+(M_T+KK*theta*dt)

      BelT=BelT+matmul(M_T-KK*(1.d0-theta)*dt,temp(1:4)) + F*dt 
      if (ALA .or. TALA) then 
        Belt = Belt - matmul(Kc,temp_ref(1:4))*dt !Represents convection of reference temp
      end if

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
    if(MUMPS_T) then
      do k1=1,m   
      m1=icon(k1,iel) 
      do k2=1,m   
         !if (k2>=k1) then 
         counter_mumps=counter_mumps+1   
         idT%A_ELT(counter_mumps)=AelT(k2,k1) 
         !end if
      end do    
      idT%RHS(m1)=idT%RHS(m1)+BelT(k1)    
      end do  
        
    else

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
    end if

end do ! end of loop over cells for building temp matrices

!Before solving system, but after constructing matrices
!Store old temperature 

T_prev = T

!==============================================!
!=====[solve system]===========================!
!==============================================!

if(MUMPS_T) then

idT%ICNTL(5) = 1                               ! elemental format

idT%JOB = 1                                    ! analysis phase 
CALL DMUMPS(idT)

idT%JOB = 2                                    ! factorisation phase 
CALL DMUMPS(idT)

idT%JOB = 3                                    ! solve phase
CALL DMUMPS(idT)

T=idT%RHS                                      ! Transfer solution
else

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

end if
write(*,*) 'min/max T',minval(T),maxval(T), sum(T)/np

open(unit=123,file='OUT/solution_T.dat',status='replace')
do i=1,np
   write(123,'(5f20.10)') x(i),y(i),u(i),T(i),T(i)-(1-y(i))
end do
close(123)
call cpu_time(time_end_temp)

write(2000,*) time_end_vel-time_start_vel, time_end_temp-time_start_temp ; call flush(2000)

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
write(1009,*) time,sum(T)/np ; call flush(1009)


vsurf=0.d0
do i=1,np
  if (y(i) .lt. eps) then
    vsurf = vsurf + abs(u(i))
  end if 
end do
vsurf = vsurf / nnx

write(*,*) 'Ra=',Ra,'Di=',Di
write(*,*) 'Nu=',Nu,'Vrms=',vrms
write(*,*) "Vsurf = ", vsurf
write(*,*) 'Nu_rel=',Nu/T0,'Vrms_rel=',vrms/vr 

!==============================================!
!=====[output data in vtu format]==============!
!==============================================!

if (istep==1 .or. mod(istep,output_freq)==0) &
call output_for_paraview (np,nel,x,y,u,v,press_el,T,density,qx,qy,icon,istep,dxu_nodal,dyv_nodal,dxv_nodal,dyu_nodal,dxu_elemental,dyv_elemental,dxv_elemental,dyu_elemental,phi_nodal,phi_elemental,adiabatic_heating)

!==============================================!
!====[is steady state reached ?]===============!
!==============================================!

chi_Nu=abs(Nu-Nu_prev)/Nu
chi_vrms=abs(vrms-vrms_prev)/vrms 

write(1006,*) time,chi_Nu   ; call flush(1006)
write(1007,*) time,chi_vrms ; call flush(1007)

 write(*,*) 'chi_Nu  =',chi_Nu,'tol=',tol_nu
 write(*,*) 'chi_vrms=',chi_vrms,'tol=',tol_nu

if ((chi_Nu<tol_nu .and. chi_vrms<tol_nu) .or. (istep .eq. nstep)) then


  ! if (.not.(istep==1 .or. mod(istep,output_freq)==0)) &
  ! call output_for_paraview (np,nel,x,y,u,v,press_el,T,density,qx,qy,icon,istep,dxu_nodal,dyv_nodal,dxv_nodal,dyu_nodal,dxu_elemental,dyv_elemental,dxv_elemental,dyu_elemental,phi_nodal,phi_elemental,adiabatic_heating)


   ! output avrg temperature profile and temp array
   open(unit=9999,file='OUT/temp_final.dat')
   counter=0
   do j=1,nny
   do i=1,nnx
      counter=counter+1
      Tavrg(j)=Tavrg(j)+T(counter)
      write(9999,"(f10.5)",advance='no') T(counter)
   end do
   write(9999,*) new_line("x")
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

