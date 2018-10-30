subroutine solve_uzawa_MUMPS(idV,G,Nfem,nel,rhs_h,Vsol,Psol)
implicit none

include 'mpif.h'                               !   
include 'dmumps_struc.h'      

integer, intent(in) :: Nfem, nel
real(8), dimension(Nfem,Nfem) :: KKK,KKKmem
real(8), dimension(Nfem,Nfem) :: G(Nfem,nel)
real(8), dimension(Nfem) :: Vsol,rhs_f,Vsolmem
real(8), dimension(nel) :: Psol,Psolmem
real(8), dimension(nel) :: rhs_h,q

real(8), dimension(Nfem) :: B,Rv,h,phi
!real(8), dimension(Nfem,1) :: GP

! real(8), dimension(:),   allocatable :: work   ! work array needed by the solver
! integer, dimension(:), allocatable :: ipvt     ! work array needed by the solver 
integer iter,job
real(8) rcond,P_diff,V_diff

real(8) :: alpha
real(8), parameter :: tol=1.d-6
integer, parameter :: niter=250
type(dmumps_struc) idV,idV_const

! open(unit=123,file='conv_uzawa2.dat')

Psol=0.d0
Vsol=0.d0
Psolmem=0.d0
Vsolmem=0.d0
   
idV_const = idV
!write(*,*) maxval(idV%A_elt), maxval(idV_const%A_elt)
! allocate(work(Nfem))
! allocate(ipvt(Nfem))
! call DGECO (KKK, Nfem, Nfem, ipvt, rcond, work)
!write(*,*) maxval(idV%rhs), minval(idV%rhs)
!compute u1
!B=rhs_f-matmul(G,Psol)
idV%rhs = idV%rhs - matmul(G,Psol)
!write(*,*) maxval(idV%rhs), minval(idV%rhs)

! call DGESL (KKK, Nfem, Nfem, ipvt, B, job) 
idV%ICNTL(5) = 1                               ! elemental format

idV%JOB = 1                                    ! analysis phase 
CALL DMUMPS(idV)

idV%JOB = 2                                    ! factorisation phase 
CALL DMUMPS(idV)

idV%JOB = 3                                    ! solve phase
CALL DMUMPS(idV)

Vsol=idV%rhs
! write(*,*) "Start Loop"
!write(*,*) maxval(Vsol)
do iter=1,niter
   idV%A_elt = idV_const%A_elt
   !compute qk

   q=rhs_h-matmul(transpose(G),Vsol)

   !compute pk

   phi=matmul(G,q) 

   !compute hk

   idV%rhs=phi

   ! call DGESL (KKK, Nfem, Nfem, ipvt, B, job) 
   idV%ICNTL(5) = 1                               ! elemental format
   ! write(*,*) 1
   idV%JOB = 1                                    ! analysis phase 
   CALL DMUMPS(idV)

   idV%JOB = 2                                    ! factorisation phase 
   CALL DMUMPS(idV)

   idV%JOB = 3                                    ! solve phase
   CALL DMUMPS(idV)

   h=idV%rhs
   !write(*,*) 
   !compute alpha

   alpha=dot_product(q,q)/dot_product(phi,h) !; print *,alpha

   !update pressure

   Psol=Psol-alpha*q

   !update velocity

   Vsol=Vsol+alpha*h

   !check for convergence

   !Rv=matmul(KKKmem,Vsol)+matmul(G,Psol)-rhs_f

   V_diff=maxval(abs(Vsol-Vsolmem))/maxval(abs(Vsol))
   P_diff=maxval(abs(Psol-Psolmem))/maxval(abs(Psol))

    !write(*,*) iter,V_diff,P_diff,maxval(abs(Vsol-Vsolmem)),maxval(abs(Psol-Psolmem))

   Psolmem=Psol
   Vsolmem=Vsol

   if (max(V_diff,P_diff)<tol) exit

end do
! write(*,*) "End loops"
! deallocate(ipvt)
! deallocate(work)

! close(123)

end subroutine



