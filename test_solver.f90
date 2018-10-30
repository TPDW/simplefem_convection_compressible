program test_solver
implicit none
real(8), dimension(25,25) :: A 
real(8), dimension(25) :: b,x
integer i,np,j,job
real(8), dimension(:), allocatable :: work,ipvt
real(8) rcond 
external solve_linpack
np=size(A,1)
A=0
do i=1,25
A(i,i)=2.d0
B(i)=1.d0
do j=1,25
if (i .gt. j) A(i,j)=1.d0
end do
end do
write(*,*) size(A,1)
call solve_linpack(A,b,np,x)


! job=0
! allocate(work(Np))
! allocate(ipvt(Np))
! call DGECO (A, np, np, ipvt, rcond, work)
! write(*,*) "rcond"
! write(*,*) rcond

! call DGESL (A, np, np, ipvt, B, job)
! deallocate(ipvt)
! deallocate(work)

do i=1,25
write(*,*) B(i)
end do 
write(*,*) "Batman!"
end program