module variables
implicit none
integer, parameter :: ndim = 3
integer :: natoms, maxstep
real*8 :: rho, L, dt, temp, energy
real*8, allocatable ::pos(:,:), vel(:,:), force(:,:)
REAL*8, PARAMETER :: pi=4.d0*ATAN(1.d0)
end module variables