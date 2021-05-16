module gr
use variables
implicit none

contains
subroutine calculate_gr(T)
integer, intent(in) :: T
integer :: maxiter, i, j, k, y, bin, nat
real*8 :: r = 0, d(3), rjk, step
real*8, allocatable :: dum(:), g(:), xyz(:,:,:)
character(len=2) :: comment, at

step = 0.05d0
maxiter = int((L/2.d0)/step)
print *, maxiter
allocate(dum(maxiter))
allocate(g(maxiter))
dum(1:maxiter) = 0.d0
allocate(xyz(T, natoms, 3))
open(15,file="traj.xyz",status='unknown')
 do i = 1,T
    read(15,*)nat
    read(15,*)comment
    do j = 1,nat
       read(15,*) at, xyz(i,j,1:3)
    enddo
enddo

!do z = 1, maxiter
do i = 1, T
do j = 1, natoms-1
    do k = j+1,natoms
        d(1:3) = xyz(i, k, 1:3) - xyz(i, j, 1:3)
        d(1:3) = d(1:3) - L*nint(d(1:3)/L)
        rjk = sqrt(d(1)**2 + d(2)**2 + d(3)**2)
        if(rjk <0.5*L) then
        bin = int(rjk/step)
        if (bin > maxiter) exit
        dum(bin) = dum(bin) + 1.d0
        endif
      end do
   end do
end do


open(3, file="hist.txt", status="unknown", form='formatted')
do bin = 1, maxiter-1
    r = DFLOAT(bin)*step
    g(bin) = dum(bin)/(2.d0*rho*r*r*step*pi*T*natoms)
    write(3,*)r, g(bin)
end do
close(3)
close(15)
end subroutine
end module gr