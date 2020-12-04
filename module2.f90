module gr
use variables
implicit none

contains
subroutine calculate_gr(T)
integer, intent(in) :: T
integer :: maxiter, i, j, k, y, bin, nat
real*8 :: r = 0, d(3), rjk
real*8, allocatable :: dum(:), g(:), xyz(:,:,:)
character(len=2) :: comment, at

maxiter = int((L/2)/0.1d0)+1
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
       read(15,'(a2,3f15.6)') at, xyz(i,j,1:3)
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
        bin = int(rjk/0.1d0) + 1
        if (bin > maxiter) exit
        dum(bin) = dum(bin) + 1.d0
        !write(11,*) bin
        endif
      end do
   end do
end do


open(3, file="hist.txt", status="unknown", form='formatted')
do bin = 1, maxiter-1
    r = DFLOAT(bin)*0.1d0
    g(bin) = real(dum(bin))/(2.d0*rho*r*r*0.1d0*pi*T*natoms)
    write(3,*)r, g(bin)
end do
close(15)
end subroutine
end module gr