module gr
use variables
implicit none

contains
subroutine calculate_gr(T)
integer, intent(in) :: T
integer :: maxiter, i, j, k, x, y, z, bin,nij = 0, nat
real*8 :: r = 0, d(3), rjk, r_lower, r_upper,zero
real*8, allocatable :: dum(:), g(:), xyz(:,:,:)
character(len=2) :: comment, at

maxiter = int((L/2)/0.1d0)
!print *, maxiter
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
do y = 1, T
   do j = 1, natoms-1
      do k = j+1,natoms
            d(1:3) = xyz(y, k, 1:3) - xyz(y, j, 1:3)
            d(1:3) = d(1:3) - L*nint(d(1:3)/L)
            rjk = dot_product(d(1:3), d(1:3))
            rjk = sqrt(rjk)
            if(rjk < 0.5d0*L) then
            bin = int(rjk/0.1d0) + 1
            if (bin > maxiter) exit
            dum(bin) = dum(bin) + 1.d0
            write(11,*) bin
            endif
      end do
   end do
end do


open(3, file="hist.txt", status="unknown", form='formatted')
do bin = 1, maxiter-1
    r = DFLOAT(bin)*0.1d0
    g(bin) = dum(bin)/(pi*r*0.1d0*63*T*rho)
    write(3,*)r, g(bin)
end do

end subroutine
end module gr