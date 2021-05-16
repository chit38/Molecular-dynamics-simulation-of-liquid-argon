module rmsdisplacement
use variables
implicit none

contains
subroutine rmsd(T)
integer, intent(in) :: T
real*8, allocatable :: xyz(:,:,:)
real*8 :: rdum(3), rj, r2=0
integer :: i, j, k, nat, t0_max
character(len=2):: comment, at

t0_max = T/2
allocate(xyz(T,natoms,3))
open(13,file="traj.xyz",status='unknown')
 do i = 1,T
    read(13,*)nat
    read(13,*)comment
    do j = 1,nat
       read(13,*) at, xyz(i,j,1:3)
    enddo
enddo
    
open(17, file = "r2.txt", status = 'unknown', form = 'formatted')
    do i = 1, t0_max-1
        if(mod(i,50).eq.0)then
            do j = 1, natoms
                do k = 1, t0_max
                    rdum(1:3) = xyz(i, j, 1:3) - xyz(i+k, j, 1:3)
                    rdum(1:3) = rdum(1:3) - L*nint(rdum(1:3)/L)
                    rj = sqrt(rdum(1)**2 + rdum(2)**2 + rdum(3)**2)
                    r2 = r2 + rj
                end do
            end do
            write(17,*)i,r2/(natoms*t0_max)
        end if
    end do
    
end subroutine rmsd


end module rmsdisplacement