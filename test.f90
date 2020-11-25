program test
implicit none
integer :: i
real :: a

open(1,file = "test.dat", status = "old")
do i = 1, 1000
call random_number(a)
write(1,*)a
end do
end program test