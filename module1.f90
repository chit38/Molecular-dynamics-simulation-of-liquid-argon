module md
use variables
implicit none

contains
!----------------------------------------------------------------------------------------------------------
subroutine create_lattice()
    implicit none

    integer :: index, i, j, k
    real*8 :: dx
    

    print *,"Please enter number of atoms"
    read(*,*)natoms
    print *,"Please enter density"
    read(*,*)rho
    print *,"Please enter Temperature"
    read(*,*)temp
    print *,"Please enter step size and number of steps"
    read(*,*)dt, maxstep

    t0 = temp
    L = (dfloat(natoms)*2.d0/rho)**(1.d0/3.d0)
    dx = L/12.d0
    print *,L, dx
    allocate(pos(ndim,natoms))
    index = 0
    do i = 1, 12
        do j = 1, 12
            do k = 1,6
                index = (i-1)*72 + (j-1)*6 + k
                pos(1, index) = 0.5 + (i-1)*dx
                pos(2, index) = 0.5 + (j-1)*dx
                pos(3, index) = 0.5 + (k-1)*dx
            end do
        end do
    end do
end subroutine create_lattice
!---------------------------------------------------------------------------------

SUBROUTINE get_initial_velocities
    IMPLICIT NONE
    REAL*8, PARAMETER :: pi=4.d0*ATAN(1.d0)
    REAL*8 :: dum(2), sigma
    INTEGER :: ii, i, k

    allocate(vel(3,natoms))
    sigma=DSQRT(temp)
    ii=0
    DO i=1,natoms
      DO k=1,ndim
        CALL RANDOM_NUMBER(dum(1:2))
        ii=ii+1
        IF(MOD(ii,2)==0)THEN
          vel(k,i)=SIGMA*DSQRT(-2.D0*DLOG(dum(2)))*DCOS(2.D0*pi*dum(1))
        ELSE
          vel(k,i)=SIGMA*DSQRT(-2.D0*DLOG(dum(2)))*DSIN(2.D0*pi*dum(1))
        END IF
      END DO
    END DO
    !CALL temperature(temp)
    CALL rescale_velocities(temp)
    WRITE(6,*) "Initial velocities are assigned for T=", TEMP
  END SUBROUTINE get_initial_velocities
!---------------------------------------------------------------------------------
SUBROUTINE rescale_velocities(t)
    IMPLICIT NONE

    REAL*8,INTENT(IN)  :: t

    REAL*8  :: t0,scal

    CALL temperature(t0)
    print *, t0
    scal=DSQRT(t/t0)
    vel(:,:)=vel(:,:)*scal

  END SUBROUTINE
!-----------------------------------------------------------------------------------------
SUBROUTINE temperature(t)
    IMPLICIT NONE
    REAL*8, intent(out) :: t

    INTEGER :: i
    REAL*8 :: mv2
  
    mv2=0.D0
    DO i=1,natoms
      mv2=mv2+DOT_PRODUCT(vel(1:ndim,i),vel(1:ndim,i)) 
    END DO
    t=(1.d0/3.d0)*mv2/DFLOAT(natoms)
  END SUBROUTINE temperature
!-----------------------------------------------------------------------------------------
  SUBROUTINE velocity_verlet_r
    IMPLICIT NONE
 
    INTEGER :: i
!
    DO i=1,natoms
      pos(1:ndim,i)=pos(1:ndim,i)+vel(1:ndim,i)*dt
    END DO
  END SUBROUTINE velocity_verlet_r
!-----------------------------------------------------------------------------------------------
  SUBROUTINE velocity_verlet_v
    IMPLICIT NONE

    INTEGER :: i
!
    DO i=1,natoms
      vel(1:ndim,i)=vel(1:ndim,i)+0.5d0*dt*force(1:ndim,i)
    END DO
  END SUBROUTINE velocity_verlet_v
!-----------------------------------------------------------------------------
SUBROUTINE print_geo(filen,iacc)
    IMPLICIT NONE

    CHARACTER (LEN=*), INTENT(IN)   :: filen
    INTEGER, INTENT(IN)             :: iacc

    INTEGER :: i,j

    IF(iacc.EQ.0)THEN
      OPEN(1,FILE=FILEN,STATUS='UNKNOWN',FORM='FORMATTED')
    ELSE
      OPEN(1,FILE=FILEN,STATUS='UNKNOWN',FORM='FORMATTED',ACCESS='APPEND')
    END IF

    !IF(ndim/=2)STOP 'print_geo() is not implemented for ndim/=2'
    WRITE(1,'(I10)')natoms
    WRITE(1,*)"comment"
    DO i=1,natoms
      DO j=1,ndim
        IF(pos(j,i)>L)pos(j,i)=pos(J,I)-L
        IF(pos(j,i)<0.d0)pos(j,i)=pos(j,i)+L
      END DO
      !if(pos(3,i)>L/2)then 
      !pos(3,i) = pos(3,i)-L/2
      !end if
      WRITE(1,"(A,3F16.6)") "Ar", pos(1:3,i)
    END DO
    CLOSE(1)
  END SUBROUTINE print_geo
!----------------------------------------------------------------------------------
  SUBROUTINE energy_forces
    IMPLICIT NONE
    INTEGER  :: i, j
    REAL*8   :: dd(ndim), d, d2, d6, d12, fjk(ndim), fik(ndim)
!
    
    force(:,:)=0.d0
    energy=0.D0

!   Loop of all the atom pairs (i,j)
    DO i=1,natoms-1
      DO j=i+1,natoms

        dd(1:ndim)=pos(1:ndim,i)-pos(1:ndim,j)
        dd(1:ndim)=dd(1:ndim)-L*NINT(dd(1:ndim)/L)

!       Calculate pair distance
        d=DSQRT(DOT_PRODUCT(dd(1:ndim),dd(1:ndim)))

       IF(d<1.0e-4 )THEN  ! if the distance is too small, STOP
         PRINT *, "i=", i, " j=", j, " d=", d
         PRINT *, "Atom i = ", pos(1:3,i)
         PRINT *, "Atom j = ", pos(1:3,j)
         STOP 'ERROR! Atoms i and j are overlapping!'
       END IF

!       Calculate Energy
        d2=1.d0/(d**2)
        d6=d2*d2*d2
        d12=d6*d6
        energy=energy+d12-d6

!       Calculate Force
        fik(1:ndim)=-(-2.D0*d12+d6)*dd(1:ndim)*d2
        fjk(1:ndim)=-fik(1:ndim)

        force(1:ndim,i)=force(1:ndim,i)+fik(1:ndim)
        force(1:ndim,j)=force(1:ndim,j)+fjk(1:ndim)

      END DO
    END DO
    energy=energy*4.D0
    !print *, "energy =", energy
    force(:,:)=force(:,:)*24.d0
  END SUBROUTINE energy_forces
!---------------------------------------------------------------------------------
  SUBROUTINE print_energy(istep)
    IMPLICIT NONE
    INTEGER :: istep
    REAL*8  :: KE, TE
    INTEGER :: ICALL=0
    SAVE    :: ICALL
!
    icall=icall+1 
    IF(ICALL.EQ.1)THEN
      OPEN(11,FILE='energy.dat',STATUS='UNKNOWN')
      WRITE(*, '(5A12)')'ISTEP', 'TEMP.', 'K.E.', 'P.E.', 'T.E'
    END IF
!
    ke=1.5d0*dfloat(natoms)*temp
    te=ke+energy
 !
    WRITE(*,"(I12, 4F24.4)")istep, temp, ke, energy, te
    WRITE(11,*)istep, temp, ke, energy, te
END SUBROUTINE print_energy
!------------------------------------------------------------------------------

end module md