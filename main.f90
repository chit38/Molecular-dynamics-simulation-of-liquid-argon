program main
use variables
use md
use gr
implicit none

integer :: istep, T= 0

call create_lattice()
call get_initial_velocities
CALL print_geo('traj.xyz',0)
allocate(force(ndim,natoms))
CALL energy_forces
CALL print_energy(0)
DO istep=1,maxstep
!   Velocity Verlet Integrator
    CALL velocity_verlet_v
    CALL velocity_verlet_r
    CALL energy_forces
    CALL velocity_verlet_v
!   Properties: Temperature 
    CALL temperature(temp)
    IF(MOD(istep,50).EQ.0)THEN 
     CALL print_geo('traj.xyz',istep)
     CALL print_energy(istep)
     T = T + 1
    END IF
  END DO

call calculate_gr(T-1)
end program main