!	evolve all ptcicles (only free streaming at this moment)
	subroutine shift_particle_single_evt(N_ptc, mean_x, mean_y)

	implicit none

	include "parameters_common_blocks.inc"
			
	double precision x, y, z

	integer i_ptc
	integer N_ptc

	double precision mean_x, mean_y, mean_z

!	calculate mean_x, mean_y
	mean_x = 0.d0
	mean_y = 0.d0
	mean_z = 0.d0

	do i_ptc = 1, N_ptc

		x = x_ptc(i_ptc)
		y = y_ptc(i_ptc)
		z = z_ptc(i_ptc)

		mean_x = mean_x + x
		mean_y = mean_y + y
		mean_z = mean_z + z 

	end do

	mean_x = mean_x / N_ptc
	mean_y = mean_y / N_ptc
	mean_z = mean_z / N_ptc

!	shift to center of mass frame
	do i_ptc = 1, N_ptc

		x = x_ptc(i_ptc)
		y = y_ptc(i_ptc)
		z = z_ptc(i_ptc)

		x_ptc(i_ptc) = x - mean_x
		y_ptc(i_ptc) = y - mean_y
!               do not shift z axis
!                z_ptc(i_ptc) = z - mean_z
	end do

	end subroutine
