!	evolve all particles (only free streaming at this moment)
	subroutine evolve_particle_single_evt_glauber_taueta(N_ptc)

	implicit none

	include "parameters_common_blocks.inc"
			
	double precision x, y, z
	double precision px, py, pz, pxyz
	double precision vx, vy, vz
	double precision dt
	double precision rpd, eta
        integer ieta
       
	integer i_ptc
	integer N_ptc
 
!       here we use free-streaming approximation
!       so each particle contribute to a certain eta no matter what the evolution time is (eta = rpd)
!       For each particle, for a fixed tau, we have different evolution time since we have fixed tau & eta
!	dt = t_fs

!	evolve all particles
	do i_ptc = 1, N_ptc

		x = x_ptc(i_ptc)
		y = y_ptc(i_ptc)
		z = z_ptc(i_ptc)

		px = px_ptc(i_ptc)
		py = py_ptc(i_ptc)
		pz = pz_ptc(i_ptc)
		rpd = rpd_ptc(i_ptc)

		pxyz = sqrt(px*px + py*py + pz*pz)
!		massless particles
		vx = px/pxyz
		vy = py/pxyz
		vz = pz/pxyz

!		this assumes z = 0 when t=0
!               for free-streaming, eta = rpd
!		rpd = 0.5d0 * log((pxyz+pz)/(pxyz-pz))
                eta = rpd
                dt = tau0_hydro * cosh(eta)
!                dt = tau0_hydro / sqrt(1.d0 - vz**2)

		x_ptc(i_ptc) = x + vx * dt
		y_ptc(i_ptc) = y + vy * dt
		z_ptc(i_ptc) = z + vz * dt

!		only free streaming, no momentum change
!		px_ptc(i_ptc) = px
!		py_ptc(i_ptc) = py
!		pz_ptc(i_ptc) = pz
!		rpd_ptc(i_ptc) = rpd

		x = x_ptc(i_ptc)
		y = y_ptc(i_ptc)
		z = z_ptc(i_ptc)

!		now calculate space-time rapidity
		eta = 0.5d0 * log((dt+z)/(dt-z)) 
                eta_ptc(i_ptc) = eta

!		write (*, *) i_ptc, px, py, pz, dt, z, rpd, eta
!		write (12, *) i_ptc, eta, rpd

                ieta = nint(eta/deta)
                nevt_eta(ieta) = nevt_eta(ieta) + 1.d0
                nevt_eta_tot = nevt_eta_tot + 1.d0
                if (abs(eta) .lt. 0.5d0) nevt_eta_unit = nevt_eta_unit + 1.d0

	end do

	end subroutine
