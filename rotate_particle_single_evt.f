!	rotate the distribution (only free streaming at this moment)
	subroutine rotate_particle_single_evt(N_ptc, mean_x, mean_y, eps_n, phi_n, v_n)

	implicit none

	include "parameters_common_blocks.inc"
			
	double precision x, y, z
	double precision px, py, pz

	integer i_ptc, i_n1, i_n2
	integer N_ptc

	double precision mean_x, mean_y, mean_z
	double precision mean_px, mean_py, mean_pz

	integer i_mm
	double precision eps_n(10) 
	double precision phi_n(10) 
	double precision v_n(10) 

	integer i_rt

	double precision rt, phi, phi_xy
	double precision mean_rtn(10), mean_rtn_cosn(10), mean_rtn_sinn(10)

	double precision pt, phi_p, phi_pxpy
	double precision mean_cosn_phip_phin(10)
	double precision mean_ptn(10), mean_ptn_cosn(10), mean_ptn_sinn(10)

	double precision phi2
	double precision x2, y2, z2

!	the event plane angle for 2nd moments
	phi2 = phi_n(2) + PI/2.d0

	do i_ptc = 1, N_ptc

		x = x_ptc(i_ptc)
		y = y_ptc(i_ptc)
		z = z_ptc(i_ptc)

		x2 = x * cos(phi2) + y * sin(phi2)
		y2 = -x * sin(phi2) + y * cos(phi2)
		z2 = z

		x_ptc(i_ptc) = x2
		y_ptc(i_ptc) = y2
		z_ptc(i_ptc) = z2

	end do

!       rotate positions (x,y) of original binary collisions
        do i_n1 = 1, atom_A
            do i_n2 = 1, atom_B

                x = x_coll(i_n1, i_n2)
                y = y_coll(i_n1, i_n2)

                x2 = x * cos(phi2) + y * sin(phi2)
                y2 = -x * sin(phi2) + y * cos(phi2)

                x_coll(i_n1, i_n2) = x2
                y_coll(i_n1, i_n2) = y2

            end do
        end do

	return


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
                z_ptc(i_ptc) = z - mean_z
	end do

!	calculate spatial moments in center of mass frame
	do i_mm = 1, 10
		phi_n(i_mm) = 0.d0
	end do

	do i_mm = 1, 10
		mean_rtn(i_mm) = 0.d0
		mean_rtn_cosn(i_mm) = 0.d0
		mean_rtn_sinn(i_mm) = 0.d0
	end do

	do i_ptc = 1, N_ptc

		x = x_ptc(i_ptc)
		y = y_ptc(i_ptc)
		z = z_ptc(i_ptc)

		rt = sqrt(x*x + y*y)
		phi_xy = atan2(y,x)

		px = px_ptc(i_ptc)
		py = py_ptc(i_ptc)
		pz = pz_ptc(i_ptc)

		pt = sqrt(px*px + py*py)
		phi_pxpy = atan2(py, px)

		do i_mm = 1, 10

			i_rt = i_mm
			if (i_mm .eq. 1) i_rt = 3

			mean_rtn(i_mm) = mean_rtn(i_mm) + rt**i_rt * pt**index_pt
			phi = phi_xy - phi_n(i_mm)
			mean_rtn_cosn(i_mm) = mean_rtn_cosn(i_rt) + rt**i_mm * cos(i_mm*phi) * pt**index_pt
			mean_rtn_sinn(i_mm) = mean_rtn_sinn(i_rt) + rt**i_mm * sin(i_mm*phi) * pt**index_pt
		end do

	end do

	do i_mm = 1, 10
		mean_rtn(i_mm) = mean_rtn(i_mm) / N_ptc
		mean_rtn_cosn(i_mm) = mean_rtn_cosn(i_mm) / N_ptc
		mean_rtn_sinn(i_mm) = mean_rtn_sinn(i_mm) / N_ptc
	end do

!	if only one particle produced, then perfect symmetry, no higher spatial moments: rt^n=0 in shifted frame
	if (N_ptc .gt. 1) then

		do i_mm = 1, 10

			if (mean_rtn(i_mm) .eq. 0.d0) then
!				make sure N_ptc particles do not come from the same place
				if (x_ptc(1) .ne. x_ptc(2)) then
					if (i_mm .eq. N_ptc) then
						eps_n(i_mm) = 1.d0
!						phi_n(i_mm) = 2.d0*PI/i_mm
						write (*,*) "perfect spatial symmetry for", N_ptc, "particles." 
						write (*,*) N_ptc, mean_rtn_cosn(i_mm), mean_rtn_sinn(i_mm), mean_rtn(i_mm), eps_n(i_mm), phi_n(i_mm)
					end if
				end if
				cycle
			end if

			eps_n(i_mm) = sqrt(mean_rtn_cosn(i_mm)**2 + mean_rtn_sinn(i_mm)**2) / mean_rtn(i_mm)
			phi_n(i_mm) = (atan2(mean_rtn_sinn(i_mm), mean_rtn_cosn(i_mm))) / i_mm
		end do

	end if

!	i_mm = 2
!	write (*,*) N_ptc, mean_rtn_cosn(i_mm), mean_rtn_sinn(i_mm), mean_rtn(i_mm), eps_n(i_mm), phi_n(i_mm)

!	calculate momentum moments in center of mass frame with respect to tilted plane (phi_n)
	do i_mm = 1, 10
		mean_cosn_phip_phin(i_mm) = 0.d0
	end do

	do i_ptc = 1, N_ptc

		px = px_ptc(i_ptc)
		py = py_ptc(i_ptc)
		pz = pz_ptc(i_ptc)

		pt = sqrt(px*px + py*py)
		phi_pxpy = atan2(py, px)

		do i_mm = 1, 10
			phi_p = phi_pxpy - (phi_n(i_mm) + PI / i_mm)
			mean_cosn_phip_phin(i_mm) = mean_cosn_phip_phin(i_mm) + cos(i_mm*phi_p)
		end do

	end do

	do i_mm = 1, 10
		mean_cosn_phip_phin(i_mm) = mean_cosn_phip_phin(i_mm) / N_ptc
		v_n(i_mm) = mean_cosn_phip_phin(i_mm)
	end do

	return

	end subroutine
