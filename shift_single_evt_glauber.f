!	shift the density profile to center of mass frame
	subroutine shift_single_evt_glauber(impact_b, N_coll, N_part, N_charge, mean_x, mean_y, ecc_std, ecc_evp, eps_n, phi_n)

	implicit none

	include "parameters_common_blocks.inc"
			
	double precision impact_b
	
	double precision x, y, z
	
	integer i_atom, i_n1, i_n2 
	
	double precision N_coll, N_part, N_charge

	double precision part_x, part_xsq, part_y, part_ysq, part_xy
	double precision coll_x, coll_xsq, coll_y, coll_ysq, coll_xy

	double precision mean_x, mean_xsq, mean_y, mean_ysq, mean_xy
	double precision sigma_xsq, sigma_ysq, sigma_xy

	double precision ecc_std, ecc_evp

	double precision x1, y1, z1
	double precision x2, y2, z2
	double precision distance_sq, dpeevp_sq

	double precision frac, n_ch

	double precision dr

	integer i_mm
	double precision eps_n(10) 
	double precision phi_n(10) 

	integer i_rt

	double precision rt_sq, rt, phi, phi_xy
	double precision part_rtn(10), part_rtn_cosn(10), part_rtn_sinn(10)
	double precision coll_rtn(10), coll_rtn_cosn(10), coll_rtn_sinn(10)
	double precision mean_rtn(10), mean_rtn_cosn(10), mean_rtn_sinn(10)
	
	dr = 2.d0*nucleon_dr 	! nucleon size

	do i_mm = 1, 10
		phi_n(i_mm) = 0.d0
	end do

!	shift for nucleus 1
	do i_atom = 1, atom_A
		x_n1(i_atom) = x_n1(i_atom) - mean_x
		y_n1(i_atom) = y_n1(i_atom) - mean_y
		z_n1(i_atom) = z_n1(i_atom)
	end do
	
!	shift for nucleus 2 
	do i_atom = 1, atom_B
		x_n2(i_atom) = x_n2(i_atom) - mean_x
		y_n2(i_atom) = y_n2(i_atom) - mean_y 
		z_n2(i_atom) = z_n2(i_atom)
	end do

!	intialization before colliding		
	N_coll = 0.d0
	N_part = 0.d0

	part_x = 0.d0
	part_y = 0.d0
	part_xsq = 0.d0
	part_ysq = 0.d0
	part_xy = 0.d0
		
	coll_x = 0.d0
	coll_y = 0.d0
	coll_xsq = 0.d0
	coll_ysq = 0.d0
	coll_xy = 0.d0

	mean_x = 0.d0
	mean_y = 0.d0
	mean_xsq = 0.d0
	mean_ysq = 0.d0
	mean_xy = 0.d0
		
	sigma_xsq = 0.d0
	sigma_ysq = 0.d0
	sigma_xy = 0.d0

	ecc_std = 0.d0
	ecc_evp = 0.d0

	do i_mm = 1, 10
		part_rtn(i_mm) = 0.d0
		part_rtn_cosn(i_mm) = 0.d0
		part_rtn_sinn(i_mm) = 0.d0

		coll_rtn(i_mm) = 0.d0
		coll_rtn_cosn(i_mm) = 0.d0
		coll_rtn_sinn(i_mm) = 0.d0
	end do

!	now count participants and calculate moments of participant distribution
	do i_n1 = 1, atom_A
		x1 = x_n1(i_n1)
		y1 = y_n1(i_n1)
		z1 = z_n1(i_n1)

		if (part_n1(i_n1) .EQ. 1) then
			
			N_part = N_part + 1.d0
!			write (*,*) x1, y1, z1

			part_x = part_x + x1 
			part_y = part_y + y1
			part_xsq = part_xsq + x1**2 
			part_ysq = part_ysq + y1**2 
			part_xy = part_xy + x1*y1

			rt_sq = x1*x1 + y1*y1
			phi_xy = atan2(y1, x1)
			rt = sqrt(rt_sq)

			do i_mm = 1, 10

				i_rt = i_mm
				if (i_mm .eq. 1) i_rt = 3

				part_rtn(i_mm) = part_rtn(i_mm) + rt**i_rt
				phi = phi_xy - phi_n(i_mm)
				part_rtn_cosn(i_mm) = part_rtn_cosn(i_mm) + rt**i_rt * cos(i_mm*phi)
				part_rtn_sinn(i_mm) = part_rtn_sinn(i_mm) + rt**i_rt * sin(i_mm*phi)
			end do

		end if
	end do

	do i_n2 = 1, atom_B
		x2 = x_n2(i_n2)
		y2 = y_n2(i_n2)
		z2 = z_n2(i_n2)

		if (part_n2(i_n2) .EQ. 1) then

			N_part = N_part + 1.d0
!			write (*,*) x2, y2, z2
			
			part_x = part_x + x2 
			part_y = part_y + y2
			part_xsq = part_xsq + x2**2 
			part_ysq = part_ysq + y2**2 
			part_xy = part_xy + x2*y2

			rt_sq = x2*x2 + y2*y2
			phi_xy = atan2(y2, x2)
			rt = sqrt(rt_sq)

			do i_mm = 1, 10

				i_rt = i_mm
				if (i_mm .eq. 1) i_rt = 3

				part_rtn(i_mm) = part_rtn(i_mm) + rt**i_rt
				phi = phi_xy - phi_n(i_mm)
				part_rtn_cosn(i_mm) = part_rtn_cosn(i_mm) + rt**i_rt * cos(i_mm*phi)
				part_rtn_sinn(i_mm) = part_rtn_sinn(i_mm) + rt**i_rt * sin(i_mm*phi)
			end do

		end if
	end do

!	count binary collisions and calculate moments of binary collision distribution
	do i_n1 = 1, atom_A
		x1 = x_n1(i_n1)
		y1 = y_n1(i_n1)
		z1 = z_n1(i_n1)

		do i_n2 = 1, atom_B
			x2 = x_n2(i_n2)
			y2 = y_n2(i_n2)
			z2 = z_n2(i_n2)

			x = (x1+x2)/2.d0
			y = (y1+y2)/2.d0
			z = (z1+z2)/2.d0

			x_coll(i_n1, i_n2) = x
			y_coll(i_n1, i_n2) = y
			z_coll(i_n1, i_n2) = z

			if (coll_n1_n2(i_n1, i_n2) .EQ. 1) then

				N_coll = N_coll + 1.d0

				coll_x = coll_x + x
				coll_y = coll_y + y
				coll_xsq = coll_xsq + x**2
				coll_ysq = coll_ysq + y**2
				coll_xy = coll_xy + x*y

				rt_sq = x*x + y*y
				phi_xy = atan2(y, x)
				rt = sqrt(rt_sq)

				do i_mm = 1, 10

					i_rt = i_mm
					if (i_mm .eq. 1) i_rt = 3

					coll_rtn(i_mm) = coll_rtn(i_mm) + rt**i_rt
					phi = phi_xy - phi_n(i_mm)
					coll_rtn_cosn(i_mm) = coll_rtn_cosn(i_mm) + rt**i_rt * cos(i_mm*phi)
					coll_rtn_sinn(i_mm) = coll_rtn_sinn(i_mm) + rt**i_rt * sin(i_mm*phi)
				end do

			end if
		end do
	end do

!	now come to different models

!	N_part scaling with centrality dependence
	if (flag_glauber .eq. 11) then

		frac = frac_coll 
		n_ch = mean_npp

		if (N_part .gt. 0.d0) then

			N_charge = (frac * N_coll/N_part + (1.d0-frac)/2.d0) * N_part * n_ch

			mean_x = (frac * N_coll/N_part + (1.d0-frac)/2.d0) * part_x * n_ch
			mean_y = (frac * N_coll/N_part + (1.d0-frac)/2.d0) * part_y * n_ch
			mean_xsq = (frac * N_coll/N_part + (1.d0-frac)/2.d0) * part_xsq * n_ch
			mean_ysq = (frac * N_coll/N_part + (1.d0-frac)/2.d0) * part_ysq * n_ch
			mean_xy = (frac * N_coll/N_part + (1.d0-frac)/2.d0) * part_xy * n_ch 

			do i_mm = 1, 10
				mean_rtn(i_mm) = (frac * N_coll/N_part + (1.d0-frac)/2.d0) * part_rtn(i_mm) * n_ch
				mean_rtn_cosn(i_mm) = (frac * N_coll/N_part + (1.d0-frac)/2.d0) * part_rtn_cosn(i_mm) * n_ch
				mean_rtn_sinn(i_mm) = (frac * N_coll/N_part + (1.d0-frac)/2.d0) * part_rtn_sinn(i_mm) * n_ch
			end do

		end if

	end if

!	2-component scaling
	if (flag_glauber .eq. 21) then

		frac = frac_coll
		n_ch = mean_npp

		N_charge = (frac * N_coll + (1.d0-frac)/2.d0 * N_part) * n_ch

		mean_x = (frac * coll_x + (1.d0-frac)/2.d0 * part_x) * n_ch
		mean_y = (frac * coll_y + (1.d0-frac)/2.d0 * part_y) * n_ch
		mean_xsq = (frac * coll_xsq + (1.d0-frac)/2.d0 * part_xsq) * n_ch
		mean_ysq = (frac * coll_ysq + (1.d0-frac)/2.d0 * part_ysq) * n_ch
		mean_xy = (frac * coll_xy + (1.d0-frac)/2.d0 * part_xy) * n_ch 

		do i_mm = 1, 10
			mean_rtn(i_mm) = (frac * coll_rtn(i_mm)  + (1.d0-frac)/2.d0 * part_rtn(i_mm)) * n_ch
			mean_rtn_cosn(i_mm) = (frac * coll_rtn_cosn(i_mm) + (1.d0-frac)/2.d0 * part_rtn_cosn(i_mm)) * n_ch
			mean_rtn_sinn(i_mm) = (frac * coll_rtn_sinn(i_mm) + (1.d0-frac)/2.d0 * part_rtn_sinn(i_mm)) * n_ch
		end do

	end if

!	N_coll scaling with centrality dependence 
	if (flag_glauber .eq. 31) then

		frac = frac_coll
		n_ch = mean_npp

		if (N_coll .gt. 0) then

			N_charge = (frac + (1.d0-frac)/2.d0 * N_part/N_coll) * N_coll * n_ch

			mean_x = (frac + (1.d0-frac)/2.d0 * N_part/N_coll) * coll_x * n_ch
			mean_y = (frac + (1.d0-frac)/2.d0 * N_part/N_coll) * coll_y * n_ch
			mean_xsq = (frac + (1.d0-frac)/2.d0 * N_part/N_coll) * coll_xsq * n_ch
			mean_ysq = (frac + (1.d0-frac)/2.d0 * N_part/N_coll) * coll_ysq * n_ch
			mean_xy = (frac + (1.d0-frac)/2.d0 * N_part/N_coll) * coll_xy * n_ch 

			do i_mm = 1, 10
				mean_rtn(i_mm) = (frac * N_coll/N_part + (1.d0-frac)/2.d0) * coll_rtn(i_mm) * n_ch
				mean_rtn_cosn(i_mm) = (frac * N_coll/N_part + (1.d0-frac)/2.d0) * coll_rtn_cosn(i_mm) * n_ch
				mean_rtn_sinn(i_mm) = (frac * N_coll/N_part + (1.d0-frac)/2.d0) * coll_rtn_sinn(i_mm) * n_ch
			end do

		end if
	end if

!	calculate the mean, variace, and eccentricities
	if (N_part .GT. 0.d0 .AND. N_charge .GT. 0.d0) then

		mean_x = mean_x / N_charge
		mean_y = mean_y / N_charge
		mean_xsq = mean_xsq / N_charge
		mean_ysq = mean_ysq / N_charge
		mean_xy = mean_xy / N_charge
		
		sigma_xsq = mean_xsq - mean_x**2
		sigma_ysq = mean_ysq - mean_y**2
		sigma_xy = mean_xy - mean_x * mean_y

		ecc_std = (mean_ysq - mean_xsq) / (mean_ysq + mean_xsq)
		ecc_evp = sqrt((sigma_xsq - sigma_ysq)**2 + 4.d0*sigma_xy**2) / (sigma_xsq + sigma_ysq)

!		if there is only a single collision point or only one of collisions produces charged particles
		if (mean_xsq .eq. 0.d0 .AND. mean_ysq .eq. 0.d0) ecc_std = 0.d0
		if (mean_xsq + mean_ysq .eq. 0.d0) ecc_std = 0.d0
		if (sigma_xsq .eq. 0.d0 .AND. sigma_ysq .eq. 0.d0) ecc_evp = 0.d0
		if (sigma_xsq + sigma_ysq .eq. 0.d0) ecc_evp = 0.d0
 
!		write (*,*) mean_x, mean_y, mean_xsq, mean_ysq, mean_xy 
!		write (*,*) sigma_xsq, sigma_ysq, sigma_xy

!		write (*,*) ecc_std, ecc_evp

		do i_mm = 1, 10
			mean_rtn(i_mm) = mean_rtn(i_mm) / N_charge
			mean_rtn_cosn(i_mm) = mean_rtn_cosn(i_mm) / N_charge
			mean_rtn_sinn(i_mm) = mean_rtn_sinn(i_mm) / N_charge

!			write (*,*) i_mm, mean_rtn_cosn(i_mm), mean_rtn_sinn(i_mm), mean_rtn(i_mm)
		end do

		do i_mm = 1, 10
!			this gives the standard eps_n
!			eps_n(i_mm) = mean_rtn_cosn(i_mm) / mean_rtn(i_mm)

!			this gives eps_n after rotation
			eps_n(i_mm) = sqrt(mean_rtn_cosn(i_mm)**2 + mean_rtn_sinn(i_mm)**2) / mean_rtn(i_mm)

			phi_n(i_mm) = (atan2(mean_rtn_sinn(i_mm), mean_rtn_cosn(i_mm))) / i_mm
		end do

	end if

	end subroutine
