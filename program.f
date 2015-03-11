	
	program main

	implicit none

	include "parameters_common_blocks.inc"
		
	double precision tmp 

	double precision impact_b

	integer i
	integer ix, iy, ieta

	double precision x, y, z, eta
	double precision x2, y2, z2, eta2

	double precision r
	double precision rho_A
	integer atom

	double precision dr
	
	integer i_atom, i_n1, i_n2 

	double precision N_coll, N_part, N_charge, N_ptc
	double precision Nsq_coll, Nsq_part, Nsq_charge, Nsq_ptc
	double precision sigma_Ncoll, sigma_Npart, sigma_Ncharge, sigma_Nptc

	double precision Ncoll_evt, Npart_evt, Ncharge_evt
	integer intNptc_evt, intNptc_rpdcut_evt
	double precision Edptc_evt, Edptc_rpdcut_evt

	double precision Nptc_tot, Nptc_rpdcut_tot

	integer i_evt, ntot_evt, ntrptc_evt
	
	double precision sigma_xsq, sigma_ysq, sigma_xy
	
	double precision ecc_std, ecc_evp
	double precision eccsq_std, eccsq_evp
	double precision sigma_eccstd, sigma_eccevp
	double precision eccstd_evt, eccevp_evt

	double precision meanx_evt, meany_evt

	integer i_mm
	double precision eps_n(10)
	double precision epsn_evt(10)
	double precision phin_evt(10)
	double precision v_n(10)
	double precision vn_evt(10)

	double precision pp_multiplicity_dist, n_ch
	double precision part_pp_multiplicity_dist, part_nch
	double precision coll_pp_multiplicity_dist, coll_nch
	integer n

	double precision cdf(0:1000), part_cdf(0:1000), coll_cdf(0:1000)

	double precision eps_mom_lab(0:3)
	double precision mtg(0:3, 0:3)
	double precision emt(0:3, 0:3) 
 
        integer i_rpd, n_rpd
        double precision rpd_max, d_rpd

	double precision b_min, b_max

	integer iosend

	character (len = 100) filename

!       please the array sizes are big enough
	double precision CDFntrevt_NEptc(0:10000)
	double precision ntrevt_NEptc(0:10000), ecc_NEptc(0:10000,10)
	double precision NEptc, NEptc_min, NEptc_max, d_NEptc
	integer n_NEptc, i_NEptc, iii
	
        Nptc_tot = 0.d0
        Nptc_rpdcut_tot = 0.d0

	NEptc_min = 0.d0
	NEptc_max = 10000.d0
	d_NEptc = 1.d0
	n_NEptc = nint((NEptc_max-NEptc_min)/d_NEptc)

!	qin
9999	open (99999,  file="/dev/urandom", 
     $		access="stream", form="unformatted", iostat = iosend)

	if (iosend .eq. 0) then
		read (99999) idum
		write (*,*) "random seed from noise!"
	else
		idum = time()
		write (*,*) "random seed from time!"
	end if
	if (idum .eq. 0) go to 9999
	if (idum .gt. 0) idum = -idum
!        idum = -1111
	write (*,*) idum

!	beginning of code
!	beginning of code
!	beginning of code

!	the _beam rapdity
!	y_beam = 8.58d0
!	y_beam = 5.36d0
	E_beam = 0.5*sqrt_sNN
	M_beam = 0.938d0
	P_beam = sqrt(E_beam*E_beam - M_beam*M_beam)
	y_beam = 0.5d0 * log((E_beam + P_beam)/(E_beam - P_beam))
	write (*,*) "rapidity:", E_beam, P_beam, y_beam, eta_0, sigma_eta

	dr = 2.d0*nucleon_dr 	! nucleon size (hard sphere for nucleons)

!	initialize the nuclear density distribution
	atom = max(atom_A, atom_B)
	r = 0.d0	
	tmp = rho_A(r, atom)
!	write (*,*) "nuclear distribution function:", atom, r, tmp

!!	test
!	atom = atom_A
!	do n = 1, 100000
!		call mc_rho_A(x, y, z, atom_A)
!		x = 2.d0*x
!		y = 2.d0*y
!		z = 2.d0*z
!		r = sqrt(x*x + y*y + z*z)
!		write (*,*) x, y, z
!	end do
!
!	return

!	initialize the cdf of the multiplicity distribution
	call cdf_pp_multiplicity_dist(cdf_nb, part_cdf_nb, coll_cdf_nb)

!!	check pp multiplicity distribution
!	do n = 0, 100
!		cdf(n) = 0.d0
!	end do
!
!	open (9, file="dt_ebe_nch_pp.dat")
!	do i_evt = 1, n_evt
!		call mc_pp_multiplicity_dist(n_ch, part_nch)
!		write (9,*) n_ch
!		n = anint(n_ch)
!		cdf(n) = cdf(n) + 1
!	end do
!	close (9)
	
	open (9, file="dt_pp_multiplicity.dat")
	do n = 0, 100
		n_ch = n
!		write (*,*) n, cdf(n)/n_evt, pp_multiplicity_dist(n_ch)
		write (9,*) n, pp_multiplicity_dist(n_ch)
	end do
	close (9)

!	get the maximum value for fun_eta
	call get_fun_eta_max()
        write (*,*) fun_eta_max_right, fun_eta_max_left, fun_eta_max_lrboth

!	initialize initial conditions
	do ix = -nx, nx
	do iy = -ny, ny
	do ieta = -neta, neta
	
		do i = 0, 3
			ep_density(i, ix, iy, ieta) = 0.d0
			ep_density_tot(i, ix, iy, ieta) = 0.d0

			v_flow(i, ix, iy, ieta) = 0.d0
			v_flow_tot(i, ix, iy, ieta) = 0.d0
		end do
	
	end do
	end do
	end do

        do ieta = -neta, neta
        
        	nevt_eta(ieta) = 0.d0
                nevt_eta_tot = 0.d0
                nevt_eta_unit = 0.d0
        
        end do

!	beginning nucleus-nucleus collisions
	b_min = 0.d0
	b_max = 2.5d0*radius
!	b_max = 2.d0*radius

	open (10, file="dt_ebe_b_npart_ncoll_ncharge_nptc.dat")
        open (15, file="dt_all_binary_collisions.dat")

!	do loop in impact parameters
!	impact_b = 0.d0
!	do while (impact_b .lt. 2.1d0*radius)

!		initialization before calling glauber
	
		ntot_evt = 0
		ntrptc_evt = 0

		N_coll = 0.d0
		N_part = 0.d0
		N_charge = 0.d0
		N_ptc = 0.d0

		Nsq_coll = 0.d0
		Nsq_part = 0.d0
		Nsq_charge = 0.d0
		Nsq_ptc = 0.d0
		
		ecc_std = 0.d0
		ecc_evp = 0.d0

		eccsq_std = 0.d0
		eccsq_evp = 0.d0
	
		do i_mm = 1, 10

			eps_n(i_mm) = 0.d0
			v_n(i_mm) = 0.d0

		end do

		do i_NEptc = 0, n_NEptc
		
			ntrevt_NEptc(i_NEptc) = 0.d0
			CDFntrevt_NEptc(i_NEptc) = 0.d0
			do i_mm = 1, 10
				ecc_NEptc(i_NEptc,i_mm) = 0.d0
			end do
		
		end do

!		do loop in events
100		do i_evt = 1, n_evt

			if (ntrptc_evt  .ge. ntrptcevt_max) exit

			ntot_evt = ntot_evt + 1

			if (b_min .eq. b_max) then
				impact_b = b_min
			else
				tmp = ran2(idum)
				impact_b = sqrt(b_min**2 + (b_max**2 - b_min**2) * tmp)
			end if

!			write (*,*) "one glauber event:", i_evt, impact_b
			call single_evt_glauber(impact_b, Ncoll_evt, Npart_evt, Ncharge_evt, meanx_evt, meany_evt, eccstd_evt, eccevp_evt)
!			write (*,*) Ncoll_evt, Npart_evt, Ncharge_evt, meanx_evt, meany_evt, eccstd_evt, eccevp_evt

			if (Npart_evt .lt. 1.d0) cycle

!!			output nucleon positions and if they participate the collision
!			open (11, file="dt_nucleus1.dat")
!			open (12, file="dt_part_nucleus1.dat")
!			do i_atom = 1, atom_A
!				if (i_evt .eq. 1) then
!
!					write (11, *) x_n1(i_atom), y_n1(i_atom)
!					if (part_n1(i_atom) .eq. 1) write (12, *) x_n1(i_atom), y_n1(i_atom)
!
!				end if
!			end do
!			close (11)
!			close (12)

!			open (13, file="dt_nucleus2.dat")
!			open (14, file="dt_part_nucleus2.dat")
!			do i_atom = 1, atom_B
!				if (i_evt .eq. 1) then
!
!					write (13, *) x_n2(i_atom), y_n2(i_atom)
!					if (part_n2(i_atom) .eq. 1) write (14, *) x_n2(i_atom), y_n2(i_atom)
!
!				end if
!			end do
!			close (13)
!			close (14)

!			record every event for later use
!			do i_atom = 1, atom_A
!
!				xxx_n1(i_evt, i_atom) = x_n1(i_atom)
!				yyy_n1(i_evt, i_atom) = y_n1(i_atom)
!				zzz_n1(i_evt, i_atom) = z_n1(i_atom)
!
!				participant_n1(i_evt, i_atom) = part_n1(i_atom)
!
!			end do
!
!			do i_atom = 1, atom_B
!
!				xxx_n2(i_evt, i_atom) = x_n2(i_atom)
!				yyy_n2(i_evt, i_atom) = y_n2(i_atom)
!				zzz_n2(i_evt, i_atom) = z_n2(i_atom)
!
!				participant_n2(i_evt, i_atom) = part_n2(i_atom)
!
!			end do
!
!			do i_n1 = 1, atom_A
!				do i_n2 = 1, atom_B
!
!				xxx_coll(i_evt, i_n1, i_n2) = x_coll(i_n1, i_n2)
!				yyy_coll(i_evt, i_n1, i_n2) = y_coll(i_n1, i_n2)
!				zzz_coll(i_evt, i_n1, i_n2) = z_coll(i_n1, i_n2)
!
!				collision_n1_n2(i_evt, i_n1, i_n2) = coll_n1_n2(i_n1, i_n2)
!
!				end do
!			end do


!			write (*,*) "shift glauber event to center of mass frame:"
			call shift_single_evt_glauber(impact_b, Ncoll_evt, Npart_evt, Ncharge_evt, meanx_evt, meany_evt, eccstd_evt, eccevp_evt, epsn_evt, phin_evt)
!			write (*,*) Ncoll_evt, Npart_evt, Ncharge_evt, meanx_evt, meany_evt, eccstd_evt, eccevp_evt
!			write (*,*) "eps_n:", epsn_evt(1), epsn_evt(2), epsn_evt(3), epsn_evt(4)

!			write (*,*) "get particles and geometry:"
			call get_particle_single_evt_glauber(intNptc_evt, Edptc_evt, intNptc_rpdcut_evt, Edptc_rpdcut_evt)
!			write (*,*) intNptc_evt, Edptc_evt, meanx_evt, meany_evt

!			select ultra-central collisions (a_coll = 0., nrpdy = 1)
			if (abs(n_rpd_unit - 2.d0*rpd_cut) .lt. 0.01d0) then

!			if (intNptc_evt. lt. 42) cycle		! 0-40%
!			if (intNptc_evt. lt. 65) cycle		! 0-20%
!			if (intNptc_evt. lt. 83) cycle		! 0-10%
!			if (intNptc_evt. lt. 100) cycle		! 0-5%
!			if (intNptc_evt. lt. 119) cycle		! 0-2%
!			if (intNptc_evt. lt. 132) cycle		! 0-1%
!			if (intNptc_evt. lt. 172) cycle		! 0-0.1%
!			if (intNptc_evt. lt. 208) cycle		! 0-0.01%
!			if (intNptc_evt. lt. 242) cycle		! 0-0.001%


!			if (intNptc_evt.ge.25.or.intNptc_evt.lt.6) cycle		! 60-90%
!			if (intNptc_evt.ge.42.or.intNptc_evt.lt.25) cycle		! 40-60%
!			if (intNptc_evt.ge.52.or.intNptc_evt.lt.42) cycle		! 30-40%
!			if (intNptc_evt.ge.65.or.intNptc_evt.lt.52) cycle		! 20-30%
!			if (intNptc_evt.ge.83.or.intNptc_evt.lt.65) cycle		! 10-20%
!			if (intNptc_evt.ge.100.or.intNptc_evt.lt.83) cycle		! 5-10%
!			if (intNptc_evt.ge.132.or.intNptc_evt.lt.100) cycle		! 1-5%
			if (intNptc_evt.lt.132) cycle		! 0-1%
 
                        else

        			stop

			end if

                        Nptc_tot = Nptc_tot + intNptc_evt
                        Nptc_rpdcut_tot = Nptc_rpdcut_tot + intNptc_rpdcut_evt

!			output all participants positions (only for accepted events)
!				
!			open (11, file="dt_nucleus1.dat")
!			open (12, file="dt_part_nucleus1.dat")
!			do i_atom = 1, atom_A
!				write (11, *) x_n1(i_atom), y_n1(i_atom)
!				if (part_n1(i_atom) .eq. 1) write (12, *) x_n1(i_atom), y_n1(i_atom)
!			end do
!			close (11)
!			close (12)
!
!			open (13, file="dt_nucleus2.dat")
!			open (14, file="dt_part_nucleus2.dat")
!			do i_atom = 1, atom_B
!				write (13, *) x_n2(i_atom), y_n2(i_atom)
!				if (part_n2(i_atom) .eq. 1) write (14, *) x_n2(i_atom), y_n2(i_atom)
!			end do
!			close (13)
!			close (14)
!

!			write (*,*) "shift to center of mass frame:"
!			call shift_particle_single_evt(intNptc_evt, meanx_evt, meany_evt)
!			write (*,*) intNptc_evt, Edptc_evt, meanx_evt, meany_evt

!			write (*,*) "get the geometry:"
!			call get_geometry_particle(intNptc_evt, meanx_evt, meany_evt, epsn_evt, phin_evt, vn_evt)
!			write (*,*) "eps_n:", epsn_evt(1), epsn_evt(2), epsn_evt(3), epsn_evt(4)
!			write (*,*) "phi_n:", phin_evt(1), phin_evt(2), phin_evt(3), phin_evt(4)

!			write (*,*) "evolve particles and geometry:"
			call evolve_particle_single_evt_glauber_taueta(intNptc_evt)

!			write (*,*) "shift to center of mass frame:"
!			call shift_particle_single_evt(intNptc_evt, meanx_evt, meany_evt)
!			write (*,*) intNptc_evt, Edptc_evt, meanx_evt, meany_evt

!			write (*,*) "get the geometry after free-streaming:"
			call get_geometry_particle(intNptc_evt, meanx_evt, meany_evt, epsn_evt, phin_evt, vn_evt)
!			write (*,*) "eps_n:", epsn_evt(1), epsn_evt(2), epsn_evt(3), epsn_evt(4)
!			write (*,*) "phi_n:", phin_evt(1), phin_evt(2), phin_evt(3), phin_evt(4)

!			write (*,*) "get the geometry after rotation:"
			call rotate_particle_single_evt(intNptc_evt, meanx_evt, meany_evt, epsn_evt, phin_evt, vn_evt)
!			write (*,*) "eps_n:", epsn_evt(1), epsn_evt(2), epsn_evt(3), epsn_evt(4)
!			write (*,*) "phi_n:", phin_evt(1), phin_evt(2), phin_evt(3), phin_evt(4)


			N_coll = N_coll + Ncoll_evt
			N_part = N_part + Npart_evt
			N_charge = N_charge + Ncharge_evt
			N_ptc = N_ptc + intNptc_evt

			Nsq_coll = Nsq_coll + Ncoll_evt**2
			Nsq_part = Nsq_part + Npart_evt**2
			Nsq_charge = Nsq_charge + Ncharge_evt**2
			Nsq_ptc = Nsq_ptc + intNptc_evt**2
	
			if (intNptc_evt .GT. 0.d0) then

				if (intNptc_evt .eq. 2) then
					write (*,*) "Hey:", intNptc_evt, meanx_evt, meany_evt, epsn_evt(1), epsn_evt(2), epsn_evt(3)
					write (*,*) x_ptc(1), x_ptc(2), y_ptc(1), y_ptc(2)

				end if

				ntrptc_evt = ntrptc_evt + 1
                                write (*,*) ntrptc_evt

				ecc_std = ecc_std + eccstd_evt
				ecc_evp = ecc_evp + eccevp_evt

				eccsq_std = eccsq_std + eccstd_evt**2
				eccsq_evp = eccsq_evp + eccevp_evt**2

!                               we first check Nptc distribution without rpd_cut (i.e., total number of ptc), to determine mean value and slope parameter
!                               we then apply the rpd_cut and increases the mean number by a factor, and check Nptc_rpdcut distribution (to determine this factor)
				i_NEptc = int((dble(intNptc_evt) - NEptc_min) / d_NEptc)
!				i_NEptc = int((dble(Edptc_evt) - NEptc_min) / d_NEptc)
!				i_NEptc = int((dble(intNptc_rpdcut_evt) - NEptc_min) / d_NEptc)
!				i_NEptc = int((dble(Edptc_rpdcut_evt) - NEptc_min) / d_NEptc)
                                if (i_NEptc .gt. n_NEptc) then
                                write (*,*) "NEptc array size 2 small!"
                                stop
                                end if

				ntrevt_NEptc(i_NEptc) = ntrevt_NEptc(i_NEptc) + 1.d0

				do iii = 0, i_NEptc
					CDFntrevt_NEptc(iii) = CDFntrevt_NEptc(iii) + 1.d0
				end do

				do i_mm = 1, 10

					eps_n(i_mm) = eps_n(i_mm) + epsn_evt(i_mm)
					v_n(i_mm) = v_n(i_mm) + vn_evt(i_mm)

					ecc_NEptc(i_NEptc, i_mm) = ecc_NEptc(i_NEptc, i_mm) +  epsn_evt(i_mm)

				end do

			!	write (10,*) ntrptc_evt, impact_b, Ncoll_evt, Npart_evt, Ncharge_evt, intNptc_evt, Edptc_evt

                                !skip the following when doing centrality determination
                                !cycle
                                
				do ix = -nx, nx
				do iy = -ny, ny
				do ieta = -neta, neta
	
					x = ix * dx
					y = iy * dy
					eta = ieta * deta

					call get_lab_dist_parameters_taueta(intNptc_evt, x, y, eta, eps_mom_lab, emt)

					do i = 0, 3
						ep_density(i, ix, iy, ieta) = eps_mom_lab(i)
						ep_density_tot(i, ix, iy, ieta) = ep_density_tot(i, ix, iy, ieta) + ep_density(i, ix, iy, ieta)
					end do

!                                       simplest way to calculate initial flow (for testing)
!                                       the sophisticated way is to solve the eigen-value equations
                                        if (eps_mom_lab(0) .gt. 0.d0) then
                                               do i = 1, 3
                                                        v_flow(i, ix, iy, ieta) = ep_density(i, ix, iy, ieta) / ep_density(0, ix, iy, ieta)
                                               end do
                                        end if

			        end do
			        end do
!			        write (*,*) 100*(ix+nx)/(2*nx), "% finished!", ntrptc_evt
			        end do

!	    		        output binary collisions positions (only for accepted events)
                                do i_n1 = 1, atom_A
                                do i_n2 = 1, atom_B
                                   if (coll_n1_n2(i_n1, i_n2) .EQ. 1) then
                                      x = x_coll(i_n1, i_n2)
                                      y = y_coll(i_n1, i_n2)
                                      write (15,*) x, y
                                   end if
                                end do
                                end do

			end if

		end do

!		average over events
		
		N_coll = N_coll / ntrptc_evt
		N_part = N_part / ntrptc_evt
		N_charge = N_charge / ntrptc_evt
		N_ptc = N_ptc / ntrptc_evt

		Nsq_coll = Nsq_coll / ntrptc_evt
		Nsq_part = Nsq_part / ntrptc_evt
		Nsq_charge = Nsq_charge / ntrptc_evt
		Nsq_ptc = Nsq_ptc / ntrptc_evt

		sigma_Ncoll = sqrt(Nsq_coll - N_coll**2)
		sigma_Npart = sqrt(Nsq_part - N_part**2)
		sigma_Ncharge = sqrt(Nsq_charge - N_charge**2)
		sigma_Nptc = sqrt(Nsq_ptc - N_ptc**2)

!		for eccenetricity we need to average over true events (if no collision happens, no medium, eccentricity is undetermined)
		ecc_std = ecc_std / ntrptc_evt
		ecc_evp = ecc_evp / ntrptc_evt

		eccsq_std = eccsq_std / ntrptc_evt
		eccsq_evp = eccsq_evp / ntrptc_evt

		sigma_eccstd = sqrt(eccsq_std - ecc_std**2)
		sigma_eccevp = sqrt(eccsq_evp - ecc_evp**2)

		do i_mm = 1, 10

			eps_n(i_mm) = eps_n(i_mm) / ntrptc_evt
			v_n(i_mm) = v_n(i_mm) / ntrptc_evt

		end do
		
	        open (41, file="dt_eccentricity_vs_NEptc.dat")

		do i_NEptc = 0, n_NEptc
		
			do i_mm = 1, 10
				ecc_NEptc(i_NEptc,i_mm) = ecc_NEptc(i_NEptc,i_mm) / ntrevt_NEptc(i_NEptc) 
			end do

			NEptc = (i_NEptc+0.5d0)*d_NEptc
			NEptc = i_NEptc*d_NEptc

			if (ntrevt_NEptc(i_NEptc) .gt. 0) then
				write (41,*) NEptc, dble(ntrevt_NEptc(i_NEptc))/ntrptc_evt/d_NEptc, dble(CDFntrevt_NEptc(i_NEptc))/ntrptc_evt, ecc_NEptc(i_NEptc,2), ecc_NEptc(i_NEptc,3)
			end if

		end do
                close (41)

	        open (42, file="dt_avr_Ntot_rpdcut_to_Ntot.dat")
                write (*,*) Nptc_tot/ntrptc_evt, Nptc_rpdcut_tot/ntrptc_evt, Nptc_rpdcut_tot/Nptc_tot, Nptc_tot/Nptc_rpdcut_tot
                write (42,*) Nptc_tot/ntrptc_evt, Nptc_rpdcut_tot/ntrptc_evt, Nptc_rpdcut_tot/Nptc_tot, Nptc_tot/Nptc_rpdcut_tot
                close (42)

!                return

!		finalize/average initial conditions

		open (50, file="dt_hydro_input.dat")

		open (55, file="dt_lab_2dhydro_ic_ed.dat")
		open (56, file="dt_lab_2dhydro_ic_vflow.dat")
		open (57, file="dt_lab_1dhydro_ic_x.dat")
		open (58, file="dt_lab_1dhydro_ic_y.dat")
		open (59, file="dt_lab_1dhydro_ic_eta.dat")

		open (65, file="dt_lab_2dhydro_xeta_ic_ed.dat")
		open (66, file="dt_lab_2dhydro_xeta_ic_vflow.dat")

		open (75, file="dt_lab_2dhydro_yeta_ic_ed.dat")
		open (76, file="dt_lab_2dhydro_yeta_ic_vflow.dat")

		do ix = -nx, nx
		do iy = -ny, ny
		do ieta = -neta, neta
		
			x = ix * dx
			y = iy * dy
			eta = ieta * deta
	
			do i = 0, 3
				ep_density(i, ix, iy, ieta) = ep_density_tot(i, ix, iy, ieta) / ntrptc_evt
			end do
	
!                       simplest way to calculate initial flow (for testing)
!                       the sophisticated way is to solve the eigen-value equations
                        if (ep_density(i, ix, iy, ieta) .gt. 0.d0) then
                                do i = 1, 3
                                        v_flow(i, ix, iy, ieta) = ep_density(i, ix, iy, ieta) / ep_density(0, ix, iy, ieta)
                                end do
                        end if

!	                write (50,"(3i5,5e14.6)") ix, iy, ieta, ep_density(0, ix, iy, ieta), v_flow(1, ix, iy, ieta), v_flow(2, ix, iy, ieta), v_flow(3, ix, iy, ieta), 0.d0
!	                write (50,*) ix, iy, ieta, ep_density(0, ix, iy, ieta), v_flow(1, ix, iy, ieta), v_flow(2, ix, iy, ieta), v_flow(3, ix, iy, ieta), 0.d0
!                       what are needed: T^{\tau \tau}, T^{\tau x}, T^{\tau y}, T^{\tau \eta}
	                write (50,*) ix, iy, ieta, ep_density(0, ix, iy, ieta), ep_density(1, ix, iy, ieta), ep_density(2, ix, iy, ieta), ep_density(3, ix, iy, ieta), 0.d0

!	                if (ieta .eq. 0) then
	                if (ieta .eq. Neta/2) then
	                	write (55,*) x, y, ep_density(0, ix, iy, ieta)
	                	write (56,*) x, y, v_flow(1, ix, iy, ieta), v_flow(2, ix, iy, ieta), v_flow(3, ix, iy, ieta)
	                end if

	                if (ieta .eq. 0 .AND. iy .eq. 0) then
	                	write (57,*) x, ep_density(0, ix, iy, ieta), v_flow(1, ix, iy, ieta), v_flow(2, ix, iy, ieta), v_flow(3, ix, iy, ieta)
	                end if

	                if (ieta .eq. 0 .AND. ix .eq. 0) then
	                	write (58,*) y, ep_density(0, ix, iy, ieta), v_flow(1, ix, iy, ieta), v_flow(2, ix, iy, ieta), v_flow(3, ix, iy, ieta)
	                end if

	                if (ix .eq. 0 .AND. iy .eq. 0) then
	                	write (59,*) eta, ep_density(0, ix, iy, ieta), v_flow(1, ix, iy, ieta), v_flow(2, ix, iy, ieta), v_flow(3, ix, iy, ieta)
	                end if

!	                if (iy .eq. 0) then
	                if (iy .eq. Ny/2) then
	                	write (65,*) x, eta, ep_density(0, ix, iy, ieta)
	                	write (66,*) x, eta, v_flow(1, ix, iy, ieta), v_flow(2, ix, iy, ieta), v_flow(3, ix, iy, ieta)
	                end if

!	                if (ix .eq. 0) then
	                if (ix .eq. Nx/2) then
	                	write (75,*) y, eta, ep_density(0, ix, iy, ieta)
	                	write (76,*) y, eta, v_flow(1, ix, iy, ieta), v_flow(2, ix, iy, ieta), v_flow(3, ix, iy, ieta)
	                end if

		end do
		end do
		end do

		close (50)

		close (55)
		close (56)
		close (57)
		close (58)
		close (59)

                close (15)

	        open (16, file="dt_avr_initial_rapidity.dat")
                do ieta = -neta, neta
                
                        eta = ieta * deta
                        write (16, *) eta, nevt_eta(ieta)/nevt_eta_tot/deta, nevt_eta_unit/nevt_eta_tot, nevt_eta_tot/nevt_eta_unit
                
                end do
                close (16)

		write (*,*) "code works"
!		compare to average N_charge and N_ptc, we can check that
!		if we increase mean_npp by a factor of 10
!		then we need to increase k_pp by a factor of 10, to get N_pt = N_charge
		write (*,*) impact_b, N_part, N_coll, N_charge, N_ptc, N_charge/(N_part/2.d0), N_ptc/(N_part/2.d0)

!	        output npart, ncoll, ncharge, eccentricity
	        open (20, file="dt_avr_ncoll_npart_ncharge_nptc.dat")
	        open (21, file="dt_avr_eccentricity.dat")
	        open (22, file="dt_avr_ncharge_nptc_per_pair.dat")

		write (20,*) impact_b, N_part, sigma_Npart, N_coll, sigma_Ncoll, N_charge, sigma_Ncharge
		write (21,*) impact_b, N_part, N_coll, ecc_std, sigma_eccstd, ecc_evp, sigma_eccevp
		write (22,*) impact_b, N_part, N_coll, N_charge/(N_part/2.d0), N_ptc/(N_part/2.d0)

                close (20)
                close (21)
                close (22)

!	        output for geometry moments
	        open (30, file="dt_avr_ptc_spatial_moments.dat")
	        open (31, file="dt_avr_ptc_sep_momentum_moments.dat")

		write (30,*) impact_b, eps_n(1), eps_n(2), eps_n(3), eps_n(4), eps_n(5)
		write (31,*) impact_b, v_n(1), v_n(2), v_n(3), v_n(4), v_n(5)

                close (30)
                close (31)

!		impact_b = impact_b + 5.d0*dr
!
!	end do

        close (10)

	end program
