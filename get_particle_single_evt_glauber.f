!	create all particles in the center of mass frame
	subroutine get_particle_single_evt_glauber(N_ptc, ET_ptc, Nptc_rpdcut, ETptc_rpdcut)

	implicit none

	include "parameters_common_blocks.inc"
			
	double precision x, y, z
	
	integer i_atom, i_n1, i_n2 

	double precision x1, y1, z1
	double precision x2, y2, z2
	
	double precision frac, n_ch, part_nch, coll_nch

	double precision tmp, rdn

	integer id
	double precision fun_eta, fun_max

	integer i_ptc
	integer N_ptc, Nptc_rpdcut
	double precision ET_ptc, ETptc_rpdcut

	double precision mu_z

	double precision b_pt, c_pt
	double precision pt, phi, rpd

	double precision mean_pt

	N_ptc = 0
	mean_pt = 0.d0
	ET_ptc = 0.d0

	Nptc_rpdcut = 0
	ETptc_rpdcut = 0.d0

!	for momentum distribution dN/dpt=1/(1+pt/b)^c, its cdf is 1-(1+pt/b)^{1-c}=rdn, so pt = b[(1-rdn)^{1/(1-c)}-1] 	
!	b_pt = 0.83d0 
!	c_pt = 10.17d0
!
!	b_pt = 1.237d0
!	c_pt = 8.428d0

!	we use the second one
!	for momentum distribution dN/dpt^2=1/(1+pt^2/b^2)^c, its cdf is 1-(1+pt^2/b^2)^{1-c)=rdn, so pt^2 = b^2[(1-rdn)^{1/(1-c)}-1] 	

!	b_pt = 0.8d0
!	c_pt = 5.d0

!	b_pt = 0.88d0
!	c_pt = 4.d0

!	05/25/2013
!	fit pp @ 200GeV pion data, hep-ex/0304038 (fit log(dN/d2ptdy))
!	b_pt = 1.173d0
!	c_pt = 4.204d0

!	05/25/2013
!	fit AuAu @ 200GeV pion minimum data, nucl-ex/0306021, (fit log(dN/d2ptdy))
!	b_pt = 1.042d0
!	c_pt = 4.179d0

!	05/25/2013
!	fit d+Au @ 200GeV data, nucl-ex/0306021, (fit log(dN/d2ptdy))
!	b_pt = 1.176d0
!	c_pt = 4.087d0

!	05/25/2013 (http://hepdata.cedar.ac.uk/reaction)
!	fit p+Pb data, ARXIV:1210.4520, (fit log(dN/d2ptdy))
	b_pt = 1.116d0
	c_pt = 3.076d0

!	different models: this code only works for 2-component (21) model at this moment
	if (flag_glauber .eq. 11) frac = 0.d0
	if (flag_glauber .eq. 21) frac = frac_coll
	if (flag_glauber .eq. 31) frac = 1.d0

!	generate particles from participants
	do i_n1 = 1, atom_A
		x1 = x_n1(i_n1)
		y1 = y_n1(i_n1)
		z1 = z_n1(i_n1)

		if (part_n1(i_n1) .EQ. 1) then

!			pp multiplicity fluctuations (get how many particle produced)

			call mc_pp_multiplicity_dist(n_ch, part_nch, coll_nch)

!			accept or reject particles for 2-compoent particle, and assign momenta to accepted particles
			do i_ptc = 1, nint(n_ch)
!			do i_ptc = 1, nint(part_nch)

				rdn = ran2(idum)

!				condition to accept the particles
!				if (rdn .le. (1.d0-frac)/2.d0) then
				if (rdn .le. 1.d0) then

					id = 1
					call get_rpd_particle(rpd, id)

					N_ptc = N_ptc + 1
                                        if (abs(rpd) .lt. rpd_cut) Nptc_rpdcut = Nptc_rpdcut + 1

!					coordinates
					x_ptc(N_ptc) = x1
					y_ptc(N_ptc) = y1
!					z_ptc(N_ptc) = z1 / gamma_pancake
					z_ptc(N_ptc) = 0.d0
!!					one might need a Gaussian distribution for z
!!					mu_z = z_ptc(N_ptc)
!					mu_z = 0.d0
!					call mc_gaussian(z1, mu_z, sigma_pancake)
!					z_ptc(N_ptc) = z1
				
!					momentum	
100					rdn = ran2(idum)
					if (rdn .eq. 0.d0 .OR. rdn .eq. 1.d0) go to 100
!					pt = b_pt * ((1.d0-rdn)**(1.d0/(1.d0-c_pt)) - 1.d0)
					pt = sqrt(b_pt*b_pt * ((1.d0-rdn)**(1.d0/(1.d0-c_pt)) - 1.d0))
					phi = ran2(idum) * 2.d0*PI

!					write (11,*) pt, rpd

					mean_pt = mean_pt + pt
					ET_ptc = ET_ptc + pt !* cosh(rpd)
                                        if (abs(rpd) .lt. rpd_cut) ETptc_rpdcut = ETptc_rpdcut + pt

					px_ptc(N_ptc) = pt * cos(phi)
					py_ptc(N_ptc) = pt * sin(phi)
!					massless particles
					pz_ptc(N_ptc) = pt * sinh(rpd)
					rpd_ptc(N_ptc) = rpd

				end if

			end do

		end if
	end do

!	write (*,*) "particles from A", N_ptc

	do i_n2 = 1, atom_B
		x2 = x_n2(i_n2)
		y2 = y_n2(i_n2)
		z2 = z_n2(i_n2)

		if (part_n2(i_n2) .EQ. 1) then

!			pp multiplicity fluctuations (get how many particle produced)
			call mc_pp_multiplicity_dist(n_ch, part_nch, coll_nch)

!			accept or reject particles for 2-compoent particle, and assign momenta to accepted particles
			do i_ptc = 1, nint(n_ch)
!			do i_ptc = 1, nint(part_nch)

				rdn = ran2(idum)

!				condition to accept the particles
!				if (rdn .le. (1.d0-frac)/2.d0) then
				if (rdn .le. 1.d0) then

					id = -1
					call get_rpd_particle(rpd, id)

					N_ptc = N_ptc + 1
                                        if (abs(rpd) .lt. rpd_cut) Nptc_rpdcut = Nptc_rpdcut + 1

!					coordinates
					x_ptc(N_ptc) = x2
					y_ptc(N_ptc) = y2
!					z_ptc(N_ptc) = z2 / gamma_pancake
					z_ptc(N_ptc) = 0.d0
!!					one might need a Gaussian distribution for z
!!					mu_z = z_ptc(N_ptc)
!					mu_z = 0.d0
!					call mc_gaussian(z2, mu_z, sigma_pancake)
!					z_ptc(N_ptc) = z2

!					momenta	
200					rdn = ran2(idum)
					if (rdn .eq. 0.d0 .OR. rdn .eq. 1.d0) go to 200
!					pt = b_pt * ((1.d0-rdn)**(1.d0/(1.d0-c_pt)) - 1.d0)
					pt = sqrt(b_pt*b_pt * ((1.d0-rdn)**(1.d0/(1.d0-c_pt)) - 1.d0))
					phi = ran2(idum) * 2.d0*PI

!					write (11,*) pt, rpd

					mean_pt = mean_pt + pt
					ET_ptc = ET_ptc + pt !* cosh(rpd)
                                        if (abs(rpd) .lt. rpd_cut) ETptc_rpdcut = ETptc_rpdcut + pt

					px_ptc(N_ptc) = pt * cos(phi)
					py_ptc(N_ptc) = pt * sin(phi)
!					massless particles
					pz_ptc(N_ptc) = pt * sinh(rpd)
					rpd_ptc(N_ptc) = rpd

				end if

			end do

		end if
	end do

!	write (*,*) "particles from A + B", N_ptc

!	generate particle from binary collisions
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

			if (coll_n1_n2(i_n1, i_n2) .EQ. 1) then

!				pp multiplicity fluctuations (get how many particle produced)
				call mc_pp_multiplicity_dist(n_ch, part_nch, coll_nch)

!				accept or reject particles for 2-compoent particle, and assign momenta to accepted particles
				do i_ptc = 1, nint(n_ch)

					rdn = ran2(idum)

!					condition to accept the particles
!					if (rdn .le. frac) then
					if (rdn .le. frac*2.d0/(1.d0-frac)) then

						id = 0
						call get_rpd_particle(rpd, id)

						N_ptc = N_ptc + 1
                                                if (abs(rpd) .lt. rpd_cut) Nptc_rpdcut = Nptc_rpdcut + 1

!						coordinates
						x_ptc(N_ptc) = x
						y_ptc(N_ptc) = y
!						z_ptc(N_ptc) = z / gamma_pancake
						z_ptc(N_ptc) = 0.d0
!!						one might need a Gaussian distribution for z
!!						mu_z = z_ptc(N_ptc)
!						mu_z = 0.d0
!						call mc_gaussian(z, mu_z, sigma_pancake)
!						z_ptc(N_ptc) = z

!						momenta	
300						rdn = ran2(idum)
						if (rdn .eq. 0.d0 .OR. rdn .eq. 1.d0) go to 300
!						pt = b_pt * ((1.d0-rdn)**(1.d0/(1.d0-c_pt)) - 1.d0)
						pt = sqrt(b_pt*b_pt * ((1.d0-rdn)**(1.d0/(1.d0-c_pt)) - 1.d0))
						phi = ran2(idum) * 2.d0*PI

!						write (11,*) pt, rpd

						mean_pt = mean_pt + pt
						ET_ptc = ET_ptc + pt ! * cosh(rpd)
                                                if (abs(rpd) .lt. rpd_cut) ETptc_rpdcut = ETptc_rpdcut + pt

						px_ptc(N_ptc) = pt * cos(phi)
						py_ptc(N_ptc) = pt * sin(phi)
!						massless particles
						pz_ptc(N_ptc) = pt * sinh(rpd)
					        rpd_ptc(N_ptc) = rpd

					end if

				end do

			end if
		end do
	end do

!	write (*,*) "particles from A + B + C", N_ptc

	mean_pt = mean_pt / N_ptc

!	write (*,*) "mean pt for distribution:", N_ptc, mean_pt, ET_ptc

	if (N_ptc .gt. Nptc_max) then

		write (*,*) "too many particle produced:", N_ptc

                stop

	end if

	end subroutine
