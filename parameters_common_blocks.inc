!	constant parameters
	
	double precision PI	
	parameter (PI = 3.14159265358979324d0)

	double precision ran2

	integer idum
	common /idum/ idum

	integer n_evt
!	parameter (n_evt = 1000) 
!	parameter (n_evt = 100000) 
!	parameter (n_evt = 1000000) 
!	parameter (n_evt = 10000000) 
!	parameter (n_evt = 100000000) 
!	parameter (n_evt = 200000000) 
	parameter (n_evt = 1000000000) 

	integer ntrptcevt_max
!	parameter (ntrptcevt_max = 2) 
	parameter (ntrptcevt_max = 200) 
!	parameter (ntrptcevt_max = 1000) 
!	parameter (ntrptcevt_max = 10000) 
!	parameter (ntrptcevt_max = 100000) 
!	parameter (ntrptcevt_max = 1000000) 
!	parameter (ntrptcevt_max = 2000000) 
!	parameter (ntrptcevt_max = 10000000) 

!	nucleon size
	double precision nucleon_dr
	parameter (nucleon_dr = 0.2d0)

	integer atom_A, atom_B
!	parameter (atom_A = 1, atom_B = 1)
!	parameter (atom_A = 2, atom_B = 197)
!	parameter (atom_A = 1, atom_B = 208)
	parameter (atom_A = 208, atom_B = 1)
!	parameter (atom_A = 208, atom_B = 208)
!	parameter (atom_A = 197, atom_B = 197)
!	parameter (atom_A = 63, atom_B = 63)

	double precision radius
	parameter (radius = 6.62d0) ! for Pb
!	parameter (radius = 6.38d0) ! for Au
!	parameter (radius = 4.20641d0) ! for Cu

!	sigma_NN = 42mb for sqrt(s) = 200GeV
!	sigma_NN = 41mb for sqrt(s) = 130GeV
!	sigma_NN = 36mb for sqrt(s) = 62.4GeV
!	sigma_NN = 32mb for sqrt(s) = 19.6GeV

	double precision sqrt_sNN, sigma_NN
!	parameter (sqrt_sNN = 200.d0, sigma_NN = 4.2d0)
!	parameter (sqrt_sNN = 7000.d0, sigma_NN = 6.d0)
!	parameter (sqrt_sNN = 2760.d0, sigma_NN = 5.8d0)
!	parameter (sqrt_sNN = 2760.d0, sigma_NN = 6.4d0)
	parameter (sqrt_sNN = 5020.d0, sigma_NN = 6.77d0)

!	we use frac*rho_coll + (1-frac)/2*rho_part
!	if we use frac*rho_coll + (1-frac)*rho_part, we need to change the code
!	it is better to have two options
	double precision frac_coll
	parameter (frac_coll = 0.d0)
!	parameter (frac_coll = 0.11d0)
!	parameter (frac_coll = 0.13d0)
!	we have set the code to make sure every participant produces particles
!	so 2*frac/(1-frac) < 1, that is frac<1/3, which is almost always true

!	parameters for the eta distribution
	double precision eta_0, sigma_eta
!	common /eta_0/ eta_0, sigma_eta
	parameter (eta_0 = 2.5d0, sigma_eta = 1.4d0)
!	parameter (eta_0 = 2.5d0, sigma_eta = 2.0d0)

!	the _beam rapdity
	double precision E_beam, M_beam, P_beam, y_beam
	common /y_beam/ E_beam, M_beam, P_beam, y_beam

!	this cut is set as experiment
	double precision rpd_cut
	parameter (rpd_cut = 2.4d0)

!	this is determined from the f(eta) function (see mathematic file)
!	for N_ch/(dN_ch/deta) = n_rpd_unit
	double precision n_rpd_unit
!	parameter (n_rpd_unit = 1.d0)
	parameter (n_rpd_unit = rpd_cut * 2.d0)

!	This factor maybe should include in the particle production
	double precision rescale_factor
!	parameter (rescale_factor = 1.d0)
!	parameter (rescale_factor = 1.77d0)
!	parameter (rescale_factor = 2.d0)
!	parameter (rescale_factor = 5.d0)
!	parameter (rescale_factor = 10.d0)
	parameter (rescale_factor = 7.6d0)

!	this is the factor to tune the initial total particle number so that the final multiplicity can match
!	then the above rescale_factor can be absorbed into this factor
	double precision resc_fac
	parameter (resc_fac = 1.d0)
!	parameter (resc_fac = 1.77d0)

!	the parameters in negative binomial distribution for pp multiplicity
!	we fix these two parameters by setting res_fac to be 1, tuen the total mulplicity distribution to experimental data
	double precision mean_npp, k_pp
!	parameter (mean_npp = 2.35d0, k_pp = 1.9d0)
!	parameter (mean_npp = 2.34d0, k_pp = 2.d0)	! STAR: Jon Gan's thesis
!	parameter (mean_npp = 2.55d0, k_pp = 1.86d0)	! A.K. Dash & B. Mohanty: arxi:0908.0888
!	parameter (mean_npp = 2.48d0, k_pp = 2.3d0)	! UA5 p-pbar analysis
!	parameter (mean_npp = 0.32d0, k_pp = 1.88d0)	! PHENIX paper: PRC 78, 044902(2008), with limited acceptance

!	LHC: 7TeV
!	parameter (mean_npp = 6.0d0, k_pp = 1.1d0)
!	LHC: 5.02TeV
!	parameter (mean_npp = 5.36d0, k_pp = 0.9d0)
!	parameter (mean_npp = 5.36d0*resc_fac, k_pp = 0.6d0*resc_fac)
	parameter (mean_npp = 5.0d0*resc_fac, k_pp = 0.63d0*resc_fac)
!	parameter (mean_npp = 5.0d0*resc_fac, k_pp = 0.63d0)
!	parameter (mean_npp = 5.1d0*resc_fac, k_pp = 0.64d0*resc_fac)
!	LHC: 2.36TeV
!	parameter (mean_npp = 4.50d0, k_pp = 1.3d0)
!	LHC: 0.9TeV
!	parameter (mean_npp = 3.6d0, k_pp = 1.7d0)

!	the flag of Glauber model
	integer flag_glauber
	parameter (flag_glauber = 21)

!	we first get the non-contracted profile, then transform to pancake-like profile by shrinking z-direction
!	since we make all particles produce at z=0, we don't need this and the width below
	double precision gamma_pancake
!	parameter (gamma_pancake = 3500.d0)
!	parameter (gamma_pancake = 2510.d0)
!	parameter (gamma_pancake = 1380.d0)
!	parameter (gamma_pancake = 100.d0)
	parameter (gamma_pancake = 1.d0)

!	the gaussian width for shrinked z-direction
	double precision sigma_pancake
!	parameter (sigma_pancake = 0.2d0)
	parameter (sigma_pancake = 0.d0)

	double precision tau0_hydro
	parameter (tau0_hydro = 0.6d0)

	double precision t_fs
	common /t_fs/ t_fs

!	the gaussian smearing width for getting smooth/continuous profile
	double precision smear_xyz, smear_xy, smear_z, smear_eta
!	common /smear_xyz/ smear_xyz, smear_xy, smear_z, smear_eta
!	parameter(smear_xy = tau0_hydro, smear_eta = 0.7d0)
	parameter(smear_xy = 0.7d0, smear_eta = 0.7d0)
!	the value of smear_xy is fixed such that for 1000 event, the average initial conditions (in x & y directions) are pretty smooth
!	after fixing the value of smear_xy, we may fix the value of smear_eta to get the final rapidity dN/deta distribution

	integer index_pt
	parameter (index_pt = 0)

	integer npower_feta
	parameter (npower_feta = 1)

!!!!!!!!!!!!!!!!!!!!!!!

! 	commmon blocks

	double precision x_n1, y_n1, z_n1, x_n2, y_n2, z_n2
	common /xyz_part/ x_n1(atom_A), y_n1(atom_A), z_n1(atom_A), x_n2(atom_B), y_n2(atom_B), z_n2(atom_B)

	integer part_n1, part_n2
	common /if_part/ part_n1(atom_A), part_n2(atom_B)

	double precision nch_n1, nch_n2
	common /nch_part/ nch_n1(atom_A), nch_n2(atom_B)

!	for each event
!	double precision xxx_n1, yyy_n1, zzz_n1, xxx_n2, yyy_n2, zzz_n2	
!	common /participant_positions/ xxx_n1(n_evt, atom_A), yyy_n1(n_evt, atom_A), zzz_n1(n_evt, atom_A), xxx_n2(n_evt, atom_B), yyy_n2(n_evt, atom_B), zzz_n2(n_evt, atom_B)
!	integer participant_n1, participant_n2
!	common /if_participant/ participant_n1(n_evt, atom_A), participant_n2(n_evt, atom_B)

	double precision x_coll, y_coll, z_coll
	common /xyz_coll/ x_coll(atom_A, atom_B), y_coll(atom_A, atom_B), z_coll(atom_A, atom_B)

	integer coll_n1_n2
	common /if_coll/ coll_n1_n2(atom_A, atom_B)

	double precision nch_n1_n2
	common /nch_coll/ nch_n1_n2(atom_A, atom_B)

!	for each event
!	double precision xxx_coll, yyy_coll, zzz_coll
!	common /collision_positions/ xxx_coll(n_evt, atom_A, atom_B), yyy_coll(n_evt, atom_A, atom_B), zzz_coll(n_evt, atom_A, atom_B)
!	integer collision_n1_n2
!	common /if_collision/ collision_n1_n2(n_evt, atom_A, atom_B)

	double precision cdf_nb
	common /cdf_nb/ cdf_nb(0:1000)

	double precision part_cdf_nb
	common /part_cdf_nb/ part_cdf_nb(0:1000)

	double precision coll_cdf_nb
	common /coll_cdf_nb/ coll_cdf_nb(0:1000)

	integer Nptc_max
	parameter (Nptc_max = 100000)

	double precision x_ptc, y_ptc, z_ptc, eta_ptc
	common /xyz_particle/ x_ptc(Nptc_max), y_ptc(Nptc_max), z_ptc(Nptc_max), eta_ptc(Nptc_max)

	double precision px_ptc, py_ptc, pz_ptc, rpd_ptc
	common /pxpypz_particle/ px_ptc(Nptc_max), py_ptc(Nptc_max), pz_ptc(Nptc_max), rpd_ptc(Nptc_max)

!	for 3D hydro initial conditions

	double precision eps_0
	parameter (eps_0 = 0.146517d0)

	integer nx, ny, nz, neta
	parameter (nx = 81, ny = 81, nz = 81, neta=51)
!	parameter (nx = 31, ny = 31, nz = 31, neta=31)
!	parameter (nx = 71, ny = 71, nz = 71, neta=41)
	double precision dx, dy, dz, deta
	parameter (dx = 0.1d0, dy = 0.1d0, dz = 0.1d0, deta = 0.2d0)

	double precision ep_density, v_flow
	common /ep_density/ ep_density(0:3, -nx:nx, -ny:ny, -neta:neta), v_flow(1:3, -nx:nx, -ny:ny, -neta:neta)

	double precision ep_density_tot, v_flow_tot
	common /ep_density_tot/ ep_density_tot(0:3, -nx:nx, -ny:ny, -neta:neta),v_flow_tot(0:3, -nx:nx, -ny:ny, -neta:neta)

        double precision nevt_eta, nevt_eta_tot, nevt_eta_unit
	common /nevt_eta/ nevt_eta(-100:100), nevt_eta_tot, nevt_eta_unit

	double precision fun_eta_max_right, fun_eta_max_left, fun_eta_max_lrboth
	common /fun_eta_max/ fun_eta_max_right, fun_eta_max_left, fun_eta_max_lrboth

