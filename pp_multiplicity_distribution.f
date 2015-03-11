!	The charged/total particle multiplicity distribution is well described by negative binormial (NB) distribution: mean_n and k are two parameters
	function pp_multiplicity_dist(n_ch) 

	implicit none

	include "parameters_common_blocks.inc"

	double precision pp_multiplicity_dist, n_ch

	double precision NB_dist, mean_n, k

	double precision gammln

	mean_n = mean_npp
	k = k_pp

	pp_multiplicity_dist = exp(gammln(n_ch+k)) / exp(gammln(n_ch + 1.d0)) / exp(gammln(k)) * (mean_n/k)**n_ch / (mean_n/k+1.d0)**(n_ch+k)
!	pp_multiplicity_dist = NB_dist(n_ch, mean_n, k)

!	write (*,*) mean_n, k, pp_multiplicity_dist

	end function


!	The charged/total particle multiplicity distribution is well described by negative binormial (NB) distribution: mean_n and k are two parameters
!	participant scaling
	function part_pp_multiplicity_dist(n_ch) 

	implicit none

	include "parameters_common_blocks.inc"

	double precision part_pp_multiplicity_dist, n_ch

	double precision pp_multiplicity_dist, NB_dist, mean_n, k

	double precision gammln

	mean_n = mean_npp * (1.d0 - frac_coll) / 2.d0
	k = k_pp * (1.d0 - frac_coll) / 2.d0

!	constructed directly from NB distribution
	part_pp_multiplicity_dist = exp(gammln(n_ch+k)) / exp(gammln(n_ch + 1.d0)) / exp(gammln(k)) * (mean_n/k)**n_ch / (mean_n/k+1.d0)**(n_ch+k)
!	part_pp_multiplicity_dist = NB_dist(n_ch, mean_n, k)
	return

!	constructed from pdf of binary scaling (the odds are gone, redistributed to evens)
	if (n_ch .eq. 0) then
		part_pp_multiplicity_dist = NB_dist(2.d0*n_ch, mean_n, k) + NB_dist(2.d0*n_ch+1.d0, mean_n, k)/2.d0
	else
		part_pp_multiplicity_dist = NB_dist(2.d0*n_ch-1.d0, mean_n, k)/2.d0 + NB_dist(2.d0*n_ch, mean_n, k) + NB_dist(2.d0*n_ch+1.d0, mean_n, k)/2.d0
	end if

!	write (*,*) mean_n, k, pp_multiplicity_dist

	end function


!	The charged/total collicle multiplicity distribution is well described by negative binormial (NB) distribution: mean_n and k are two parameters
!	binary collision scaling
	function coll_pp_multiplicity_dist(n_ch) 

	implicit none

	include "parameters_common_blocks.inc"

	double precision coll_pp_multiplicity_dist, n_ch

	double precision pp_multiplicity_dist, NB_dist, mean_n, k

	double precision gammln

	mean_n = mean_npp * frac_coll
	k = k_pp * frac_coll

!	constructed directly from NB distribution
	coll_pp_multiplicity_dist = exp(gammln(n_ch+k)) / exp(gammln(n_ch + 1.d0)) / exp(gammln(k)) * (mean_n/k)**n_ch / (mean_n/k+1.d0)**(n_ch+k)
!	coll_pp_multiplicity_dist = NB_dist(n_ch, mean_n, k)
	return

!	constructed from pdf of binary scaling (the odds are gone, redistributed to evens)
	if (n_ch .eq. 0) then
		coll_pp_multiplicity_dist = NB_dist(2.d0*n_ch, mean_n, k) + NB_dist(2.d0*n_ch+1.d0, mean_n, k)/2.d0
	else
		coll_pp_multiplicity_dist = NB_dist(2.d0*n_ch-1.d0, mean_n, k)/2.d0 + NB_dist(2.d0*n_ch, mean_n, k) + NB_dist(2.d0*n_ch+1.d0, mean_n, k)/2.d0
	end if

!	write (*,*) mean_n, k, pp_multiplicity_dist

	end function


	function NB_dist(n_ch, mean_n, k)

	implicit none

	double precision NB_dist, n_ch, mean_n, k

	double precision gammln

	NB_dist = exp(gammln(n_ch+k)) / exp(gammln(n_ch + 1.d0)) / exp(gammln(k)) * (mean_n/k)**n_ch / (mean_n/k+1.d0)**(n_ch+k)

	end function

