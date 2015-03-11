!	generate the number of produced particles from one nucleon-nucleon collision
	subroutine mc_pp_multiplicity_dist(n_ch, part_nch, coll_nch)

	implicit none	

	include "parameters_common_blocks.inc"

	double precision n_ch, part_nch, coll_nch

	integer n

	double precision u

	u = ran2(idum)

	do n = 0, 500

!		cdf_nb(0:500) need to be called somewhere before
		if (u .le. cdf_nb(n)) exit

	end do

	n_ch = n

	do n = 0, 250

!		cdf_nb(0:250) need to be called somewhere before
		if (u .le. part_cdf_nb(n)) exit

	end do

	part_nch = n

	do n = 0, 250

!		cdf_nb(0:250) need to be called somewhere before
		if (u .le. coll_cdf_nb(n)) exit

	end do

	coll_nch = n

	end subroutine


!	the cdf for the NB distribution
	subroutine cdf_pp_multiplicity_dist(cdf, part_cdf, coll_cdf)

	implicit none

	double precision cdf(0:1000), part_cdf(0:1000), coll_cdf(0:1000)

	double precision pp_multiplicity_dist, n_ch
	double precision part_pp_multiplicity_dist
	double precision coll_pp_multiplicity_dist

	integer n

	double precision temp

	double precision beta, betacf, betai

	temp = 0.d0
	do n = 0, 500
		
		n_ch = n

		temp = temp + pp_multiplicity_dist(n_ch)

		cdf(n) = temp

	end do

!	for double check (using pdf of participant scaling)
	temp = 0.d0
	do n = 0, 250
		
		n_ch = n

		temp = temp + part_pp_multiplicity_dist(n_ch)

		part_cdf(n) = temp

	end do

!	for double check (using pdf of binary collision scaling)
	temp = 0.d0
	do n = 0, 250
		
		n_ch = n

		temp = temp + coll_pp_multiplicity_dist(n_ch)

		coll_cdf(n) = temp

	end do

!	write (*,*) "2:", coll_cdf(0), coll_cdf(1), coll_cdf(2), coll_cdf(3)

	end subroutine

