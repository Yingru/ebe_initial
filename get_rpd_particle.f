
	subroutine get_rpd_particle(eta, id)

	implicit none

	include "parameters_common_blocks.inc"

	double precision eta
	integer id

	double precision rdn, rdn_eta, rdn_fun

	double precision fun

	double precision fun_eta

	double precision eta_max, fun_max

!	maximum values for x&y directions
	eta_max = y_beam
	fun_max = 1.d0 + (eta_0 + eta)/y_beam
!	fun_max = (1.d0 + (eta_0 + eta)/y_beam)**2

	if (id .eq. 1) fun_max = fun_eta_max_right
	if (id .eq. -1) fun_max = fun_eta_max_left
	if (id .eq. 0) fun_max = fun_eta_max_lrboth

!	get eta from the desired distribution by throwing dices (simple rejection method)
100	rdn = ran2(idum)
	rdn_eta = (rdn - 0.5d0) *  2.d0 * eta_max

!       for symmetric one, choose two asymmetric distribution with equal probability
!        if (id .eq. 0) then
!                rdn = ran2(idum)
!                if (rdn .gt. 0.5d0) id = 1
!                if (rdn .le. 0.5d0) id = -1
!        end if

	eta = rdn_eta
	fun = fun_eta(eta, id)

	rdn = ran2(idum)
	rdn_fun = rdn * fun_max

	if (rdn_fun .gt. fun) goto 100

	end

!	define the eta distribution function

	function fun_eta(eta, id)

	implicit none

	include "parameters_common_blocks.inc"

	double precision fun_eta, eta
	
	integer id

	double precision fun_eta_left, fun_eta_right, fun_eta_lrboth

!	moving in z direction (for forward moving participants)
	if (id .eq. 1) fun_eta = fun_eta_right(eta)

!	moving in -z direction (for backward moving participants)
	if (id .eq. -1) fun_eta = fun_eta_left(eta)

!	either moving in z or -z direction (for binary collisions)
	if (id .eq. 0) fun_eta = fun_eta_lrboth(eta)

	return

	end

!	get the max values of eta distribution function

	subroutine get_fun_eta_max()

	implicit none

	include "parameters_common_blocks.inc"

	double precision eta
	double precision fun_eta_right, fun_eta_left, fun_eta_lrboth

	double precision right, left, lrboth

	integer ieta, neta_here
	double precision eta_max, deta_here

	neta_here = 1000
	eta_max =  y_beam
	deta_here = eta_max / dble(neta_here)

	fun_eta_max_right = 0.d0
	fun_eta_max_left = 0.d0
	fun_eta_max_lrboth = 0.d0

	do ieta = -neta_here, neta_here, 1

		eta = ieta * deta_here

		right = fun_eta_right(eta)
		left = fun_eta_left(eta)
		lrboth = fun_eta_lrboth(eta)

		if (right .gt. fun_eta_max_right) fun_eta_max_right = right
		if (left .gt. fun_eta_max_left) fun_eta_max_left = left
		if (lrboth.gt. fun_eta_max_lrboth) fun_eta_max_lrboth = lrboth 

	end do

	end

!	define the eta distribution function (right)

	function fun_eta_right(eta)

	implicit none

	include "parameters_common_blocks.inc"

	double precision fun_eta_right, eta

	double precision unitstep

!	moving in z direction (for forward moving participants)
!	fun_eta_right = exp(-(abs(eta) - eta_0)**2 / (2.d0*sigma_eta**2) * unitstep(abs(eta) - eta_0)) * (1.d0 + eta/y_beam) * unitstep(y_beam - abs(eta))
	fun_eta_right = exp(-(abs(eta) - eta_0)**2 / (2.d0*sigma_eta**2) * unitstep(abs(eta) - eta_0)) * (1.d0 + eta/y_beam)**npower_feta * unitstep(y_beam - abs(eta))

	return

	end

!	define the eta distribution function (left)

	function fun_eta_left(eta)

	implicit none

	include "parameters_common_blocks.inc"

	double precision fun_eta_left, eta

	double precision unitstep

!	moving in -z direction (for backward moving participants)
!	fun_eta_left = exp(-(abs(eta) - eta_0)**2 / (2.d0*sigma_eta**2) * unitstep(abs(eta) - eta_0)) * (1.d0 - eta/y_beam) * unitstep(y_beam - abs(eta))
	fun_eta_left = exp(-(abs(eta) - eta_0)**2 / (2.d0*sigma_eta**2) * unitstep(abs(eta) - eta_0)) * (1.d0 - eta/y_beam)**npower_feta * unitstep(y_beam - abs(eta))

	return

	end

!	define the eta distribution function (left&right)

	function fun_eta_lrboth(eta)

	implicit none

	include "parameters_common_blocks.inc"

	double precision fun_eta_lrboth, eta

	double precision right, left
	double precision fun_eta_right, fun_eta_left

	double precision unitstep

!	either moving in z or -z direction (for binary collisions)
!	fun_eta_lrboth = exp(-(abs(eta) - eta_0)**2 / (2.d0*sigma_eta**2) * unitstep(abs(eta) - eta_0)) * unitstep(y_beam - abs(eta))
	fun_eta_lrboth = exp(-(abs(eta) - eta_0)**2 / (2.d0*sigma_eta**2) * unitstep(abs(eta) - eta_0)) * unitstep(y_beam - abs(eta))

	right = fun_eta_right(eta)
	left = fun_eta_left(eta)
	fun_eta_lrboth = 0.5d0 * (right + left)

	return

	end

!	define step function

	function unitstep(x)

	implicit none

	double precision x, unitstep

	unitstep = 0.d0
	if (x. ge. 0.d0) unitstep = 1.d0

	end 
