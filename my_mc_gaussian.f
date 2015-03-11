!	generate two independent random numbers according to gaussian distribution (polar box-mueller transformation)
	
	subroutine mc_gaussian(x, mean, sigma)

	implicit none

	include "parameters_common_blocks.inc"

	double precision x, y, mean, sigma
	double precision rsq, w

	save y

	integer i
	data i /0/ 
	save i

	if (sigma .eq. 0.d0) then
		x = mean
		return
	end if

	if (i .eq. 1) then

		x = y
		i = 0

!		go to 100

	else

100		x = 2.d0*ran2(idum) - 1.d0
		y = 2.d0*ran2(idum) - 1.d0

		rsq = x*x + y*y

		if (rsq .GT. 1.d0) go to 100
		
		w = sqrt((-2.d0*log(rsq)) / rsq)	
		x = x * w
		y = y * w

		i = 1

	end if
		
! 	transform to non-zero mean and non-1 width, note here "+ mean"
	x = x * sigma + mean
	
!	write (*,*) i, x

	end subroutine
