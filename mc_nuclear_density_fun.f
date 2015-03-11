
	subroutine mc_rho_A(x, y, z, atom)	! generate nucleon position (x, y , z) according to rho_A(r) distribution (nuclear density function)
	
	implicit none

	include "parameters_common_blocks.inc"
	
	double precision x, y, z
	double precision cos_theta, sin_theta, phi
	integer atom

	double precision rho_A, r
	
	double precision rmax_ws

	double precision rdn1, rdn2, rdn3, rdn4
	double precision accept_prob 

!	for proton
	if (atom .eq. 1) then

		x = 0.d0
		y = 0.d0
		z = 0.d0
		
		return

	end if


!	for deutron
	if (atom .eq. 2) then

		call mc_deutron(x, y, z)

		return

	end if

!	for others

	r = 0.d0
	rmax_ws = 2.d0 * radius

!100	rdn1 = ran2(idum)
!	r = rmax_ws * rdn1
!	rdn2 = ran2(idum)
!	accept_prob =  r**2 * rho_A(r, atom) / 5.d0
!
!	if (rdn2 .GT. accept_prob) go to 100

100	rdn1 = ran2(idum)
	r = rmax_ws * rdn1**(1.d0/3.d0)
	rdn2 = ran2(idum)
	accept_prob =  rho_A(r, atom) / rho_A(0.d0, atom)

	if (rdn2 .GT. accept_prob) go to 100

!	do 
!		rdn1 = ran2(idum)
!		r = rmax_ws * rdn1**(1.d0/3.d0)
!		rdn2 = ran2(idum)
!		accept_prob =  rho_A(r, atom) / rho_A(0.d0, atom)
!		if (rdn2 .LE. accept_prob) exit
!	end do

!	rdn2 = 1.d0
!	accept_prob = 0.d0
!
!	do while (rdn2 .GT. accept_prob) 
!		rdn1 = ran2(idum)
!		r = rmax_ws * rdn1**(1.d0/3.d0)
!		rdn2 = ran2(idum)
!		accept_prob =  rho_A(r, atom) / rho_A(0.d0, atom)
!	end do

	rdn3 = ran2(idum)
	cos_theta = 1.d0 - 2.d0 * rdn3
	sin_theta = sqrt(1.d0 - cos_theta**2)

	rdn4 = ran2(idum)
	phi = 2.d0 * PI * rdn4

	x = r * sin_theta * cos(phi)
	y = r * sin_theta * sin(phi)
	z = r * cos_theta

!	write (100, *) r, r
	
	end subroutine


!	the subroutine for deutron

	subroutine mc_deutron(x, y, z)

	implicit none

	include "parameters_common_blocks.inc"

	double precision x, y, z

	double precision a, b
	double precision r, H, P
	double precision rdn, rnd2, rdn3, rdn4
	double precision cos_theta, sin_theta, phi

	double precision x1, y1, z1
	double precision x2, y2, z2

	save x2, y2, z2

	integer i
	data i /0/ 
	save i

	a = 0.228d0
	b = 1.18d0

	if (i .eq. 1) then

		x = x2
		y = y2
		z = z2

		i = 0

	else

!		the maximum for deutron r is around 15fm, r/2 = 7.5fm
100		r = ran2(idum) * 15.d0
		if (r .eq. 0.d0) goto 100
!		H = (exp(-a*r) - exp(-b*r))/r
!		P =  r**2 * H**2
		P = (exp(-a*r) - exp(-b*r))**2

!		the maximum value for P is <0.3
		rdn = ran2(idum) * 0.3d0

		if (rdn .gt. P) goto 100
	
		rdn3 = ran2(idum)
		cos_theta = 1.d0 - 2.d0 * rdn3
		sin_theta = sqrt(1.d0 - cos_theta**2)
		
		rdn4 = ran2(idum)
		phi = 2.d0 * PI * rdn4
		
		x1 = 0.5d0 * r * sin_theta * cos(phi)
		y1 = 0.5d0 * r * sin_theta * sin(phi)
		z1 = 0.5d0 * r * cos_theta

		x2 = -x1
		y2 = -y1
		z2 = -z1

		x = x1
		y = y1
		z = z1

!		this sentence determine wheter we use the saved value next time
		i = 1

	end if

	end
