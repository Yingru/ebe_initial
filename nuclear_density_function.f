
	function rho_A(r, atom)	! nuclear density function rho_A(r)
	
	implicit none

	include "parameters_common_blocks.inc"
	
	double precision rho_A, r
	integer atom
	
	double precision rho_0, R_A, omega, diff
	
	rho_A = 0D0

!	for Cu
	if (atom .eq. 63 .OR. atom .eq. 64) then 
		rho_0 = 0.1693D0
		R_A = 4.20641D0
		omega = 0D0
		diff = 0.5977D0

		rho_A = rho_0 / (1 + exp((r - R_A) / diff))
	end if

!	for Au
	if (atom .eq. 197) then  
		rho_0 = 0.1693D0
		R_A = 6.38D0
		omega = 0D0
!		diff = 0.535D0
		diff = 0.54D0

		rho_A = rho_0 / (1 + exp((r - R_A) / diff))
	end if

!	for Pb
	if (atom .eq. 207 .OR. atom .eq. 208) then  
		rho_0 = 0.16D0
		R_A = 6.62D0
		omega = 0D0
		diff = 0.546D0

		rho_A = rho_0 / (1 + exp((r - R_A) / diff))
	end if

	end function



