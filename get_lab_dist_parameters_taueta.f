
	subroutine get_lab_dist_parameters_taueta(N_ptc, x0, y0, eta0, eps_mom_lab, emt)

	implicit none

	include "parameters_common_blocks.inc"

	integer i_ptc
	integer N_ptc

	double precision fxyz_gauss
	double precision x0, y0, z0, eta0
	double precision xi, yi, zi, etai

	double precision pxi, pyi, pzi

	double precision p0i, pti, pxyzi

        double precision mti, rpdi, ptaui, petai

	double precision pp(0:3)
	double precision eps_mom_lab(0:3), eps_lab, mom_lab(3)

	double precision mtg(0:3, 0:3)
	double precision emt(0:3, 0:3) 

	integer m, n, i

	double precision tmp

!	define the metric tensor
	do m = 0, 3
		do n = 0, 3
			mtg(m, n) = 0.d0
		end do
	end do

	mtg(0,0) = 1.d0
	mtg(1,1) = -1.d0
	mtg(2,2) = -1.d0
	mtg(3,3) = -1.d0

!	get energy momentum density (flow) in lab frame
	do m = 0, 3

		tmp = 0.d0
		do i_ptc = 1, N_ptc

			xi = x_ptc(i_ptc)
			yi = y_ptc(i_ptc)
			zi = z_ptc(i_ptc)
			etai = eta_ptc(i_ptc)

			pxi = px_ptc(i_ptc)
			pyi = py_ptc(i_ptc)
			pzi = pz_ptc(i_ptc)
			rpdi = rpd_ptc(i_ptc)
			pti = sqrt(pxi*pxi + pyi*pyi)
			pxyzi = sqrt(pxi*pxi + pyi*pyi + pzi*pzi)

!			massless
			p0i = pxyzi
                        mti = pti

!                        write (*,*) i_ptc, xi, yi, zi, etai
!                        write (*,*) i_ptc, xi, yi, zi, etai, ptaui, petai
                        
			if ((x0-xi)**2 + (y0-yi)**2 .lt. 10.d0*smear_xy**2) then
!			if ((eta0-etai)**2 .lt. 10.d0*smear_eta**2) then
	
                                ptaui = mti * cosh(rpdi - etai)
                                petai = mti * sinh(rpdi - etai) / tau0_hydro

!				pp(0) = p0i
!				pp(1) = pxi
!				pp(2) = pyi
!				pp(3) = pzi

				pp(0) = ptaui
				pp(1) = pxi
				pp(2) = pyi
				pp(3) = petai

!				fxyz_gauss = 1.d0/(sqrt(2.d0*PI*smear_xyz*smear_xyz))**3 * exp(-((x0-xi)**2 + (y0-yi)**2 + (z0-zi)**2)/(2.d0*smear_xyz*smear_xyz))
				fxyz_gauss = 1.d0/(2.d0*PI*smear_xy**2) * exp(-((x0-xi)**2 + (y0-yi)**2)/(2.d0*smear_xy**2))
!				fxyz_gauss = 1.d0
				fxyz_gauss = fxyz_gauss * 1.d0/sqrt(2.d0*PI*smear_eta**2) * exp(-(eta0-etai)**2/(2.d0*smear_eta**2))

				if (ptaui .gt. 0.d0) tmp = tmp + fxyz_gauss * pp(m) / tau0_hydro 
!				if (ptaui .gt. 0.d0) tmp = tmp + fxyz_gauss 

                        end if
!			end if

		end do

		eps_mom_lab(m) = tmp

!		rescale the final transverse energy density (as we got everything from final state information) to the time at \tau=0.6fm/c
		eps_mom_lab(m) = eps_mom_lab(m) * rescale_factor

	end do

	return

	end subroutine
	
