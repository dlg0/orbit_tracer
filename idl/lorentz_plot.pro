function dlg_interpB, r, bStruct, bMag = bMag

	bR  = interpolate ( bStruct.bR, $
			( r[0] - bStruct.rleft ) / bStruct.rdim * (bStruct.nW-1.0), $
        ( r[2] - min ( bStruct.z ) ) / bStruct.zdim * (bStruct.nH-1.0), $
		cubic = -0.5 )
    bPhi  = interpolate ( bStruct.bPhi, $
			( r[0] - bStruct.rleft ) / bStruct.rdim * (bStruct.nW-1.0), $
        ( r[2] - min ( bStruct.z ) ) / bStruct.zdim * (bStruct.nH-1.0), $
		cubic = -0.5 )
    bz  = interpolate ( bStruct.bz, $
			( r[0] - bStruct.rleft ) / bStruct.rdim * (bStruct.nW-1.0), $
        ( r[2] - min ( bStruct.z ) ) / bStruct.zdim * (bStruct.nH-1.0), $
		cubic = -0.5 )

	bMag	= sqrt ( bR^2 + bPhi^2 + bz^2 )

	return, [ bR, bPhi, bz ]

end


function dlg_vxB, r, v, bStruct

	b   = dlg_interpB ( r, bStruct, bMag = bMag )

	vxB_R	= v[1] * b[2] - v[2] * b[1]
	vxB_phi	= -1.0 * ( v[0] * b[2] - v[2] * b[0] )
	vxB_z	= v[0] * b[1] - v[1] * b[0]

	return, [ vxB_R, vxB_phi, vxB_z ]

end

function dlg_vPar, r, u, bStruct

    bDotGradB  = interpolate ( bStruct.bDotGradB, $
			( r[0] - bStruct.rleft ) / bStruct.rdim * (bStruct.nW-1.0), $
        ( r[2] - min ( bStruct.z ) ) / bStruct.zdim * (bStruct.nH-1.0), $
		cubic = -0.5 )
	
   return, -u / bStruct.mi * bDotGradB

end

function dlg_vPer, r, u, bStruct

    bR  = interpolate ( bStruct.bR, $
			( r[0] - bStruct.rleft ) / bStruct.rdim * (bStruct.nW-1.0), $
        ( r[2] - min ( bStruct.z ) ) / bStruct.zdim * (bStruct.nH-1.0), $
		cubic = -0.5 )
    bPhi  = interpolate ( bStruct.bPhi, $
			( r[0] - bStruct.rleft ) / bStruct.rdim * (bStruct.nW-1.0), $
        ( r[2] - min ( bStruct.z ) ) / bStruct.zdim * (bStruct.nH-1.0), $
		cubic = -0.5 )
    bz  = interpolate ( bStruct.bz, $
			( r[0] - bStruct.rleft ) / bStruct.rdim * (bStruct.nW-1.0), $
        ( r[2] - min ( bStruct.z ) ) / bStruct.zdim * (bStruct.nH-1.0), $
		cubic = -0.5 )
    
    bMag    = sqrt ( bR^2 + bPhi^2 + bz^2 )

    return, sqrt ( 2.0 * u * bMag / bStruct.mi )

end

function dlg_checkIn, r, bStruct

	inOut  = interpolate ( bStruct.mask, $
			( r[0] - bStruct.rleft ) / bStruct.rdim * (bStruct.nW-1.0), $
        ( r[2] - min ( bStruct.z ) ) / bStruct.zdim * (bStruct.nH-1.0) )

	if inOut gt 0 then return, 1 else $
		return, 0

end

function dlg_gc_velocity, vPer, vPar, r, bStruct 

    grad_R = interpolate ( bStruct.grad_R, $
			( r[0] - bStruct.rleft ) / bStruct.rdim * (bStruct.nW-1.0), $
        ( r[2] - min ( bStruct.z ) ) / bStruct.zdim * (bStruct.nH-1.0), $
		cubic = -0.5 ) 
    grad_phi = interpolate ( bStruct.grad_phi, $
			( r[0] - bStruct.rleft ) / bStruct.rdim * (bStruct.nW-1.0), $
        ( r[2] - min ( bStruct.z ) ) / bStruct.zdim * (bStruct.nH-1.0), $
		cubic = -0.5 ) 
    grad_z = interpolate ( bStruct.grad_z, $
			( r[0] - bStruct.rleft ) / bStruct.rdim * (bStruct.nW-1.0), $
        ( r[2] - min ( bStruct.z ) ) / bStruct.zdim * (bStruct.nH-1.0), $
		cubic = -0.5 ) 

    curv_R = interpolate ( bStruct.curv_R, $
			( r[0] - bStruct.rleft ) / bStruct.rdim * (bStruct.nW-1.0), $
        ( r[2] - min ( bStruct.z ) ) / bStruct.zdim * (bStruct.nH-1.0), $
		cubic = -0.5 ) 
    curv_phi = interpolate ( bStruct.curv_phi, $
			( r[0] - bStruct.rleft ) / bStruct.rdim * (bStruct.nW-1.0), $
        ( r[2] - min ( bStruct.z ) ) / bStruct.zdim * (bStruct.nH-1.0), $
		cubic = -0.5 ) 
    curv_z = interpolate ( bStruct.curv_z, $
			( r[0] - bStruct.rleft ) / bStruct.rdim * (bStruct.nW-1.0), $
        ( r[2] - min ( bStruct.z ) ) / bStruct.zdim * (bStruct.nH-1.0), $
		cubic = -0.5 ) 

    unitb_R = interpolate ( bStruct.unitb_R, $
			( r[0] - bStruct.rleft ) / bStruct.rdim * (bStruct.nW-1.0), $
        ( r[2] - min ( bStruct.z ) ) / bStruct.zdim * (bStruct.nH-1.0), $
		cubic = -0.5 ) 
    unitb_phi = interpolate ( bStruct.unitb_phi, $
			( r[0] - bStruct.rleft ) / bStruct.rdim * (bStruct.nW-1.0), $
        ( r[2] - min ( bStruct.z ) ) / bStruct.zdim * (bStruct.nH-1.0), $
		cubic = -0.5 ) 
    unitb_z = interpolate ( bStruct.unitb_z, $
			( r[0] - bStruct.rleft ) / bStruct.rdim * (bStruct.nW-1.0), $
        ( r[2] - min ( bStruct.z ) ) / bStruct.zdim * (bStruct.nH-1.0), $
		cubic = -0.5 ) 

    vgc_R   = vPar * unitb_R + vPer^2 * grad_R + vPar^2 * curv_R 
    vgc_phi   = vPar * unitb_phi + vPer^2 * grad_phi + vPar^2 * curv_phi
    vgc_z   = vPar * unitb_z + vPer^2 * grad_z + vPar^2 * curv_z

    return, [ vgc_R, vgc_phi, vgc_z ]

end 

pro lorentz_plot, $
	curv = curv, $
	grad = grad

	eqdsk_fName	= '../eqdsk/g120740.00275.EFIT02.mds.uncorrected.qscale_1.00000'
	eqdsk_fName	= '../eqdsk/g129x129_1051206002.01120'
	eqdsk_fName	= '../eqdsk/eqdsk.122993'
	;eqdsk_fName	= '../eqdsk/JET_75500B19_trxpl_15_plasma_state.geq'
	;eqdsk_fName	= '../eqdsk/Scen4_bn2.57_129x129'

	eqdsk = readGEqdsk ( eqdsk_fName )	

	eqdsk.bR	= eqdsk.bR*0
	eqdsk.bz	= eqdsk.bz*0
	eqdsk.bPhi	= -eqdsk.bPhi
	eqdsk.bMag	= sqrt ( eqdsk.bR^2+eqdsk.bPhi^2+eqdsk.bz^2)

	q	= 1.602e-19
	mi	= 1.672e-27 
	
	r   = 2.1
	z   = 0.0 
	phi	= 0.0
	
	vR   	= 0e6
	vz    	= 0e6
	vPhi	= 3.1e7 

	en_	= mi * ( vR^2 + vz^2 + vPhi^2 ) / 2.0 / 1.602e-19 * 1e-3; [keV]
	print, en_, 'keV'
	
	dt = 0.02e-7 
	
	rr	= [ R, phi, z ]
	vv	= [ vR, vPhi, vz ]


	;	Calculate the Guiding Center terms

	omega   = q * eqdsk.bMag / mi

	bPhi_B   = eqdsk.bPhi / eqdsk.bMag
	bR_B     = eqdsk.bR / eqdsk.bMag
	bz_B     = eqdsk.bz / eqdsk.bMag

	lnB = alog ( eqdsk.bMag )
	lnB_dR  = dlg_pDeriv ( lnB, 1, eqdsk.rStep )
	lnB_dz  = dlg_pDeriv ( lnB, 2, eqdsk.zStep )

	bR_B_dz  	= dlg_pDeriv ( bR_B, 2, eqdsk.zStep )
	bPhi_B_dz  	= dlg_pDeriv ( bPhi_B, 2, eqdsk.zStep )
	bz_B_dR  	= dlg_pDeriv ( bz_B, 1, eqdsk.rStep )
	bPhi_B_dR  	= dlg_pDeriv ( bPhi_B, 1, eqdsk.rStep )
	bR_B_dR  	= dlg_pDeriv ( bR_B, 1, eqdsk.rStep )
	bz_B_dz  	= dlg_pDeriv ( bz_B, 2, eqdsk.zStep )

	gradB_R 	= dlg_pDeriv ( eqdsk.bMag, 1, eqdsk.rStep )
	gradB_z 	= dlg_pDeriv ( eqdsk.bMag, 2, eqdsk.zStep )

	bDotGradB   = bR_B * gradB_R + bz_B * gradB_z

	B_dR  = dlg_pDeriv ( eqdsk.bMag, 1, eqdsk.rStep )
	B_dz  = dlg_pDeriv ( eqdsk.bMag, 2, eqdsk.zStep )
	
	;	Gradient Drift Terms
	;	--------------------
	;	Pick from either of these gradient term
	;	formulations, they both work.

	if not keyword_set ( grad ) then grad = 1

	if grad eq 1 then begin

		bxGradLnB_R = 1.0 / ( 2.0 * Omega ) * ( bPhi_B * lnB_dz )
		bxGradLnB_phi   = -1.0 / ( 2.0 * Omega ) * ( bR_B * lnB_dz - bz_B * lnB_dR )
		bxGradLnB_z = -1.0 / ( 2.0 * Omega ) * ( bPhi_B * lnB_dR )

		grad_R		= bxGradLnB_R
		grad_phi	= bxGradLnB_phi
		grad_z		= bxGradLnB_z

	endif else begin
	
		BxGradB_B2_R    = eqdsk.bPhi * B_dz / ( 2.0 * omega * eqdsk.bMag^2 )
		BxGradB_B2_phi  = -1.0 * (eqdsk.bR * B_dz - eqdsk.bz * B_dR ) / ( 2.0 * omega * eqdsk.bMag^2 )
		BxGradB_B2_z    = -1.0 * eqdsk.bPhi * B_dR / ( 2.0 * omega * eqdsk.bMag^2 )
	
		grad_R		= BxGradB_B2_R
		grad_phi	= BxGradB_B2_phi
		grad_z		= BxGradB_B2_z

	endelse

	;	Curvature Drift Terms
	;	---------------------
	;	Select from any of the three formulations, they all
	;	work, i had to code them all to prove i had the right
	;	answer dammit ;-p

	if not keyword_set ( curv ) then curv = 1

	if curv eq 1 then begin

		bDotGradB_R     = bR_B * bR_B_dR + bz_B * bR_B_dz $
					- bPhi_B^2 / rebin(eqdsk.r,eqdsk.nW,eqdsk.nH)  
		bDotGradB_phi   = bPhi_B * bR_B / rebin(eqdsk.r,eqdsk.nW,eqdsk.nH) $
					+ bR_B * bPhi_B_dR + bz_B * bPhi_B_dz
		bDotGradB_z     = br_B * bz_B_dR + bz_B * bz_B_dz

		bxbDotGradB_R   = 1.0 / Omega * ( bPhi_B * bDotGradB_z - bz_B * bDotGradB_phi )
		bxbDotGradB_phi = -1.0 / Omega * ( br_B * bDotGradB_z - bz_B * bDotGradB_R )
		bxbDotGradB_z   = 1.0 / Omega * ( br_B * bDotGradB_phi - bPhi_B * bDotGradB_R ) 

		curv_R		= bxbDotGradB_R
		curv_phi	= bxbDotGradB_phi
		curv_z		= bxbDotGradB_z

	endif else if curv eq 2 then begin
	
		KR  = -1.0 * ( bPhi_B * ( bPhi_B_dR + 1.0 / rebin ( eqdsk.r, eqdsk.nW, eqdsk.nH ) * bPhi_B ) $
				- bz_B * ( bR_B_dz - bz_B_dR ) )
		KPhi    = bR_B * ( bPhi_B_dR + 1.0 / rebin ( eqdsk.r, eqdsk.nW, eqdsk.nH ) * bPhi_B ) $
				+ bz_B * bPhi_B_dz
		Kz  = -1.0 * ( bR_B * ( bR_B_dz - bz_B_dR ) + bPhi_B * bPhi_B_dz )
		
		bxK_R   = 1.0 / Omega * ( bPhi_B * Kz - bz_B * KPhi )
		bxK_phi = -1.0 / Omega * ( bR_B * Kz - bz_B * KR )
		bxK_z   = 1.0 / Omega * ( bR_B * KPhi - bPhi_B * KR )
	
		curv_R		= bxK_R
		curv_phi	= bxK_phi
		curv_z		= bxK_z

	endif else begin
	
		gradPerpB_R = B_dR - br_B^2 * B_dR - bz_B * br_B * B_dz
		gradPerpB_phi   = -br_B * bPhi_B * B_dR - bz_B * bPhi_B * B_dz
		gradPerpB_z = B_dz - br_B * bz_B * B_dR - bz_B^2 * B_dz

		KR_   = gradPerpB_R / eqdsk.bMag
		KPhi_   = gradPerpB_phi / eqdsk.bMag
		Kz_   = gradPerpB_z / eqdsk.bMag

		bxK_R_   = 1.0 / Omega * ( bPhi_B * Kz_ - bz_B * KPhi_ )
		bxK_phi_ = -1.0 / Omega * ( bR_B * Kz_ - bz_B * KR_ )
		bxK_z_   = 1.0 / Omega * ( bR_B * KPhi_ - bPhi_B * KR_ )

		curv_R		= bxK_R_
		curv_phi	= bxK_phi_
		curv_z		= bxK_z_

	endelse

	bStruct = { unitb_R : bR_B, $
            unitb_phi : bPhi_B, $
            unitb_z : bz_B, $
            grad_R : grad_R, $
            grad_phi : grad_phi, $
            grad_z : grad_z, $
            curv_R : curv_R, $
            curv_phi : curv_phi, $
            curv_z : curv_z, $
            nW : eqdsk.nW, $
            nH : eqdsk.nH, $
            rDim : eqdsk.rdim, $
            rLeft : eqdsk.rleft, $
            zDim : eqdsk.zdim, $
            bR : eqdsk.bR, $
            bPhi : eqdsk.bPhi, $
            bz : eqdsk.bz, $
            z : eqdsk.z, $
            R : eqdsk.r, $
            bDotGradB : bDotGradB, $
			mask : eqdsk.mask, $
			mi : mi }


	vv_array	= transpose ( vv )
	rr_array	= transpose ( rr )
	tt_array	= 0

	nSteps = 5000
	for i = 0, nSteps - 2 do begin
	
		vxB	= dlg_vxB ( rr, vv, bStruct )
		k1_v	= dt * q / mi * vxB	
		k1_r	= dt * ( vv )
	
		vxB	= dlg_vxB ( rr + k1_r / 2.0, vv + k1_v / 2.0, bStruct )
		k2_v	= dt * q / mi * vxB	
		k2_r	= dt * ( vv + k1_v / 2.0 )
	
		vxB	= dlg_vxB ( rr + k2_r / 2.0, vv + k2_v / 2.0, bStruct )
		k3_v	= dt * q / mi * vxB	
		k3_r	= dt * ( vv + k2_v / 2.0 )
	
		vxB	= dlg_vxB ( rr + k3_r, vv + k3_v, bStruct )
		k4_v	= dt * q / mi * vxB	
		k4_r	= dt * ( vv + k3_v )
	
		vv	= vv + ( k1_v + 2.0 * k2_v + 2.0 * k3_v + k4_v ) / 6.0
		rr	= rr + ( k1_r + 2.0 * k2_r + 2.0 * k3_r + k4_r ) / 6.0
	
		vv_array	= [ vv_array, transpose ( vv ) ]
		rr_array	= [ rr_array, transpose ( rr ) ]
		tt_array	= [ tt_array, i * dt ]
	
	endfor


	;	Find mean guiding center start point from 
	;	approx. one gyration

	;	Use magnetic field here for a guess at the 
	;	gyro frequency, that's all
	
	pos	= [ rr[0], rr[1], rr[2] ]
	bHere_g   = dlg_interpB ( pos, bStruct, bMag = g_bMag )

	omega0	= abs ( 2.0 * q * g_bMag / mi)
	iiFirstGyration	= where ( tt_array lt 4.0 * !pi / omega0 )

	r_gc	= mean ( rr_array[iiFirstGyration,0] )
	phi_gc	= mean ( rr_array[iiFirstGyration,1] )
	z_gc	= mean ( rr_array[iiFirstGyration,2] )

	phi_gc	= 0

	pos	= [ R_gc, phi_gc, z_gc ]

	bHere_gc   = dlg_interpB ( pos, bStruct, bMag = gc_bMag )

	gc_bR_B	= bHere_gc[0] / gc_bMag
	gc_bPhi_B	= bHere_gc[1] / gc_bMag
	gc_bz_B	= bHere_gc[2] / gc_bMag

	vMag_all	= sqrt ( vv_array[iiFirstGyration,0]^2 $
					+ vv_array[iiFirstGyration,1]^2 $
					+ vv_array[iiFirstGyration,2]^2 )
	vPar_all	= ( vv_array[iiFirstGyration,0] * gc_bR_B $
					+ vv_array[iiFirstGyration,1]* gc_bPhi_B $
					+ vv_array[iiFirstGyration,2] * gc_bz_B )
	vPer_all	= sqrt ( vMag_all^2 - vPar_all^2 )

	vPar		= mean ( vPar_all ) 
	vPer		= mean ( vPer_all )  
	
	print, sqrt ( vPer^2 + vPar^2 ), sqrt ( vR^2 + vPhi^2 + vz^2 )

	;plots, rr_array[iiFirstGyration,0], rr_array[iiFirstGyration,2], color = 200, thick = 3.0
	;plots, [r,r], [z-0.6,z+0.6], lineStyle = 1
	;plots, [r-0.3,r+0.3], [z,z], lineStyle = 1

	; Trace field line for comparison
	
	fl_pos	= [ R_gc, phi_gc, z_gc ]

	fl_rArray	= fl_pos[0]
	fl_phiArray	= fl_pos[1]
	fl_zArray	= fl_pos[2] 
	
	dPhi    = -2 * !pi / 100.0
	for fl_i = 0, 250 do begin
		
		bHere   = dlg_interpB ( fl_pos, bStruct, bMag = bMag )
		K1  = dPhi * bHere / bMag
		
		bHere   = dlg_interpB ( fl_pos + K1 / 2.0, bStruct, bMag = bMag )
		K2    = dPhi * bHere / bMag
		
		bHere   = dlg_interpB ( fl_pos + K2 / 2.0, bStruct, bMag = bMag )
		K3    = dPhi * bHere / bMag
		
		bHere   = dlg_interpB ( fl_pos + K3, bStruct, bMag = bMag )
		K4    = dPhi * bHere / bMag
		
		fl_pos    = fl_pos + ( K1 + 2.0 * K2 + 2.0 * K3 + K4 ) / 6.0
	
		fl_rArray  = [ fl_rArray, fl_pos[0] ]
		fl_phiArray	= [ fl_phiArray, fl_pos[1] ]
		fl_zArray  = [ fl_zArray, fl_pos[2] ]
	
	endFor

	bHere   = dlg_interpB ( pos, bStruct, bMag = bMag )
	u   = mi * vPer^2 / ( 2.0 * bMag )

	rTrack  	= pos[0]
	phiTrack 	= pos[1]
	zTrack  	= pos[2] 
	vPerTrack  	= vPer
	vParTrack   = vPar
	
	for i = 0, nSteps - 2 do begin
	
	    vPer   = dlg_vPer ( pos, u, bStruct ) 
	    vgc = dlg_gc_velocity ( vPer, vPar, pos, bStruct )
	    k1_vPar   = dt * dlg_vPar ( pos, u, bStruct ) 
	    k1_vgc  = dt * vgc
	
	    vPer   = dlg_vPer ( pos + k1_vgc / 2.0, u, bStruct ) 
	    vgc = dlg_gc_velocity ( vPer, vPar + k1_vPar / 2.0, pos + k1_vgc / 2.0, bStruct )
	    k2_vPar   = dt * dlg_vPar ( pos + k1_vgc / 2.0, u, bStruct ) 
	    k2_vgc  = dt * vgc
	 
	    vPer   = dlg_vPer ( pos + k2_vgc / 2.0, u, bStruct ) 
	    vgc = dlg_gc_velocity ( vPer, vPar + k2_vPar / 2.0, pos + k2_vgc / 2.0, bStruct )
	    k3_vPar   = dt * dlg_vPar ( pos + k2_vgc / 2.0, u, bStruct ) 
	    k3_vgc  = dt * vgc
	
	    vPer   = dlg_vPer ( pos + k3_vgc, u, bStruct ) 
	    vgc = dlg_gc_velocity ( vPer, vPar + k3_vPar, pos + k3_vgc, bStruct )
	    k4_vPar   = dt * dlg_vPar ( pos + k3_vgc, u, bStruct ) 
	    k4_vgc  = dt * vgc
	    
	    vPar    = vPar + ( k1_vPar + 2.0 * k2_vPar + 2.0 * k3_vPar + k4_vPar ) / 6.0
	    pos   = pos + ( k1_vgc + 2.0 * k2_vgc + 2.0 * k3_vgc + k4_vgc ) / 6.0
	
		stillIn	= dlg_checkIn ( pos, bStruct )
	    
		rTrack  = [ rTrack, pos[0] ]
		phiTrack	= [ phiTrack, pos[1] ]
	    zTrack  = [ zTrack, pos[2] ]

	    vPerTrack  = [ vPerTrack, vPer ]
	    vParTrack   = [ vParTrack, vPar ]
	
	endfor

	;	convert cylindrical coords to cartesian for plotting

	x_lorentz	= rr_array[*,0] * cos ( rr_array[*,1] )
	y_lorentz	= rr_array[*,0] * sin ( rr_array[*,1] )
	z_lorentz	= rr_array[*,2]
	
	x_gc	= rTrack * cos ( phiTrack )
	y_gc	= rTrack * sin ( phiTrack )
	z_gc	= zTrack
	
	x_fl	= fl_rArray * cos ( fl_phiArray )
	y_fl	= fl_rArray * sin ( fl_phiArray )
	z_fl	= fl_zArray

   	loadct, 12, /sil, rgb_table = ct12
	red	= transpose ( ct12[12*16-1,*] )
	blue	= transpose ( ct12[8*16-1,*] )
	green	= transpose ( ct12[1*16-1,*] )

	iPlot, x_lorentz, y_lorentz, z_lorentz, $
		rgb_table = 1, $
		vert = fIndGen(nSteps)/nSteps*255, $
		/zoom_on_resize, $
		thick = 3, /iso
	;iPlot, x_gc, y_gc, z_gc, /over, $
	;	rgb_table = 3, $
	;	vert = fIndGen(nSteps)/nSteps*255, $
	;	thick = 3
	iPlot, x_fl, y_fl, z_fl, /over, $
		color = green, $
		transp = 50, $
		thick = 4

	for i = 0, 11 do begin

		bbbs_phi	= eqdsk.rbbbs * 0 + i * !pi * 2 / 12
		x_bbbs	= eqdsk.rbbbs * cos ( bbbs_phi )
		y_bbbs	= eqdsk.rbbbs * sin ( bbbs_phi )
		z_bbbs	= eqdsk.zbbbs

		iPlot, x_bbbs, y_bbbs, z_bbbs, $
			/over, $
			trans = 90, $
			thick = 6
	endfor

	;for i = 0, 60 do begin

	;	lim_phi	= eqdsk.rlim * 0 + i * !pi * 2 / 180
	;	x_lim	= eqdsk.rlim * cos ( lim_phi )
	;	y_lim	= eqdsk.rlim * sin ( lim_phi )
	;	z_lim	= eqdsk.zlim

	;	iPlot, x_lim, y_lim, z_lim, $
	;		/over, $
	;		trans = 50, $
	;		thick = 10 
	;endfor



	iPlot, fl_rArray, fl_zArray, /iso, $
		color = blue, $
		transp = 70, $
		thick = 3, $
		/zoom_on_resize
	iPlot, rTrack, zTrack, /over, $
		color = red, $
		thick = 2
	iPlot, rr_array[*,0], rr_array[*,2], /over, $
		color = green, $
		thick = 2, $
		trans = 50
	iPlot, eqdsk.rbbbs, eqdsk.zbbbs, /over
	iPlot, eqdsk.rLim, eqdsk.zLim, /over


stop
end
