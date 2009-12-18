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

	vxB_R	= v[1] * bz - v[2] * bPhi
	vxB_phi	= -1.0 * ( v[0] * bz - v[2] * bR )
	vxB_z	= v[0] * bPhi - v[1] * br

	return, [ vxB_R, vxB_phi, vxB_z ]

end

function dlg_vPar, r, u, bStruct

    bDotGradB  = interpolate ( bStruct.bDotGradB, $
			( r[0] - bStruct.rleft ) / bStruct.rdim * (bStruct.nW-1.0), $
        ( r[2] - min ( bStruct.z ) ) / bStruct.zdim * (bStruct.nH-1.0), $
		cubic = -0.5 )
	
   return, -u / bStruct.mi * bDotGradB

end

function dlg_vPerp, r, u, bStruct

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

function dlg_gc_velocity, vPerp, vPar, r, bStruct 

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

    vgc_R   = vPar * unitb_R + vPerp^2 * grad_R + vPar^2 * curv_R 
    vgc_phi   = vPar * unitb_phi + vPerp^2 * grad_phi + vPar^2 * curv_phi
    vgc_z   = vPar * unitb_z + vPerp^2 * grad_z + vPar^2 * curv_z

    return, [ vgc_R, vgc_phi, vgc_z ]

end 

pro lorentz_plot

	eqdsk_fName	= 'g120740.00275.EFIT02.mds.uncorrected.qscale_1.00000'
	eqdsk = readGEqdsk ( eqdsk_fName )	

	q	= 1.602e-19
	mi	= 1.672e-27 
	
	r   = 1.2
	z   = 0.0 
	phi	= 0.0
	
	vR   	= 8e6
	vz    	= 0e5
	vPhi	= 0.1e6 
	
	en_	= mi * ( vR^2 + vz^2 + vPhi^2 ) / 2.0 / 1.602e-19 * 1e-3; [keV]
	
	dt = 1e-8 
	
	rr	= [ R, phi, z ]
	vv	= [ vR, vPhi, vz ]


	;	generate all the GC field terms

	Omega   = q * eqdsk.bMag / mi

	bPhi_B   = eqdsk.bPhi / eqdsk.bMag
	bR_B     = eqdsk.bR / eqdsk.bMag
	bz_B     = eqdsk.bz / eqdsk.bMag

	lnB = alog ( eqdsk.bMag )
	lnB_dR  = dlg_pDeriv ( lnB, 1, eqdsk.rStep )
	lnB_dz  = dlg_pDeriv ( lnB, 2, eqdsk.zStep )

	bxGradLnB_R = 1.0 / ( 2.0 * Omega ) * ( bPhi_B * lnB_dz )
	bxGradLnB_phi   = -1.0 / ( 2.0 * Omega ) * ( bR_B * lnB_dz - bz_B * lnB_dR )
	bxGradLnB_z = -1.0 / ( 2.0 * Omega ) * ( bPhi_B * lnB_dR )

	bR_B_dz  = dlg_pDeriv ( bR_B, 2, eqdsk.zStep )
	bPhi_B_dz    = dlg_pDeriv ( bPhi_B, 2, eqdsk.zStep )
	bz_B_dR  = dlg_pDeriv ( bz_B, 1, eqdsk.rStep )
	bPhi_B_dR    = dlg_pDeriv ( bPhi_B, 1, eqdsk.rStep )
	bR_B_dR  = dlg_pDeriv ( bR_B, 1, eqdsk.rStep )
	bz_B_dz  = dlg_pDeriv ( bz_B, 2, eqdsk.zStep )

	gradB_R = dlg_pDeriv ( eqdsk.bMag, 1, eqdsk.rStep )
	gradB_z = dlg_pDeriv ( eqdsk.bMag, 2, eqdsk.zStep )

	bDotGradB_R     = bR_B * bR_B_dR + bz_B * bR_B_dz $
				- bPhi_B^2 / rebin(eqdsk.r,eqdsk.nW,eqdsk.nH)  
	bDotGradB_phi   = bPhi_B * bR_B / rebin(eqdsk.r,eqdsk.nW,eqdsk.nH) $
				+ bR_B * bPhi_B_dR + bz_B * bPhi_B_dz
	bDotGradB_z     = br_B * bz_B_dR + bz_B * bz_B_dz

	bxbDotGradB_R   = 1.0 / Omega * ( bPhi_B * bDotGradB_z - bz_B * bDotGradB_phi )
	bxbDotGradB_phi = -1.0 / Omega * ( br_B * bDotGradB_z - bz_B * bDotGradB_R )
	bxbDotGradB_z   = 1.0 / Omega * ( br_B * bDotGradB_phi - bPhi_B * bDotGradB_R ) 

	bDotGradB   = bR_B * gradB_R + bz_B * gradB_z

	bStruct = { unitb_R : bR_B, $
            unitb_phi : bPhi_B, $
            unitb_z : bz_B, $
            grad_R : bxGradLnB_R, $
            grad_phi : bxGradLnB_phi, $
            grad_z : bxGradLnB_z, $
            curv_R : bxbDotGradB_R, $
            curv_phi : bxbDotGradB_phi, $
            curv_z : bxbDotGradB_z, $
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
			mask : eqdsk.mask }


	;	Use magnetic field here as a guess at the 
	;	gyro frequency, that's all
	
	g_bR  = interpolate ( bStruct.bR, $
			( rr[0] - bStruct.rleft ) / bStruct.rdim * (bStruct.nW-1.0), $
	    ( rr[2] - min ( bStruct.z ) ) / bStruct.zdim * (bStruct.nH-1.0), $
		cubic = -0.5 )
	g_bPhi  = interpolate ( bStruct.bPhi, $
			( rr[0] - bStruct.rleft ) / bStruct.rdim * (bStruct.nW-1.0), $
	    ( rr[2] - min ( bStruct.z ) ) / bStruct.zdim * (bStruct.nH-1.0), $
		cubic = -0.5 )
	g_bz  = interpolate ( bStruct.bz, $
			( rr[0] - bStruct.rleft ) / bStruct.rdim * (bStruct.nW-1.0), $
	    ( rr[2] - min ( bStruct.z ) ) / bStruct.zdim * (bStruct.nH-1.0), $
		cubic = -0.5 )
	
	g_bMag	= sqrt ( g_bR^2 + g_bPhi^2 + g_bz^2 )

	vv_array	= transpose ( vv )
	rr_array	= transpose ( rr )
	tt_array	= 0

	tp_nSteps = 1200
	for tp_i = 0, tp_nSteps - 1 do begin
	
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
		tt_array	= [ tt_array, tp_i * dt ]
	
	
	endfor

	plots, rr_array[*,0], rr_array[*,2]
stop
;	Find mean guiding center start point from 
;	approx. one gyration

omega0	= abs(2.0 * q * g_bMag / mi)
iiFirstGyration	= where ( tt_array lt 4.0 * !pi / omega0 )
r	= mean ( rr_array[iiFirstGyration,0] )
z	= mean ( rr_array[iiFirstGyration,2] )

gc_bR  = interpolate ( bStruct.bR, ( rr_array[iiFirstGyration,0] - bStruct.rleft ) / bStruct.rdim * (bStruct.nW-1.0), $
    ( rr_array[iiFirstGyration,2] - min ( bStruct.z ) ) / bStruct.zdim * (bStruct.nH-1.0), cubic = -0.5 )
gc_bPhi  = interpolate ( bStruct.bPhi, ( rr_array[iiFirstGyration,0] - bStruct.rleft ) / bStruct.rdim * (bStruct.nW-1.0), $
    ( rr_array[iiFirstGyration,2] - min ( bStruct.z ) ) / bStruct.zdim * (bStruct.nH-1.0), cubic = -0.5 )
gc_bz  = interpolate ( bStruct.bz, ( rr_array[iiFirstGyration,0] - bStruct.rleft ) / bStruct.rdim * (bStruct.nW-1.0), $
    ( rr_array[iiFirstGyration,2] - min ( bStruct.z ) ) / bStruct.zdim * (bStruct.nH-1.0), cubic = -0.5 )

gc_bMag	= sqrt ( gc_bR^2 + gc_bPhi^2 + gc_bz^2 )

gc_bR_B	= gc_bR / gc_bMag
gc_bPhi_B	= gc_bPhi / gc_bMag
gc_bz_B	= gc_bz / gc_bMag

vMag	= sqrt ( vv_array[iiFirstGyration,0]^2 + vv_array[iiFirstGyration,1]^2 + vv_array[iiFirstGyration,2]^2 )
vPar_all	= ( vv_array[iiFirstGyration,0] * gc_bR_B + vv_array[iiFirstGyration,1]* gc_bPhi_B + vv_array[iiFirstGyration,2] * gc_bz_B )
vPar	= mean ( vPar_all )
vPerp_all	= sqrt ( vMag^2 - vPar_all^2 )
vPerp	= mean ( vPerp_all )

u   = mi * vPerp^2 / ( 2.0 * mean ( gc_bMag ) )

plots, rr_array[iiFirstGyration,0], rr_array[iiFirstGyration,2], color = 200, thick = 3.0
plots, [r,r], [z-0.6,z+0.6], lineStyle = 1
plots, [r-0.3,r+0.3], [z,z], lineStyle = 1

; Trace field line for comparison

fl_rArray	= r
fl_zArray	= z

fl_pos	= [ r, 0.0, z ]

dPhi    = 2 * !pi / 100.0
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
	fl_zArray  = [ fl_zArray, fl_pos[2] ]

endFor

;	END TEST
;-----------------------------------
;-----------------------------------

pos	= [ R[0], 0.0, z[0] ]

rTrack  = R[0]
zTrack  = z[0]
vPerpTrack  = vPerp[0]
vParTrack   = vPar[0]
dtArray	= dt

distance = 0

stepCnt	= 0
firstOrbit	= 1
tau	= 0.0
dTau	= 0.0

while stillIn AND firstOrbit do begin

    vPerp   = dlg_vPerp ( pos, u, bStruct ) 
    vgc = dlg_gc_velocity ( vPerp, vPar, pos, bStruct )
    k1_vPar   = dt * dlg_vPar ( pos, u, bStruct ) 
    k1_vgc  = dt * vgc

    vPerp   = dlg_vPerp ( pos + k1_vgc / 2.0, u, bStruct ) 
    vgc = dlg_gc_velocity ( vPerp, vPar + k1_vPar / 2.0, pos + k1_vgc / 2.0, bStruct )
    k2_vPar   = dt * dlg_vPar ( pos + k1_vgc / 2.0, u, bStruct ) 
    k2_vgc  = dt * vgc
 
    vPerp   = dlg_vPerp ( pos + k2_vgc / 2.0, u, bStruct ) 
    vgc = dlg_gc_velocity ( vPerp, vPar + k2_vPar / 2.0, pos + k2_vgc / 2.0, bStruct )
    k3_vPar   = dt * dlg_vPar ( pos + k2_vgc / 2.0, u, bStruct ) 
    k3_vgc  = dt * vgc

    vPerp   = dlg_vPerp ( pos + k3_vgc, u, bStruct ) 
    vgc = dlg_gc_velocity ( vPerp, vPar + k3_vPar, pos + k3_vgc, bStruct )
    k4_vPar   = dt * dlg_vPar ( pos + k3_vgc, u, bStruct ) 
    k4_vgc  = dt * vgc
    
    vPar    = vPar + ( k1_vPar + 2.0 * k2_vPar + 2.0 * k3_vPar + k4_vPar ) / 6.0
    pos   = pos + ( k1_vgc + 2.0 * k2_vgc + 2.0 * k3_vgc + k4_vgc ) / 6.0

	stillIn	= dlg_checkIn ( pos, bStruct )
    
	rTrack  = [ rTrack, pos[0] ]
    zTrack  = [ zTrack, pos[2] ]
    vPerpTrack  = [ vPerpTrack, vPerp[0] ]
    vParTrack   = [ vParTrack, vPar[0] ]
	dtArray	= [ dtArray, dt ]

	++stepCnt
   	
endwhile



end
