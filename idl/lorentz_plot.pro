function dlg_interpB, r, bStruct, bMag = bMag, unit = unit

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

	if keyword_set ( unit ) then begin

		return, [ bR, bPhi, bz ] / bMag

	endif else begin

		return, [ bR, bPhi, bz ]

	endelse

end


function dlg_vxB, r_car, v_car, bStruct

	xyMag	= sqrt ( r_car[0]^2 + r_car[1]^2 )

	cyl2car	= [	[ r_car[0] / xyMag, - r_car[1] / xyMag, 0 ], $
				[ r_car[1] / xyMag, r_car[0] / xyMag, 0 ], $
				[ 0,0,1] ] 

	car2cyl	= invert ( cyl2car )

	r_cyl	= [ [ sqrt ( r_car[0]^2+r_car[1]^2 ) ], $
				[ aTan ( r_car[1], r_car[0] ) ], $
				[ r_car[2] ] ]

	b_cyl   = dlg_interpB ( r_cyl, bStruct, bMag = bMag )

	b_car	= cyl2car ## b_cyl

	vxB_R	= v_car[1] * b_car[2] - v_car[2] * b_car[1]
	vxB_phi	= -1.0 * ( v_car[0] * b_car[2] - v_car[2] * b_car[0] )
	vxB_z	= v_car[0] * b_car[1] - v_car[1] * b_car[0]

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

	b	= dlg_interpB ( r, bStruct, bMag = bMag )

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

	unitb	= dlg_interpB ( r, bStruct, /unit ) 

    vgc_R   = vPar * unitb[0] + vPer^2 * grad_R + vPar^2 * curv_R 
    vgc_phi   = vPar * unitb[1] + vPer^2 * grad_phi + vPar^2 * curv_phi
    vgc_z   = vPar * unitb[2] + vPer^2 * grad_z + vPar^2 * curv_z

    return, [ vgc_R, vgc_phi, vgc_z ]

end 

pro lorentz_plot, $
	noPlot = noPlot, $
	vv_cyl = vv_cyl, $
	rr_cyl = rr_cyl, $
	nSteps = nSteps, $
	x_gc = x_gc, y_gc = y_gc, z_gc = z_gc, $
	x_lorentz = x_lorentz, y_lorentz = y_lorentz, z_lorentz = z_lorentz, $
	vx_lorentz = vx_lorentz, vy_lorentz = vy_lorentz, vz_lorentz = vz_lorentz, $
	dt = dt, $
	delay = delay, $
	distance = distance, $
	match = match, $
	fieldLine = fieldLine, $
	neutralTrack = neutralTrack, $
	vPerTrack = vPerTrack, vParTrack = vParTrack, $
	stillIn = stillIn, $
	bStruct = bStruct, $
	eqdsk = eqdsk, $
	mass = mi, $
	q = q, $
	lorentz = lorentz, $
	guiding = guiding, $
	vPerIn = vPerIn, vParIn = vParIn, $
	pitchScatter = pitchScatter

	if keyword_set ( distance ) then begin
		distance_ = distance
	endif else begin
		distance = 0
	endelse
	
	if not keyword_set ( vv_cyl ) then vv_cyl = [[0.0],[0.0],[0.0]]
	;if not keyword_set ( vv_cyl ) then begin	
	;	r   = 2.1
	;	z   = 0.0 
	;	phi	= 0.0
	;	
	;	vR   	= 0e6
	;	vz    	= 0e6
	;	vPhi	= 5e6 

	;	rr_cyl	= transpose ( [ R, phi, z ] )
	;	vv_cyl	= transpose ( [ vR, vPhi, vz ] )
	;endif

	cyl2car	=	[	[ cos ( rr_cyl[1] ), sin ( rr_cyl[1] ), 0 ], $
					[ -sin ( rr_cyl[1] ), cos ( rr_cyl[1] ), 0 ], $
					[ 0, 0, 1 ] ]

	car2cyl	= invert ( cyl2car )

	vv_car	= cyl2car ## vv_cyl 

	rr_car	= [ [ rr_cyl[0] * cos ( rr_cyl[1] ) ], $
				[ rr_cyl[0] * sin ( rr_cyl[1] ) ], $
				[ rr_cyl[2] ] ]

	rr_car_start	= rr_car
	rr_cyl_start	= rr_cyl
	vv_car_start	= vv_car
	vv_cyl_start	= vv_cyl
 
	vMag	= sqrt ( vv_cyl[0]^2 + vv_cyl[1]^2 + vv_cyl[2]^2 )
	en_	= mi * vMag^2 / 2.0 / 1.602e-19 * 1e-3; [keV]
	;print, en_, 'keV';, rr_cyl[*], vv_cyl[*]
	
	if not keyword_set ( dt ) then dt = 0.01e-7 
	if not keyword_set ( nSteps ) then nSteps = 2000


	if keyword_set ( lorentz ) then begin

		;	Here is the Lorentz eqn integration, done in 
		;	cartesian coords so we have to rotate the b field
		;	and velocity vectors to cartesian first, and also
		;	convert the position 

		vv_car_array	= vv_car
		rr_car_array	= rr_car
		tt_array	= 0

		;vv_cyl_array	= vv_cyl
		;rr_cyl_array	= rr_cyl

		for i = 0, nSteps - 2 do begin

			;if dt * i lt delay then begin
			if distance gt 0 then begin

				rr_car	= rr_car + vv_car * dt 
				distance = distance - vMag * dt

			endif else begin

				vxB	= dlg_vxB ( rr_car, vv_car, bStruct )
				k1_v	= dt * q / mi * vxB	
				k1_r	= dt * ( vv_car )
		
				vxB	= dlg_vxB ( rr_car + k1_r / 2.0, vv_car + k1_v / 2.0, bStruct )
				k2_v	= dt * q / mi * vxB	
				k2_r	= dt * ( vv_car + k1_v / 2.0 )
		
				vxB	= dlg_vxB ( rr_car + k2_r / 2.0, vv_car + k2_v / 2.0, bStruct )
				k3_v	= dt * q / mi * vxB	
				k3_r	= dt * ( vv_car + k2_v / 2.0 )
		
				vxB	= dlg_vxB ( rr_car + k3_r, vv_car + k3_v, bStruct )
				k4_v	= dt * q / mi * vxB	
				k4_r	= dt * ( vv_car + k3_v )
		
				vv_car	= vv_car + ( k1_v + 2.0 * k2_v + 2.0 * k3_v + k4_v ) / 6.0
				rr_car	= rr_car + ( k1_r + 2.0 * k2_r + 2.0 * k3_r + k4_r ) / 6.0

			endelse
		
			vv_car_array	= [ vv_car_array, vv_car ]
			rr_car_array	= [ rr_car_array, rr_car ]
			tt_array	= [ tt_array, i * dt ]

			;xyMag	= sqrt ( rr_car[0]^2 + rr_car[1]^2 )

			;cyl2car	= [	[ rr_car[0] / xyMag, - rr_car[1] / xyMag, 0 ], $
			;			[ rr_car[1] / xyMag, rr_car[0] / xyMag, 0 ], $
			;			[ 0,0,1] ] 

			;car2cyl	= invert ( cyl2car )

			;vv_cyl	= car2cyl ## vv_car
			;rr_cyl	= [	[ xyMag ], $
			;			[ aTan ( rr_car[1], rr_car[0] ) ], $
			;			[ rr_car[2] ] ]
		
			;vv_cyl_array	= [ vv_cyl_array, vv_cyl ]
			;rr_cyl_array	= [ rr_cyl_array, rr_cyl ]

			;print, sqrt ( vv_car[0]^2 + vv_car[1]^2 + vv_car[2]^2 )
			;print, sqrt ( vv_cyl[0]^2 + vv_cyl[1]^2 + vv_cyl[2]^2 )

		endfor

	endif

	if keyword_set ( match ) then begin

		;	Find mean guiding center start point from 
		;	approx. one gyration

		;	Use magnetic field here for a guess at the 
		;	gyro frequency, that's all
		
		pos	= rr_cyl
		bHere_g   = dlg_interpB ( pos, bStruct, bMag = g_bMag )

		omega0	= abs ( 2.0 * q * g_bMag / mi)
		iiFirstGyration	= where ( tt_array lt 4.0 * !pi / omega0 )

		R_gc	= mean ( rr_cyl_array[iiFirstGyration,0] )
		phi_gc	= mean ( rr_cyl_array[iiFirstGyration,1] )
		z_gc	= mean ( rr_cyl_array[iiFirstGyration,2] )

		;phi_gc = phi

		pos	= [ R_gc, phi_gc, z_gc ]

		bHere_gc   = dlg_interpB ( pos, bStruct, bMag = bMag_gc )

		vPar	= vv_cyl # bHere_gc / bMag_gc
		vPer	= sqrt ( vMag^2 - vPar^2 )

	endif else begin

		pos		= transpose ( rr_cyl_start )
		bHere_gc   = dlg_interpB ( pos, bStruct, bMag = bMag_gc )
		vPar	= vv_cyl_start # bHere_gc / bMag_gc
		vPer	= sqrt ( vMag^2 - vPar^2 )
	
	endelse

	if keyword_set ( vPerIn ) then vPer = vPerIn
	if keyword_set ( vParIn ) then vPar = vParIn


	; Trace field line for comparison
	if keyword_set ( fieldLine ) then begin
	
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

	endif


	if keyword_set ( guiding ) then begin

		;	GC integration

		bHere   = dlg_interpB ( pos, bStruct, bMag = bMag )
		u   = mi * vPer^2 / ( 2.0 * bMag )

		rTrack	= fltArr ( nSteps )
		phiTrack	= fltArr ( nSteps )
		zTrack	= fltArr ( nSteps )
		vPerTrack	= fltArr ( nSteps )
		vParTrack	= fltArr ( nSteps )

		neutralTrack	= intArr ( nSteps )
		tTrack	= fltArr ( nSteps )

		rTrack[0]  	= pos[0]
		phiTrack[0] 	= pos[1]
		zTrack[0]  	= pos[2] 
		vPerTrack[0]  	= vPer
		vParTrack[0]   = vPar
		neutralTrack[0]	= 1
		tTrack[0]	= 0

		rr_car	= rr_car_start
		vv_car	= vv_car_start

		if keyword_set ( distance ) then distance = distance_

		stillIn = 1
		for i = 0, nSteps - 2 do begin
			
			;if dt * i lt delay then neutral = 1 else neutral = 0
			if distance gt 0 then neutral = 1 else neutral = 0
			if neutral eq 1 then begin

				rr_car	= rr_car + vv_car * dt 
				rr_cyl	= [ [ sqrt ( rr_car[0]^2+rr_car[1]^2 ) ], $
					[ aTan ( rr_car[1], rr_car[0] ) ], $
					[ rr_car[2] ] ]

				pos	= transpose ( rr_cyl )

				bHere_gc   = dlg_interpB ( pos, bStruct, bMag = bMag_gc )

				vPar	= -vv_cyl_start # bHere_gc / bMag_gc
				vPer	= sqrt ( vMag^2 - vPar^2 )

				bHere   = dlg_interpB ( pos, bStruct, bMag = bMag )
				u   = mi * vPer^2 / ( 2.0 * bMag )

				distance = distance - vMag * dt

				if distance le 0 then begin

					if keyword_set ( pitchScatter ) then begin

						; only changing vper here since i want an artifically
						; large scatter but also keep the particles headed in the
						; same direction relative to the field.

						vParBak	= vPar
						vMag	= sqrt ( vPer^2 + vPar^2 )
						pitch	= !pi-aTan ( vPer, vPar )

						pitch	= pitch + pitchScatter

						vPar	= -cos ( pitch ) * vMag
						vPer	= sqrt ( vMag^2 - vPar^2 )
						vPar = vParBak
	
					endif

				endif

			endif else begin

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

			endelse
		
			stillIn	= dlg_checkIn ( pos, bStruct )
			;if stillIn eq 0 and dt * i gt delay then div = 0 else div = 1
			;if stillIn eq 0 and distance le 0 then div = 0 else div = 1
			;
			rTrack[i+1]  = pos[0] ;/ div 
			phiTrack[i+1]	= pos[1] ;/ div 
		    zTrack[i+1]  = pos[2] ;/ div 

		    vPerTrack[i+1]  = vPer ;/ div 
		    vParTrack[i+1]   = vPar ;/ div 

			neutralTrack[i+1]	= neutral ;/ div 
			tTrack[i+1]	= dt * i ;/ div 

			if stillIn eq 0 and distance le 0 then return
		
		endfor

		if i lt nSteps-2 then begin
			
			NaN	= sqrt ( -1 )
			rTrack[i:*]	= NaN
			phiTrack[i:*]	= NaN
			zTrack[i:*]	= NaN
			vPerTrack[i:*]	= NaN
			vParTrack[i:*]	= NaN
			neutralTrack[i:*]	= NaN
			tTrack[i:*]	= NaN

		endif

	endif

	;	convert cylindrical coords to cartesian for plotting

	if keyword_set ( lorentz ) then begin

		x_lorentz	= rr_car_array[*,0]
		y_lorentz	= rr_car_array[*,1]
		z_lorentz	= rr_car_array[*,2]

		vx_lorentz	= vv_car_array[*,0]
		vy_lorentz	= vv_car_array[*,1]
		vz_lorentz	= vv_car_array[*,2]

	endif
	
	if keyword_set ( guiding ) then begin
	
		x_gc	= rTrack * cos ( phiTrack )
		y_gc	= rTrack * sin ( phiTrack )
		z_gc	= zTrack

	endif 

	if keyword_set ( fieldLine ) then begin
	
		x_fl	= fl_rArray * cos ( fl_phiArray )
		y_fl	= fl_rArray * sin ( fl_phiArray )
		z_fl	= fl_zArray

	endif

	if not keyword_set ( noPlot ) then begin

   		loadct, 12, /sil, rgb_table = ct12
		red	= transpose ( ct12[12*16-1,*] )
		blue	= transpose ( ct12[8*16-1,*] )
		green	= transpose ( ct12[1*16-1,*] )

		iPlot, x_lorentz, y_lorentz, z_lorentz, $
			rgb_table = 1, $
			vert = fIndGen(nSteps)/nSteps*255, $
			/zoom_on_resize, $
			thick = 3, /iso
		iPlot, x_gc, y_gc, z_gc, /over, $
			rgb_table = 3, $
			vert = fIndGen(nSteps)/nSteps*255, $
			thick = 3
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
		iPlot, rr_cyl_array[*,0], rr_cyl_array[*,2], /over, $
			color = green, $
			thick = 2, $
			trans = 50
		iPlot, eqdsk.rbbbs, eqdsk.zbbbs, /over
		iPlot, eqdsk.rLim, eqdsk.zLim, /over

	endif

	close, /all

end
