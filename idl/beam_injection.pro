pro beam_injection, $
	plot = plot

	;eqdsk_fName	= '../eqdsk/eqdsk.122993'
	eqdsk_fName	= '../eqdsk/Scen4_bn2.57_129x129'
	eqdsk = readGEqdsk ( eqdsk_fName )	


	beamStart	= [ -399.33, 30244.7, 1442.0 ] * 1e-3
	beam2		= [ 5399.88, 0.0, -130.5 ] * 1e-3

	beamMag		= sqrt ( (beam2[0]-beamStart[0])^2 $
					+ (beam2[1]-beamStart[1])^2 $
					+ (beam2[2]-beamStart[2]) )
	beamUnitVec	= ( beam2 - beamStart ) / beamMag

	r_car	= beamStart
	
	v_car	= beamUnitVec * 3e7

	proton_mass	= 1.67e-27
	e	= 1.602e-19
	amu	= 2
	z	= 2
	m	= amu * proton_mass 
	q	= z * e
	vMag	= sqrt ( v_car[0]^2 + v_car[1]^2 + v_car[2]^2 )
	energy	= 0.5 * m * vMag^2 / e	; [ev]

	gc_terms, m, eqdsk, gc_struct = bStruct
	
	nP	= 2000L
	nSteps	= 10000L

	imsl_randomOpt, set = 12345

	dt	= 0.05e-7
	offset	= 0.30e-6; [s] 70 * dt  
	spread	= 0.015e-6 
	delay	= imsl_random ( nP, /normal ) * spread + offset
	distance	= imsl_random ( nP, /normal ) * 2.0 + (23.4+6.58) - 1.0

	deflection	= 0.0025
	rotTh1	=	imsl_random ( nP, /normal ) * deflection
	rotTh2	=	imsl_random ( nP, /normal ) * deflection
	rotTh3	=	imsl_random ( nP, /normal ) * deflection

	scatter	= 0.5
	pitchScatter	=	imsl_random ( nP, /normal ) * scatter * !pi/2

	new_v	= fltArr ( nP, 3 )
	v_cyl	= fltArr ( nP, 3 )

	r_cyl	= [ [ sqrt ( r_car[0]^2+r_car[1]^2 ) ], $
				[ aTan ( r_car[1], r_car[0] ) ], $
				[ r_car[2] ] ]

	cyl2car	=	[	[ cos ( r_cyl[1] ), sin ( r_cyl[1] ), 0 ], $
					[ -sin ( r_cyl[1] ), cos ( r_cyl[1] ), 0 ], $
					[ 0, 0, 1 ] ]

	car2cyl	= invert ( cyl2car )

	for i=0,nP-1 do begin

		;	y towards z
		rot_x	= [ [ 1, 0, 0], $
					[ 0, cos ( rotTh1[i] ), sin ( rotTh1[i] ) ], $
					[ 0, -sin ( rotTh1[i] ), cos ( rotTh1[i] ) ] ]

		;	z towards x
		rot_y	= [ [ cos ( rotTh2[i] ), 0, -sin ( rotTh2[i] ) ], $
					[ 0, 1, 0 ], $
					[ sin ( rotTh2[i] ), 0, cos ( rotTh2[i] ) ] ]

    	rot_z  = [ [ cos ( rotTh3[i] ), -sin ( rotTh3[i] ), 0 ], $
    	           [ sin ( rotTh3[i] ), cos ( rotTh3[i] ), 0 ], $
    	           [ 0, 0, 1 ] ] 


		new_v[i,*]	= rot_x ## v_car
		new_v[i,*]	= rot_y ## new_v[i,*]
		new_v[i,*]	= rot_z ## new_v[i,*]

		v_cyl[i,*]	= car2cyl ## new_v[i,*]

	endfor


	;	Calculate trajectories

	x_gc	= fltArr ( nSteps, nP )
	y_gc	= fltArr ( nSteps, nP )
	z_gc	= fltArr ( nSteps, nP )
	
	x_lorentz	= fltArr ( nSteps, nP )
	y_lorentz	= fltArr ( nSteps, nP )
	z_lorentz	= fltArr ( nSteps, nP )
	
	vx_lorentz	= fltArr ( nSteps, nP )
	vy_lorentz	= fltArr ( nSteps, nP )
	vz_lorentz	= fltArr ( nSteps, nP )
	
	time	= fIndGen ( nSteps ) * dt

	vPer	= fltArr ( nSteps, nP )
	vPar	= fltArr ( nSteps, nP )

	neutral	= intArr ( nSteps, nP )

	wall	= intArr ( nP )

	r_cyl_start	= r_cyl
	for i=0,nP-1 do begin

		lorentz_plot, /noPlot, $
			vv_cyl = v_cyl[i,*], $
			rr_cyl = r_cyl, $
			nSteps = nSteps, $
			x_gc = x_gc_tmp, $
			y_gc = y_gc_tmp, $
			z_gc = z_gc_tmp, $
			distance = distance[i], $
			dt = dt, $
			stillIn = stillIn, $
			vPerTrack = vPerTrack, $
			vParTrack = vParTrack, $
			neutralTrack = neutralTrack, $
			eqdsk = eqdsk, $
			bStruct = bStruct, $
			mass = m, q = q, $
			/guiding, $
			pitchScatter = pitchScatter[i]

		x_gc[*,i]	= x_gc_tmp
		y_gc[*,i]	= y_gc_tmp
		z_gc[*,i]	= z_gc_tmp

		vPer[*,i]	= vPerTrack
		vPar[*,i]	= vParTrack
		
		wall[i]	= stillIn

		neutral[*,i]	= neutralTrack

		r_cyl	= r_cyl_start

		print, i, stillIn

	endfor	


	;	Write netCDF file containing track data

	id = nCdf_create ( 'beamParticleData.nc', /clobber )
	nCdf_control, id, /fill

	tDim_id	= nCdf_dimDef ( id, 'time', nSteps )
	nDim_id	= nCdf_dimDef ( id, 'nP', nP )

	t_id	= nCdf_varDef ( id, 'time', [tDim_id], /float )

	xgc_id	= nCdf_varDef ( id, 'xgc', [tDim_id,nDim_id], /float )
	ygc_id	= nCdf_varDef ( id, 'ygc', [tDim_id,nDim_id], /float )
	zgc_id	= nCdf_varDef ( id, 'zgc', [tDim_id,nDim_id], /float )

	vPer_id	= nCdf_varDef ( id, 'vPer', [tDim_id,nDim_id], /float )
	vPar_id	= nCdf_varDef ( id, 'vPar', [tDim_id,nDim_id], /float )

	wall_id	= nCdf_varDef ( id, 'wall', [nDim_id], /short )
	neutral_id	= nCdf_varDef ( id, 'neutral', [tDim_id,nDim_id], /short )

	nCdf_attPut, id, t_id, 'units', 's'

	nCdf_attPut, id, xgc_id, 'units', 'm'
	nCdf_attPut, id, ygc_id, 'units', 'm'
	nCdf_attPut, id, zgc_id, 'units', 'm'

	nCdf_attPut, id, xgc_id, 'long_name', 'x position'
	nCdf_attPut, id, ygc_id, 'long_name', 'y position'
	nCdf_attPut, id, zgc_id, 'long_name', 'z position'

	nCdf_attPut, id, vPer_id, 'units', 'm/s'
	nCdf_attPut, id, vPar_id, 'units', 'm/s'
	
	nCdf_attPut, id, vPer_id, 'long_name', 'Perpendicular velocity'
	nCdf_attPut, id, vPar_id, 'long_name', 'Parallel velocity'

	nCdf_attPut, id, neutral_id, 'long_name', 'Neutral status, i.e., does the particle feel the B field[0] or not[1]'
	nCdf_attPut, id, wall_id, 'long_name', 'Wall status, i.e., does the particle hit the wall[1] or not[0]'

	nCdf_attPut, id, /global, 'Title', 'ITER test particle trajectory data'
	nCdf_attPut, id, /global, 'green_particle_data', 'green_particle_data'

	nCdf_control, id, /enDef

	nCdf_varPut, id, xgc_id, x_gc
	nCdf_varPut, id, ygc_id, y_gc
	nCdf_varPut, id, zgc_id, z_gc

	nCdf_varPut, id, vPer_id, vPer 
	nCdf_varPut, id, vPar_id, vPar

	nCdf_varPut, id, t_id, time

	neutral	= fltArr ( nSteps, nP )
	nCdf_varPut, id, neutral_id, neutral
	nCdf_varPut, id, wall_id, abs(wall-1)

	nCdf_close, id


	if keyword_set ( plot ) then begin
	
		loadct, 12, /sil, rgb_table = ct12
		red	= transpose ( ct12[12*16-1,*] )
		blue	= transpose ( ct12[8*16-1,*] )
		green	= transpose ( ct12[1*16-1,*] )

		x_plot	= x_gc
		y_plot	= y_gc
		z_plot	= z_gc

		iPlot, x_plot[*,0], y_plot[*,0], z_plot[*,0], $
			rgb_table = 3, $
			vert = 255-fIndGen(nSteps)/nSteps*255, $
			thick = 2, $
			/zoom_on_resize, $
			/iso, $
			trans = 50
		
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

		;for i = 0, 20 do begin

		;	lim_phi	= (eqdsk.rlim * 0 + i * !pi * 4 / 180 ) + 40 * !dtor
		;	x_lim	= eqdsk.rlim * cos ( lim_phi )
		;	y_lim	= eqdsk.rlim * sin ( lim_phi )
		;	z_lim	= eqdsk.zlim

		;	iPlot, x_lim, y_lim, z_lim, $
		;		/over, $
		;		trans = 50, $
		;		thick = 20 
		;endfor


		for i=1,nP-1 do begin

			iPlot, x_plot[*,i], y_plot[*,i], z_plot[*,i], $
				rgb_table = 3, $
				vert = 255-fIndGen(nSteps)/nSteps*255, $
				thick = 2, $
				/over, $
				trans = 50
		
		endfor	

	endif
stop

end
