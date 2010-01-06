pro beam_injection

	;eqdsk_fName	= '../eqdsk/eqdsk.122993'
	eqdsk_fName	= '../eqdsk/Scen4_bn2.57_129x129'
	eqdsk = readGEqdsk ( eqdsk_fName )	

	r_car	= [ [ 10.0 ], $
				[ 0.0 ], $
				[ -0.4 ] ]
	
	v_car	= [ [ -0.7e6 ], $
				[ 0.5e6 ], $
				[ 0.2e6 ] ] * 8 

	proton_mass	= 1.67e-27
	e	= 1.602e-19
	amu	= 2
	z	= 2
	m	= amu * proton_mass 
	q	= z * e
	vMag	= sqrt ( v_car[0]^2 + v_car[1]^2 + v_car[2]^2 )
	energy	= 0.5 * m * vMag^2 / e	; [ev]

	gc_terms, m, eqdsk, gc_struct = bStruct
	
	nP	= 100 
	nSteps	= 1000

	imsl_randomOpt, set = 12345

	dt	= 0.1e-7
	offset	= 10e-7; [s] 70 * dt  
	spread	= 10 
	delay	= imsl_random ( nP, /normal ) * 0.2e-7 * spread + offset

	deflection	= 0.05
	rotTh1	=	imsl_random ( nP, /normal ) * deflection
	rotTh2	=	imsl_random ( nP, /normal ) * deflection
	rotTh3	=	imsl_random ( nP, /normal ) * deflection

	new_E	= energy - imsl_random ( nP, /exp ) * fIndGen ( nP ) * 400

	iiNeg	= where ( new_E lt 0, iiNegCnt )
	if iiNegCnt gt 0 then new_E[iiNeg] = -new_E[iiNeg]

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

		new_v[i,*]	= new_v[i,*] * new_E[i] / energy
		v_cyl[i,*]	= car2cyl ## new_v[i,*]

	endfor


	;	Calculate trajectories

	x_gc	= fltArr ( nSteps, nP )
	y_gc	= fltArr ( nSteps, nP )
	z_gc	= fltArr ( nSteps, nP )
	
	time	= fIndGen ( nSteps ) * dt

	vPer	= fltArr ( nSteps, nP )
	vPar	= fltArr ( nSteps, nP )

	neutral	= intArr ( nSteps, nP )

	wall	= intArr ( nP )

	r_cyl_start	= r_cyl
	for i=0,nP-1 do begin

		print, i
		lorentz_plot, /noPlot, $
			vv_cyl = v_cyl[i,*], $
			rr_cyl = r_cyl, $
			nSteps = nSteps, $
			x_gc = x_gc_tmp, $
			y_gc = y_gc_tmp, $
			z_gc = z_gc_tmp, $
			delay = delay[i], $
			dt = dt, $
			stillIn = stillIn, $
			vPerTrack = vPerTrack, $
			vParTrack = vParTrack, $
			neutralTrack = neutralTrack, $
			eqdsk = eqdsk, $
			bStruct = bStruct, $
			mass = m, q = q

		x_gc[*,i]	= x_gc_tmp
		y_gc[*,i]	= y_gc_tmp
		z_gc[*,i]	= z_gc_tmp

		vPer[*,i]	= vPerTrack
		vPar[*,i]	= vParTrack
		
		wall[i]	= stillIn

		neutral[*,i]	= neutralTrack

		r_cyl	= r_cyl_start
	endfor	

	loadct, 12, /sil, rgb_table = ct12
	red	= transpose ( ct12[12*16-1,*] )
	blue	= transpose ( ct12[8*16-1,*] )
	green	= transpose ( ct12[1*16-1,*] )
	
	iPlot, x_gc[*,0], y_gc[*,0], z_gc[*,0], $
		rgb_table = 3, $
		vert = 255-fIndGen(nSteps)/nSteps*255, $
		thick = new_E[0] / 1e3 / 250, $
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

		iPlot, x_gc[*,i], y_gc[*,i], z_gc[*,i], $
			rgb_table = 3, $
			vert = 255-fIndGen(nSteps)/nSteps*255, $
			thick = new_E[i] / 1e3 / 250, $
			/over, $
			trans = 50
	
	endfor	

	;	Write netCDF file containing track data

	id = nCdf_create ( 'particleData.nc', /clobber )
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

	nCdf_control, id, /enDef

	nCdf_varPut, id, xgc_id, x_gc
	nCdf_varPut, id, ygc_id, y_gc
	nCdf_varPut, id, zgc_id, z_gc

	nCdf_varPut, id, vPer_id, vPer 
	nCdf_varPut, id, vPar_id, vPar

	nCdf_varPut, id, t_id, time

	nCdf_varPut, id, neutral_id, neutral
	nCdf_varPut, id, wall_id, abs(wall-1)

	nCdf_close, id
	
stop

end
