pro bulk_plasma, plot = plot

	time0	= sysTime ( /seconds )

	proton_mass	= 1.67e-27
	e	= 1.602e-19
	amu	= 2
	z	= 2
	m	= amu * proton_mass 
	q	= z * e
	
	;eqdsk_fName	= '../eqdsk/eqdsk.122993'
	eqdsk_fName	= '~/code/orbit_tracer/eqdsk/Scen4_bn2.57_129x129'
	eqdsk = readGEqdsk ( eqdsk_fName )	

	gc_terms, m, eqdsk, gc_struct = bStruct

	;	setup launch location grid in phase space

	nX	= 10L
	minX	= min ( eqdsk.rbbbs )
	maxX	= max ( eqdsk.rbbbs )
	xStep	= ( maxX - minX ) / ( nX - 1.0 )
	xCoord	= fIndGen ( nX ) * xStep + minX

	;nX = 1
	;xCoord = 7.0

	nZ	= 21L
	minZ	= min ( eqdsk.zbbbs )
	maxZ	= max ( eqdsk.zbbbs )
	zStep	= ( maxZ - minZ ) / ( nZ - 1.0 )
	zCoord	= fIndGen ( nZ ) * zStep + minZ

	;nZ = 1
	;zCoord = 1

	vScale	= 3e7
	nvPer	= 4L
	vPer	= fIndGen ( nvPer ) / (nvPer - 1) * vScale * 2 + vScale / nvPer

	nvPar	= 2L
	vPar	= [-1.0,1.0] * vScale  

	nSteps	= 10000L 

	dt	= 0.01e-7

	nP	= nX * nZ * nvPer * nvPar

	x_gc	= fltArr ( nSteps, nP )
	y_gc	= fltArr ( nSteps, nP )
	z_gc	= fltArr ( nSteps, nP )

	vPer_gc	= fltArr ( nSteps, nP )
	vPar_gc	= fltArr ( nSteps, nP )

	wall	= intArr ( nP )

	x_gc_tmp	= fltArr ( nSteps )
	y_gc_tmp	= fltArr ( nSteps )
	z_gc_tmp	= fltArr ( nSteps )
	
	vPerTrack	= fltArr ( nSteps )
	vParTrack	= fltArr ( nSteps )
	imsl_randomOpt, set = 12345

	phi	= imsl_random ( nP, /uniform ) * 2 * !pi

	pno	= 0
	for i=0,nX-1 do begin
		for j=0,nZ-1 do begin
			for k=0,nvPer-1 do begin
				for l=0,nvPar-1 do begin

					vMag	= sqrt ( vPar[l]^2 + vPer[k]^2  )
					energy	= 0.5 * m * vMag^2 / e	; [ev]

					rr_cyl	= [ [ xCoord[i] ], $
								[ phi[pno] ], $
								[ zCoord[j] ] ]

					lorentz_plot, /noPlot, $
						rr_cyl = rr_cyl, $
						nSteps = nSteps, $
						x_gc = x_gc_tmp, $
						y_gc = y_gc_tmp, $
						z_gc = z_gc_tmp, $
						dt = dt, $
						stillIn = stillIn, $
						vPerTrack = vPerTrack, $
						vParTrack = vParTrack, $
						eqdsk = eqdsk, $
						bStruct = bStruct, $
						mass = m, q = q, $
						vPerIn = vPer[k], vParIn = vPar[l], $
						/guiding

					x_gc[*,pno]	= x_gc_tmp
					y_gc[*,pno]	= y_gc_tmp
					z_gc[*,pno]	= z_gc_tmp

					vPer_gc[*,pno]	= vPerTrack
					vPar_gc[*,pno]	= vParTrack
					
					wall[pno]	= stillIn

					pno++

					print, pno, i, j, k, l, energy * 1e-3, stillIn

				endfor
			endfor
		endfor
	endfor	

	time1	= sysTime ( /seconds )

	print, 'Time Taken: ', time1-time0

	time	= fIndGen ( nSteps ) * dt

	;	Write netCDF file containing track data

	id = nCdf_create ( 'bulkParticleData.nc', /clobber )
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

;	plot

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


	for i=0, nP -1 do begin

		iiPlot	= where ( x_plot[*,i] eq x_plot[*,i], iiPlotCnt )
		if iiPlotCnt gt 0 then $
		iPlot, x_plot[iiPlot,i], y_plot[iiPlot,i], z_plot[iiPlot,i], $
			rgb_table = 3, $
			vert = 255-fIndGen(nSteps)/nSteps*255, $
			thick = 2, $
			/over, $
			trans = 50
	
	endfor	

	endif
	
end
