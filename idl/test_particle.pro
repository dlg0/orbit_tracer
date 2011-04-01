pro test_particle

    proton_mass	= 1.67e-27
	e	= 1.602e-19
	amu	= 1
	z	= 1
	m	= amu * proton_mass 
	q	= z * e
	
	;eqdsk_fName	= '../eqdsk/eqdsk.122993'
	eqdsk_fName	= '/Users/dg6/scratch/sMC+rf/test_qlSerial/data/eqdsk'
	eqdsk = readGEqdsk ( eqdsk_fName )	

	gc_terms, m, eqdsk, gc_struct = bStruct

	r   = 0.696515 
	p	= 0.0
	z   =-0.0242857 

	vr = 21933.6
	vp = -226865
	vz = 17501.1

    rr_cyl = [ r, p, z]
    vv_cyl = [ vr, vp, vz]

    nSteps = 3

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
			/guiding


end
