pro gc_terms, amu, Z, $
    eqdsk = eqdsk, $
	gc_struct = gc_struct, $
	curv = curv, $
	grad = grad, $
    ar2InputFileName = ar2InputFileName, $
    oned = oned

    @constants
    q = Z * _e
    m = amu * mi

    if keyword_set(eqdsk) then begin

	    bT = eqdsk.bPhi 
	    bR = eqdsk.bR
	    bz = eqdsk.bz
        bMag = eqdsk.bMag
        dR = eqdsk.rStep
        dZ = eqdsk.zStep
        nR = eqdsk.nW
        nZ = eqdsk.nH
        r  = eqdsk.r
        z  = eqdsk.z

    endif else if keyword_set(ar2InputFileName) then begin

        ar2_read_ar2input, ar2InputFileName, ar2 = ar2

	    bT = ar2.bT 
	    bR = ar2.bR
	    bz = ar2.bz
        r  = ar2.r
        z  = ar2.z

        bMag = sqrt(br^2+bt^2+bz^2)
        dR = r[1]-r[0]
        dZ = z[1]-z[0]
        nR = n_elements(r)
        nZ = n_elements(z)

    endif else begin

        print, 'Please use either /ar2 or /eqdsk to choose a magnetic field'
        stop

    endelse

	;	Calculate the Guiding Center terms

	bT_B = bT / bMag
	bR_B = bR / bMag
	bz_B = bz / bMag

	omega   = q * bMag / m

	lnB = alog ( bMag )
	lnB_dR  = dlg_pDeriv ( lnB, 1, dR )
	lnB_dz  = dlg_pDeriv ( lnB, 2, dZ )

	bR_B_dz  	= dlg_pDeriv ( bR_B, 2, dZ )
	bT_B_dz  	= dlg_pDeriv ( bT_B, 2, dZ )
	bz_B_dR  	= dlg_pDeriv ( bz_B, 1, dR )
	bT_B_dR  	= dlg_pDeriv ( bT_B, 1, dR )
	bR_B_dR  	= dlg_pDeriv ( bR_B, 1, dR )
	bz_B_dz  	= dlg_pDeriv ( bz_B, 2, dZ )

	gradB_R 	= dlg_pDeriv ( bMag, 1, dR )
	gradB_z 	= dlg_pDeriv ( bMag, 2, dZ )

	bDotGradB   = bR_B * gradB_R + bz_B * gradB_z

	B_dR  = dlg_pDeriv ( bMag, 1, dR )
	B_dz  = dlg_pDeriv ( bMag, 2, dZ )
	
	;	Gradient Drift Terms
	;	--------------------
	;	Pick from either of these gradient term
	;	formulations, they both work.

	if not keyword_set ( grad ) then grad = 1

	if grad eq 1 then begin

		bxGradLnB_R =  1.0 / ( 2.0 * Omega ) * ( bT_B * lnB_dz )
		bxGradLnB_T = -1.0 / ( 2.0 * Omega ) * ( bR_B * lnB_dz - bz_B * lnB_dR )
		bxGradLnB_z = -1.0 / ( 2.0 * Omega ) * ( bT_B * lnB_dR )

		grad_R	= bxGradLnB_R
		grad_T	= bxGradLnB_T
		grad_z	= bxGradLnB_z

	endif else begin
	
		BxGradB_B2_R    = bT * B_dz / ( 2.0 * omega * bMag^2 )
		BxGradB_B2_T  = -1.0 * (bR * B_dz - bz * B_dR ) / ( 2.0 * omega * bMag^2 )
		BxGradB_B2_z    = -1.0 * bT * B_dR / ( 2.0 * omega * bMag^2 )
	
		grad_R	= BxGradB_B2_R
		grad_T	= BxGradB_B2_T
		grad_z	= BxGradB_B2_z

	endelse

	;	Curvature Drift Terms
	;	---------------------
	;	Select from any of the three formulations, they all
	;	work, i had to code them all to prove i had the right
	;	answer dammt ;-p

	if not keyword_set ( curv ) then curv = 1

	if curv eq 1 then begin

		bDotGradB_R     = bR_B * bR_B_dR + bz_B * bR_B_dz $
					- bT_B^2 / rebin(r,nR,nZ)  
		bDotGradB_T   = bT_B * bR_B / rebin(r,nR,nZ) $
					+ bR_B * bT_B_dR + bz_B * bT_B_dz
		bDotGradB_z     = br_B * bz_B_dR + bz_B * bz_B_dz

		bxbDotGradB_R =  1.0 / Omega * ( bT_B * bDotGradB_z - bz_B * bDotGradB_T )
		bxbDotGradB_T = -1.0 / Omega * ( br_B * bDotGradB_z - bz_B * bDotGradB_R )
		bxbDotGradB_z =  1.0 / Omega * ( br_B * bDotGradB_T - bT_B * bDotGradB_R ) 

		curv_R	= bxbDotGradB_R
		curv_T	= bxbDotGradB_T
		curv_z	= bxbDotGradB_z

	endif else if curv eq 2 then begin
	
		KR  = -1.0 * ( bT_B * ( bT_B_dR + 1.0 / rebin ( r, nR, nZ ) * bT_B ) $
				- bz_B * ( bR_B_dz - bz_B_dR ) )
		KT  = bR_B * ( bT_B_dR + 1.0 / rebin ( r, nR, nZ ) * bT_B ) $
				+ bz_B * bT_B_dz
		Kz  = -1.0 * ( bR_B * ( bR_B_dz - bz_B_dR ) + bT_B * bT_B_dz )
		
		bxK_R =  1.0 / Omega * ( bT_B * Kz - bz_B * KT )
		bxK_T = -1.0 / Omega * ( bR_B * Kz - bz_B * KR )
		bxK_z =  1.0 / Omega * ( bR_B * KT - bT_B * KR )
	
		curv_R	= bxK_R
		curv_T	= bxK_T
		curv_z	= bxK_z

	endif else begin
	
		gradPerpB_R =  B_dR - br_B^2 * B_dR - bz_B * br_B * B_dz
		gradPerpB_T = -br_B * bT_B * B_dR - bz_B * bT_B * B_dz
		gradPerpB_z =  B_dz - br_B * bz_B * B_dR - bz_B^2 * B_dz

		KR_   = gradPerpB_R / bMag
		KT_   = gradPerpB_T / bMag
		Kz_   = gradPerpB_z / bMag

		bxK_R_  =  1.0 / Omega * ( bT_B * Kz_ - bz_B * KT_ )
		bxK_T_  = -1.0 / Omega * ( bR_B * Kz_ - bz_B * KR_ )
		bxK_z_  =  1.0 / Omega * ( bR_B * KT_ - bT_B * KR_ )

		curv_R	= bxK_R_
		curv_T	= bxK_T_
		curv_z	= bxK_z_

	endelse

    if keyword_set(eqdsk) then $
	gc_truct = { unitb_R : bR_B, $
            unitb_phi : bT_B, $
            unitb_z : bz_B, $
            grad_R : grad_R, $
            grad_phi : grad_T, $
            grad_z : grad_z, $
            curv_R : curv_R, $
            curv_phi : curv_T, $
            curv_z : curv_z, $
            nW : nR, $
            nH : nZ, $
            rDim : eqdsk.rdim, $
            rLeft : eqdsk.rleft, $
            zDim : eqdsk.zdim, $
            bR : bR, $
            bPhi : bT, $
            bz : bz, $
            z : z, $
            R : r, $
            bDotGradB : bDotGradB, $
			mask : eqdsk.mask, $
			mi : m }


    if keyword_set(ar2InputFileName) then begin

        fName = 'gc_terms_e.nc'
	    nc_id = nCdf_create ( fName, /clobber )

	    nCdf_control, nc_id, /fill
	
        if keyword_set(oned) then begin

            slice = nZ/2
            z = z[slice]

            curv_r = curv_r[*,slice]
            curv_t = curv_t[*,slice]
            curv_z = curv_z[*,slice]

            grad_r = grad_r[*,slice]
            grad_t = grad_t[*,slice]
            grad_z = grad_z[*,slice]

            bDotGradB = bDotGradB[*,slice]

	        nR_id = nCdf_dimDef ( nc_id, 'nR', nR )
	        scalar_id = nCdf_dimDef ( nc_id, 'scalar', 1 )

            r_id = nCdf_varDef(nc_id, 'r', nr_id, /float)

	        curv_r_id = nCdf_varDef ( nc_id, 'curv_r', nr_id, /float )
	        curv_t_id = nCdf_varDef ( nc_id, 'curv_t', nr_id, /float )
	        curv_z_id = nCdf_varDef ( nc_id, 'curv_z', nr_id, /float )

	        grad_r_id = nCdf_varDef ( nc_id, 'grad_r', nr_id, /float )
	        grad_t_id = nCdf_varDef ( nc_id, 'grad_t', nr_id, /float )
	        grad_z_id = nCdf_varDef ( nc_id, 'grad_z', nr_id, /float )

            bDotGradB_id = nCdf_varDef ( nc_id, 'bDotGradB', nr_id, /float )

	        nCdf_control, nc_id, /enDef
	        
	        nCdf_varPut, nc_id, r_id, r
	        nCdf_varPut, nc_id, curv_r_id, curv_r 
	        nCdf_varPut, nc_id, curv_t_id, curv_t 
	        nCdf_varPut, nc_id, curv_z_id, curv_z 
	        nCdf_varPut, nc_id, grad_r_id, grad_r 
	        nCdf_varPut, nc_id, grad_t_id, grad_t 
	        nCdf_varPut, nc_id, grad_z_id, grad_z 
	        nCdf_varPut, nc_id, bDotGradB_id, bDotGradB

            nCdf_close, nc_id

        endif else begin

	        nR_id = nCdf_dimDef ( nc_id, 'nR', nR )
	        nZ_id = nCdf_dimDef ( nc_id, 'nZ', nZ )
	        scalar_id = nCdf_dimDef ( nc_id, 'scalar', 1 )

            r_id = nCdf_varDef(nc_id, 'r', nr_id, /float)
            z_id = nCdf_varDef(nc_id, 'z', nz_id, /float)

	        curv_r_id = nCdf_varDef ( nc_id, 'curv_r', [nr_id,nz_id], /float )
	        curv_t_id = nCdf_varDef ( nc_id, 'curv_t', [nr_id,nz_id], /float )
	        curv_z_id = nCdf_varDef ( nc_id, 'curv_z', [nr_id,nz_id], /float )

	        grad_r_id = nCdf_varDef ( nc_id, 'grad_r', [nr_id,nz_id], /float )
	        grad_t_id = nCdf_varDef ( nc_id, 'grad_t', [nr_id,nz_id], /float )
	        grad_z_id = nCdf_varDef ( nc_id, 'grad_z', [nr_id,nz_id], /float )

            bDotGradB_id = nCdf_varDef ( nc_id, 'bDotGradB', [nr_id,nz_id], /float )

	        nCdf_control, nc_id, /enDef
	        
	        nCdf_varPut, nc_id, r_id, r
	        nCdf_varPut, nc_id, z_id, z
	        nCdf_varPut, nc_id, curv_r_id, curv_r 
	        nCdf_varPut, nc_id, curv_t_id, curv_t 
	        nCdf_varPut, nc_id, curv_z_id, curv_z 
	        nCdf_varPut, nc_id, grad_r_id, grad_r 
	        nCdf_varPut, nc_id, grad_t_id, grad_t 
	        nCdf_varPut, nc_id, grad_z_id, grad_z 
	        nCdf_varPut, nc_id, bDotGradB_id, bDotGradB

            nCdf_close, nc_id

        endelse

    endif
    
    stop
end
