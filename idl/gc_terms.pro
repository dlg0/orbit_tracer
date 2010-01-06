pro gc_terms, mi, eqdsk, $
	gc_struct = bStruct, $
	curv = curv, $
	grad = grad

	q	= 1.602e-19

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



end
