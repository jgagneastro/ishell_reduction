Function ishell_flat_fringing, flat_image, orders_structure, orders_mask, CORRECT_BLAZE_FUNCTION=correct_blaze_function, $
  LUMCORR_FLAT=lumcorr_flat, FRINGING_FLAT=fringing_flat, CORRECT_FRINGING=correct_fringing
  ;This program (eventually function) takes an image of a raw or median-combined flat field "flat_image",
  ; and an "order_structure" correpsonding to the order that is needed. It will correct the fringing and return a
  ; fringing-corrected flat field that has a similar structure to the input "flat_image"
  ; order_structure should be a single structure with the following keywords: 
  ;   height: Vertical height of order aperture
  ;   (x0,y0): Central position of the order (pixels)
  ;   (a0,a1,a2): Polynomial coefficients for the order trace
  ; CORRECT_BLAZE_FUNCTION: If this option is set, the blaze function and lamp spectra will be removed from the flat field.
  ; LUMCORR_FLAT_FILE: Returns a flat file where the spectral structure of the lamp and blaze function are corrected.
  ; FRINGING_FLAT: Returns a file containing the fringing pattern only
  ;
  ; It would be possible to make this code faster by skipping rows with only NaNs in horizontal_median.pro and probably elsewhere too.
  
  forward_function interpol2,horizontal_median
  
  if ~keyword_set(flat_image) or ~keyword_set(orders_structure) or ~keyword_set(orders_mask) then $
    message, ' You must input a flat image, an orders structure and an orders mask'
  
  ;Number of top/bottom pixels to cut from the fringe fitting region
  ; This is mostly temporary
  npix_cutoff_top = 3L
  npix_cutoff_bottom = 5L
  
  ;The is the cutoff for fraction flux below which the fringing won't be fitted,
  ; when looking at the order flux in the vertical direction
  flux_cutoff = 0.1
  ;This is the quartile fraction used to normalize flux
  medval_cutoff = .98
  
  ;The is the cutoff for fraction flux below which the fringing won't be fitted,
  ; when looking at the lamp spectrum
  spectral_flux_cutoff = 0.5
  
  ;This is the horizontal smoothing that is used to de-couple fringing from 
  ; the spectral structure of the lamp (in pixels)
  nhsmooth = 45L
  
  ;Box size for smoothing the fringe pattern (this should be close to the Nyquist sampling)
  fringe_nsmooth = 3L
  
  ;This is the lowest allowed fraction of flux in any pixel of the flat field. Otherwise they will be masked
  min_allowed_flux = .3
  
  ;Number of rows/columns to be masked at the edges of the flat field
  nmask_bottom_rows = 4L
  nmask_top_rows = 24L
  nmask_right_cols = 4L
  nmask_left_cols = 4L
  
  ;Read image size
  nx = (size(flat_image))[1]
  ny = (size(flat_image))[2]
  xval = dindgen(nx)#make_array(ny,value=1d0,/double)
  yval = make_array(nx,value=1d0,/double)#dindgen(ny)
  
  ;Create a final flat image
  final_flat = dblarr(nx,ny)+!values.d_nan
  lumcorr_flat = dblarr(nx,ny)+!values.d_nan
  fringing_flat = dblarr(nx,ny)+!values.d_nan
  
  ;Number of orders to be parsed
  n_orders = n_elements(orders_structure)
  
  ;Loop on orders
  for i=0L, n_orders-1L do begin
    
    print, ' Correcting flat field order #['+strtrim(i+1L,2)+'/'+strtrim(n_orders,2L)+'] ...'
    
    ;Height of order
    height = ceil(orders_structure[i].height)
    
    ;Create an image of this order only by masking everything outside of it
    order_image = flat_image
    
    ;Mask everything outside the order
    order_left_location = poly(xval,orders_structure[i].left_coeffs)
    order_right_location = poly(xval,orders_structure[i].right_coeffs)
    bad = where(yval lt order_left_location or yval gt order_right_location, nbad)
    if nbad ne 0L then order_image[bad] = !values.d_nan
    
    ;Create a straightened version of the flat in this order
    straight_flat_order = dblarr(nx,height)+!values.d_nan
    xarr = dindgen(nx)
    yarr = dindgen(ny)
    
    ;Actually do the straightening of the flat
    order_center_location = poly(xarr,orders_structure[i].mid_coeffs)
    subyarr = dindgen(height)
    for l=0L, nx-1L do begin
      ;Skip this column if the order is partially out of frame
      if (order_center_location[l]-height/2d0) le 0 then continue
      straight_flat_order[l,*] = interpol2(order_image[l,*],yarr,(order_center_location[l]-height/2d0)+subyarr,/repairnans)
    endfor
    
    ;Mask the regions outside of the trace
    hmedian_flat = median(straight_flat_order,dim=1)
    flat_border_cutoff = 0.65
    bad = where(hmedian_flat le max(hmedian_flat,/nan)*flat_border_cutoff, nbad)
    if nbad ne 0L then straight_flat_order[*,bad] = !values.d_nan
    
    ;Count how many pixels were removed on each side
    gleft = where(bad lt height/2L, ngleft)
    gright = where(bad gt height/2L, ngright)
    
    ;Mask these lines in the original flat image
    bad = where(yval lt (order_left_location+ngleft) or yval gt (order_right_location-ngright), nbad)
    if nbad ne 0L then order_image[bad] = !values.d_nan
    
    ;Create a horizontally smoothed version of the flat to bring out the lamp spectrum
    flat_hsmooth = horizontal_median(straight_flat_order,nhsmooth)
    
    ;If required mask more pixels in the vertical direction
    spatial_profile = median(flat_hsmooth,dim=1)
    gfin = where(finite(spatial_profile), ngfin)
    if npix_cutoff_top gt 1L then $
      flat_hsmooth[*,gfin[0L:npix_cutoff_top-1L]] = !values.d_nan
    if npix_cutoff_top eq 1L then $
      flat_hsmooth[*,gfin[0L]] = !values.d_nan
    if npix_cutoff_bottom gt 1L then $
      flat_hsmooth[*,gfin[-npix_cutoff_bottom:*]] = !values.d_nan
    if npix_cutoff_bottom eq 1L then $
      flat_hsmooth[*,gfin[-1L]] = !values.d_nan
    
    ;Recreate a better spectral profile, if needed
    spectral_profile = median(flat_hsmooth,dim=2)
    
    ;Create an image where fringing is more obvious
    fringing_image = straight_flat_order/flat_hsmooth
    
    ;Create a 1D version of the fringing
    fringing_1d = median(fringing_image,dim=2)
    
    ;Create a smoothed version
    fringe_smooth = median(fringing_1d,fringe_nsmooth)
    
    ;Define a 1D correction to be applied to the flat field
    correct_1d = dblarr(nx)+1d0
    
    if keyword_set(correct_blaze_function) then $
      correct_1d *= spectral_profile
    if keyword_set(correct_fringing) then $
      correct_1d *= fringe_smooth
    
    ;Replicate that correction pattern in the vertical direction
    fringing_spectrum_rep = correct_1d#make_array(ny,value=1d0,/double)
    
    ;Correct the fringing in the original curved flat
    flatfield_corrected = order_image/fringing_spectrum_rep
    
    ;Patch the corrected order back into the final flat
    g_within_order = where(orders_mask eq orders_structure[i].order_id, ng_within_order)
    if ng_within_order eq 0L then $
      message, 'No order positions were found !'
    final_flat[g_within_order] = flatfield_corrected[g_within_order]
    
    ;Create a fringing flat
    fringing_rep = fringe_smooth#make_array(ny,value=1d0,/double)
    fringing_flat[g_within_order] = fringing_rep[g_within_order]
    
    ;Create a temporary version of a fully corrected flat field to mask pixels that are too dark
    fringing_specres_rep = (fringe_smooth*spectral_profile)#make_array(ny,value=1d0,/double)
    fullcor_flat = dblarr(nx,ny)+!values.d_nan
    fullcor_flat[g_within_order] = order_image[g_within_order] / fringing_specres_rep[g_within_order]
    bad = where(fullcor_flat[g_within_order] le min_allowed_flux, nbad)
    if nbad ne 0L then begin
      final_flat[g_within_order[bad]] = !values.d_nan
      fullcor_flat[g_within_order[bad]] = !values.d_nan
      fringing_flat[g_within_order[bad]] = !values.d_nan
    endif
    
    lumcorr_flat[g_within_order] = fullcor_flat[g_within_order]
    
  endfor
  
  ;Remove edge data from the flat
  ;Mask the flat edges
  if keyword_set(nmask_bottom_rows) then $
    final_flat[*,0:nmask_bottom_rows-1L] = !values.d_nan
  if keyword_set(nmask_top_rows) then $
    final_flat[*,-nmask_top_rows:*] = !values.d_nan
  if keyword_set(nmask_left_cols) then $
    final_flat[0:nmask_left_cols-1L,*] = !values.d_nan
  if keyword_set(nmask_right_cols) then $
    final_flat[-nmask_right_cols:*,*] = !values.d_nan
  
  ;Pass the final flat to main
  return, final_flat
  
End