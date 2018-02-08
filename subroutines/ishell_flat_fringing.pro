Function ishell_flat_fringing, flat_image, orders_structure, orders_mask, CORRECT_BLAZE_FUNCTION=correct_blaze_function, $
  LUMCORR_FLAT=lumcorr_flat, FRINGING_FLAT=fringing_flat, CORRECT_FRINGING=correct_fringing, $
  FRINGING_SOLUTION_1D=fringing_solution_1d, FRINGE_NSMOOTH=fringe_nsmooth, MODEL_FRINGING=model_fringing, $
  DETECTOR_PATTERNS=detector_patterns, MODELS=models, FLAT_ILLUMINATION=flat_illumination
  ;This function takes an image of a raw or median-combined flat field "flat_image",
  ; and an "order_structure" correpsonding to the order that is needed. It will correct the fringing and return a
  ; fringing-corrected flat field that has a similar structure to the input "flat_image"
  ; order_structure should be a single structure with the following keywords: 
  ;   height: Vertical height of order aperture
  ;   (x0,y0): Central position of the order (pixels)
  ;   (a0,a1,a2): Polynomial coefficients for the order trace
  ; CORRECT_BLAZE_FUNCTION: If this option is set, the blaze function and lamp spectra will be removed from the flat field.
  ; LUMCORR_FLAT_FILE: Returns a flat file where the spectral structure of the lamp and blaze function are corrected.
  ; FRINGING_FLAT: Returns a file containing the fringing pattern only
  ; FRINGING_SOLUTION_1D: Returns a 1D median fringing in each order
  ; MODEL_FRINGING: Whether or not to model fringing with a period-dependent sine model to better remove it from flats
  ; MODELS: Output 2D model of fringing (only if MODEL_FRINGING=1)
  ; DETECTOR_PATTERNS: Output 1D array with detector patterns (only if MODEL_FRINGING=1)
  ; FLAT_ILLUMINATION: Spatial illumination of the flats measured with a very wide median filter.
  ;   This includes the spectrum of the lamp and the Blaze function of each order.
  ;  
  ;
  ; It would be possible to make this code faster by skipping rows with only NaNs in horizontal_median.pro and probably elsewhere too.
  
  forward_function interpol2, horizontal_median, ishell_fringing_1d_model
  
  if ~keyword_set(flat_image) or ~keyword_set(orders_structure) or ~keyword_set(orders_mask) then $
    message, ' You must input a flat image, an orders structure and an orders mask'
  
  ;Number of top/bottom pixels to cut from the fringe fitting region
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
  if ~keyword_set(fringe_nsmooth) then $
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
  
  ;Find out the number of orders for that particular object
  good_orders = where(orders_structure.order_id ge 0, n_orders)
  fringing_solution_1d = dblarr(nx,n_orders)+!values.d_nan
  flat_illumination = dblarr(nx,n_orders)+!values.d_nan
  
  ;Loop on orders
  for i=0L, n_orders-1L do begin
    
    print, '  Bringing out fringing in flat field order #['+strtrim(i+1L,2)+'/'+strtrim(n_orders,2L)+'] ...'
    
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
    for l=0L, nx-1L do begin & $
      ;Skip this column if the order is partially out of frame
      if (order_center_location[l]-height/2d0) le 0 then continue & $
      straight_flat_order[l,*] = interpol2(order_image[l,*],yarr,(order_center_location[l]-height/2d0)+subyarr,/repairnans) & $
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
    
    ;Apply additional smoothing to remove any fringing residuals
    spectral_profile = smooth_error(median(spectral_profile,nhsmooth),nhsmooth)
    
    flat_illumination[*,i] = spectral_profile
    
    ;Create an image where fringing is more obvious
    fringing_image = straight_flat_order/flat_hsmooth
    
    ;Create a 1D version of the fringing
    fringing_1d = median(fringing_image,dim=2)
    
    ;If fringing will be modeled, horizontal smooth will be applied at the very end
    if keyword_set(model_fringing) then begin
      fringe_smooth = fringing_1d
    endif else begin
      ;Create a horizontally smoothed version of the fringing
      if fringe_nsmooth gt 1 then $
        fringe_smooth = median(fringing_1d,fringe_nsmooth) $
      else $
        fringe_smooth = fringing_1d
    endelse
    
    ;Store in the 1D fringe correction
    fringing_solution_1d[*,i] = fringe_smooth
    
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
    
    ;Create a fringing flat (if model_fringing = 1 this will be done later)
    if model_fringing eq 0 then begin
      fringing_rep = fringe_smooth#make_array(ny,value=1d0,/double)
      fringing_flat[g_within_order] = fringing_rep[g_within_order]
    endif
    
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
  
  ;If this keyword is set, fringing will be modeled with a 1-component variable-period sine
  ; This model will then be used to separate the fringing from other column-dependent detector
  ; patterns, and will leave the detector patterns in the flat fields. 
  if keyword_set(model_fringing) then begin
    ;nrowi = n_orders
    
    ;Number of parameters in the fringing model
    n_parameters = 4L
    
    ;Number of "walkers" that will independently run a Levenberg-Marquardt least-squares fitting
    ; A large number requires more CPU but will be more robust against local minima
    nfit = 40
    
    ;Initiate arrays that will contain the best-fitting parameters and the models
    model_pars = dblarr(4L,n_orders)+!values.d_nan
    models = dblarr(nx,n_orders)+!values.d_nan
    
    ;Pixel column indices on which to perform the fit
    min_fit_index = 400L
    max_fit_index = 1500L
    
    ;Defaut parameter estimates
    fit_par_estim = [0d0,$;Amplitude
      36.5d0,$;Period (pixels)
      0d0,$;Phase (radian)
      -0.0015d0];Period slope (pixel period / pixel = no units)
    
    ;Typical variations to be explored by the N walkers in parameter space
    fit_par_scatter = [0d0,$;Amplitude
      5d0,$;Period (pixels)
      !dpi/10d0,$;Phase (radian)
      1d-4];Period slope (pixel period / pixel = no units)
    
    ;Fit fringing in each order
    for oi=0L, n_orders-1L do begin
      
      ;Data to be fitted with the fringing model
      print, '  Modeling fringing in flat field order #['+strtrim(oi+1L,2)+'/'+strtrim(n_orders,2L)+'] ...'
      
      fit_y = fringing_solution_1d[min_fit_index:max_fit_index,oi]
      fit_x = dindgen(n_elements(fit_y))+min_fit_index
      
      ;Adjust estimated amplitude
      fit_par_estim[0] = weighted_median(abs(fit_y-1d0),medval=.9)
      
      ;Adjust amplitude variations to be explored
      fit_par_scatter[0] = fit_par_estim/1d2
      
      ;Create random initial positions for the walkers 
      fit_par_noise = fit_par_estim#make_array(nfit,value=1d0,/double) + (fit_par_scatter#make_array(nfit,value=1d0,/double))*randomn(seed,n_elements(fit_par_estim),nfit)
      fit_par_noise[*,0] = fit_par_estim
      
      ;Arrays to store reduced chi squares and parameter values
      redchi2s = dblarr(nfit)+!values.d_nan & $
      fitpars = dblarr(n_elements(fit_par_estim),nfit)+!values.d_nan & $
      
      ;Perform least-squares fitting for each walker
      for fiti=0L, nfit-1L do begin
        fit_pari = mpfitfun('ishell_fringing_1d_model',fit_x,fit_y,1d0,fit_par_noise[*,fiti],YFIT=yfit,status=status,err=err,/nan,/quiet)
        ;Calculate reduced chi2 using only finite pixels
        redchi2s[fiti] = total((fit_y-yfit)^2,/nan)/double(total(finite(fit_y-yfit)))
        fitpars[*,fiti] = fit_pari
      endfor
      
      ;Idenfity the best (minimum) chi2 and select that solution
      void = min(redchi2s,wmin)
      fit_par = fitpars[*,wmin]
      
      ;Optional: Display best fit
      ;plot,fit_x,fit_y, xtitle=strtrim(oi,2) & oplot, fit_x, ishell_fringing_1d_model(fit_x,fit_par), col=255
      
      ;Store best fit parameters and model
      model_pars[*,oi] = fit_par
      models[*,oi] = ishell_fringing_1d_model(dindgen(nx),fit_par)
    endfor
    
    ;Transform negative amplitudes to a pi phase shift
    gneg = where(reform(model_pars[0,*]) lt 0, ngneg)
    if ngneg ne 0L then model_pars[0,gneg] *= -1
    if ngneg ne 0L then model_pars[2,gneg] += !dpi
    model_pars[2,*] = ((model_pars[2,*]+2*!dpi) mod (!dpi*2d0))
    
    ;Take an extremely horizontally smoothed version of the observed fringing to bring out Blaze function of this order   
    nsmooth_fringe_blaze = 100L
    fringe_blaze_function = smooth(median(max(fringing_solution_1d,dim=2,/nan),nsmooth_fringe_blaze),nsmooth_fringe_blaze)
    
    ;Normalize the Blaze function appropriately and mask the edges
    fringe_blaze_function = (fringe_blaze_function-1.)/max((fringe_blaze_function-1.),/nan)
    fringe_blaze_function[0:nsmooth_fringe_blaze] = 0.
    fringe_blaze_function[-nsmooth_fringe_blaze:*] = 0.
    
    ;Create a model that includes the Blaze function
    model_with_edges = (fringe_blaze_function#make_array(n_orders,value=1d0,/double))*(models-1.)+1.
    
    ;Bring out detector patterns by dividing observed fringing by fringing model
    ; and then taking a vertical median to average out any fringing model residual
    ; We are seeking column-dependent detector patterns so those will survive a vertical median
    ; (observed fringing contains not only true fringing but also detector patterns) 
    detector_patterns = median(fringing_solution_1d/model_with_edges,dim=2)
    
    ;Optional: Figure displaying detector patterns and detector "breaks" at 64-pixels steps
    ;plot,detector_patterns-1,yrange=[-.03,.04] & for i=0L, 32L do oplot, i*64+[0,0], [-10,10], col=255 & oplot,detector_patterns-1
    
    ;Create a 2D version of the detector patterns
    detector_patterns_2d = (detector_patterns#make_array(n_orders,value=1d0,/double))
    
    ;Smooth the fringing solution w/o detector patterns to eliminate any pixel-to-pixel
    ; sensitivity variations which should remain in the flat field
    fringing_no_detector_patterns = fringing_solution_1d/detector_patterns_2d
    for sri=0L, n_orders-1L do $
      fringing_no_detector_patterns[*,sri] = median(fringing_no_detector_patterns[*,sri],fringe_nsmooth)
    
    ;Recreate 2D image of fringing in flat
    fringing_no_detector_patterns_2d = (fringing_no_detector_patterns[*,sri]#make_array(ny,value=1d0,/double))
    for sri=0L, n_orders-1L do begin & $
      g_within_order = where(orders_mask eq orders_structure[sri].order_id, ng_within_order) & $
      if ng_within_order eq 0L then $
        message, 'No order positions were found !' & $
      fringing_flat[g_within_order] = fringing_no_detector_patterns_2d[g_within_order] & $
    endfor
    
    ;Remove detector patterns from flat and fringing solutions
    fringing_solution_1d = fringing_no_detector_patterns
    detector_patterns_ny = detector_patterns#make_array(ny,value=1d0,/double)
    fringing_flat /= detector_patterns_ny
    
    ;Put back the detector patterns in flat fields
    final_flat *= detector_patterns_ny
    lumcorr_flat *= detector_patterns_ny
    
  endif
  
  ;Pass the final flat to main
  return, final_flat
  
End