Pro ishell_extract_exposure, fits_data, tcs_obj, object_names, integration_times, $; Information from observing log
    darks_ids_uniq, flat_ids_uniq, n_orders, g_science, ng_science, sci, science_flat_ids, nx, ny, $; Other needed variables
    edge_npix_mask, edge_npix_avoid, sizelag, window_lag_poly_fit, ndegree_sigma_fit, ndegree_poly_fit, nrows_sky, quartile_fraction_for_norm, $; Adjustable parameters Pt.1
    bad_pixels_threshold, second_bad_pixels_threshold, read_noise, dark_current, gain, $; Adjustable parameters Pt.2
    generate_full_orders_spectra_figures,generate_residuals_fits_images,generate_full_orders_trace_profile_figures, $; Adjustable parameters Pt.3
    spectra_dir, output_dir, trace_profiles_dir, previews_dir, $; Paths
    correct_blaze_function_in_flatfield,correct_fringing_in_flatfield,model_fringing,code_version,remove_detector_patterns_from_data,override_flats, $;Header information
    flats_uniq_cube_corrected, orders_mask_cube, $; Needed data
    order_ids, order_heights, order_left_coeffs, order_mid_coeffs, $; Orders structure
    DO_DARK_SUBTRACTION=do_dark_subtraction, DO_TRACE_2D_FIT=do_trace_2d_fit, MASK_TRACE_EDGES=mask_trace_edges, REFINE_POLY_DROP_PIXELS=refine_poly_drop_pixels, $; Keywords Pt.1
    MODEL_REFINEMENT_CURVATURE=model_refinement_curvature, GENERATE_INDIVIDUAL_ORDERS_TRACE_PROFILE_FIGURES=generate_individual_orders_trace_profile_figures, $; Keywords Pt.2
    GENERATE_INDIVIDUAL_ORDERS_SPECTRA_FIGURES=generate_individual_orders_spectra_figures; Keywords Pt.3

  ;Select data file
  data_file = fits_data[g_science[sci]]

  ;Select the proper flat field and order structure
  g_flat_id = where(strlowcase(flat_ids_uniq) eq strlowcase(science_flat_ids[sci]), ngflat_id)
  if ngflat_id eq 0L then begin
    message, ' No flat ID corresponding to object '+object_names[g_science[sci]]+' (TCS_OBJ = '+tcs_obj[g_science[sci]]+') was found, using generic flat field if possible...'
    g_flat_id = where(flat_ids_uniq eq 'BLANK', ngflat_id)
    if ngflat_id eq 0L then $
      message, ' No flat ID corresponding to object '+object_names[g_science[sci]]+' (TCS_OBJ = '+tcs_obj[g_science[sci]]+') was found, and no generic flat was found !'
  endif

  flat_field_corrected = flats_uniq_cube_corrected[*,*,g_flat_id[0L]]
  orders_mask = orders_mask_cube[*,*,g_flat_id[0L]]
  
  ;Read correct order information for that exposure
  ord_heights = order_heights[*,g_flat_id[0L]]
  ord_ids = order_ids[*,g_flat_id[0L]]
  ord_left_coeffs = order_left_coeffs[*,*,g_flat_id[0L]]
  ord_mid_coeffs = order_mid_coeffs[*,*,g_flat_id[0L]]
  
  ;Determine the number of orders
  good_orders = where(ord_ids ge 0, n_orders)

  ;If all orders of this exposure are already present, skip this exposure altogether
  object_name = object_names[g_science[sci]]
  output_files = spectra_dir+file_basename(data_file,file_ext(data_file,/FITS))+'_OBJ_'+object_name+'_ORD_'+strtrim(ord_ids[good_orders],2)+'_spectrum.txt'
  if min(file_test(output_files)) eq 1 then begin
    print, '   Skipping extraction because it already exists...'
    return
  endif

  ;Select the proper combined dark exposure
  if do_dark_subtraction eq 1 then begin
    g_dark_id = where(strlowcase(darks_ids_uniq) eq rtrim(integration_times[g_science[sci]],3), ngdark_id)
    if ngdark_id eq 0L then $
      message, ' No dark ID corresponding to exposure time '+rtrim(integration_times[g_science[sci]],3)+'s was found, skipping the dark subtraction !', /continue
    dark_im = double(readfits(fits_data[g_dark_id[0L]]))
  endif else $
    dark_im = dblarr(nx,ny)

  ;Read the data file
  data_im = double(readfits(data_file,header,/silent))

  ;Mask edge pixels
  data_im[0L:edge_npix_mask-1L,*] = !values.d_nan
  data_im[-edge_npix_mask:*,*] = !values.d_nan

  ;Subtract darks
  data_im -= dark_im

  ;Flat-field the data
  data_im /= flat_field_corrected

  ;Create an X positions array
  nx = (size(data_im))[1]
  xarr = dindgen(nx)
  ny = (size(data_im))[2]
  yarr = dindgen(ny)

  ;Initiate necessary variables for the 2D model traces mode
  if keyword_set(do_trace_2d_fit) then begin
    ;Create directories for models and residuals
    models_dir = output_dir+'models'+path_sep()
    if ~file_test(models_dir) then file_mkdir, models_dir
    residuals_dir = output_dir+'residuals'+path_sep()
    if ~file_test(residuals_dir) then file_mkdir, residuals_dir
    ;Create an array to store the 2D models if they are gonna be created
    full_2d_model = dblarr(nx,ny)+!values.d_nan
    ;Initiate an array to store the ytrace position and seeing vs wavelength polynomials
    full_model_polynomials = dblarr(ndegree_poly_fit+ndegree_sigma_fit,n_orders)+!values.d_nan
  endif

  ;Create necessary directory for intermediate IDL files used for full-orders figures
  extraction_idl_dir = output_dir+'extraction_idl_data'+path_sep()
  if ~file_test(extraction_idl_dir) then file_mkdir, extraction_idl_dir

  ;Initiate storage arrays that will be useful for diagnostic plots
  full_orders_spectra = dblarr(nx,n_orders)+!values.d_nan
  full_orders_raw_spectra = dblarr(nx,n_orders)+!values.d_nan
  ;Trace profiles are stored as a pointer array rather than a 2D array, because they might
  ; have different sizes (not all orders have the same height)
  full_orders_trace_profile = ptrarr(n_orders,/allocate)
  full_orders_trace_xfit = ptrarr(n_orders,/allocate)
  full_orders_trace_yfit = ptrarr(n_orders,/allocate)

  ;Do a loop over orders
  ;NOTE: A split_for could go here if needed
  for i=0L, n_orders-1L do begin

    print, '  Extracting Orders ['+strtrim(i+1L,2L)+'/'+strtrim(n_orders,2L)+'] of Exposure ['+strtrim(sci+1L,2L)+'/'+strtrim(ng_science,2L)+'] ('+object_names[g_science[sci]]+') !'

    ;Skip extraction if output data already exists
    strkill, object_name, ' '
    output_file = output_files[i]
    if file_test(output_file) then begin
      print, '   Skipping extraction because it already exists...'

      ;If the partial IDL files are present, restore them to allow the code to plot
      ; full-orders figures
      extraction_data_file = extraction_idl_dir+'extraction_idl_data'+strtrim(ord_ids[good_orders[i]],2L)+'_'+object_name+'_'+file_basename(data_file,file_ext(data_file))+'.sav'
      if file_test(extraction_data_file) then begin
        ;Restore the necessary data
        opt_spectrum_norm = !NULL
        raw_spectrum_norm = !NULL
        trace_profile = !NULL
        ytrace_precise = !NULL
        gauss_fit_precise = !NULL
        restore, extraction_data_file;, opt_spectrum_norm, raw_spectrum_norm, trace_profile, ytrace_precise, gauss_fit_precise
        ;Make sure that the right data was restored
        if keyword_set(opt_spectrum_norm) and keyword_set(raw_spectrum_norm) and $
          keyword_set(trace_profile) and keyword_set(ytrace_precise) and $
          keyword_set(gauss_fit_precise) then begin

          ;Store the partial IDL results in full-orders arrays
          full_orders_spectra[*,i] = opt_spectrum_norm
          full_orders_raw_spectra[*,i] = raw_spectrum_norm
          *full_orders_trace_profile[i] = trace_profile
          *full_orders_trace_xfit[i] = ytrace_precise
          *full_orders_trace_yfit[i] = gauss_fit_precise
        endif

      endif

      ;Restore the partial 2D model if it exists and that the rest is to be created
      if keyword_set(do_trace_2d_fit) then begin
        ;Make sure the save file exists
        models_file = models_dir+'models_order'+strtrim(ord_ids[good_orders[i]],2L)+'_'+object_name+'_'+file_basename(data_file,file_ext(data_file))+'.sav'
        if file_test(models_file) then begin
          ;Make sure the save file contains the updated 2d model
          updated_2d_model = !NULL
          restore, models_file;, best_fit_parameters, functargs, updated_2d_model, badpix_mask
          if keyword_set(best_fit_parameters) then $
            full_model_polynomials[*,i] = best_fit_parameters
          if keyword_set(updated_2d_model) then begin
            ;Identify order positions and patch the 2D model at the order position
            g_order = where(orders_mask eq ord_ids[good_orders[i]], ng_order)
            full_2d_model[g_order] = updated_2d_model[g_order]
          endif
        endif
      endif
      continue
    endif

    ;Create a list of y positions for that orders trace
    y_position_order_i_estimate = poly(xarr,ord_mid_coeffs[*,good_orders[i]])

    ;Define the height of the order
    height = ceil(ord_heights[good_orders[i]])

    ;Create an image where only the relevant order is seen
    im_order = data_im
    g_order = where(orders_mask eq ord_ids[good_orders[i]], ng_order, complement=bad_order, ncomplement=nbad_order)
    if nbad_order ne 0L then im_order[bad_order] = !values.d_nan

    ;Straighten the data trace with that estimated polynomial
    straight_data_order = dblarr(nx,height)+!values.d_nan
    subyarr = dindgen(height)
    for l=0L, nx-1L do begin
      ;Skip this column if the order is partially out of frame
      if (y_position_order_i_estimate[l]-height/2d0) le 0 then continue
      straight_data_order[l,*] = interpol2(reform(im_order[l,*]),yarr,(y_position_order_i_estimate[l]-double(height)/2d0)+subyarr,/repairnans,badvalue=!values.d_nan)
    endfor

    ;Do a horizontal (spectral direction) median crunch to get a profile of the trace
    trace_profile = median(straight_data_order,dim=1)

    ;Mask the edges of the trace
    if keyword_set(mask_trace_edges) then begin
      trace_profile[0L:mask_trace_edges-1L] = !values.d_nan
      trace_profile[-mask_trace_edges:*] = !values.d_nan
    endif

    ;Re-scale the trace profile
    trace_profile -= weighted_median(trace_profile,medval=.1)
    trace_profile >= 0d0
    trace_profile /= max(trace_profile,/nan)

    ;Remove any NaNs in the trace profile
    bad = where(~finite(trace_profile), nbad)
    if nbad ne 0L then trace_profile[bad] = 0d0

    ;Determine where the trace profile is at less than 5% of its maximum
    max_profile = max(trace_profile,wmax,/nan)
    ytrace = lindgen(height)
    right_glow_flux = where(trace_profile le .05 and ytrace gt wmax, nglow_flux_right)
    right_cutoff = ytrace[min(right_glow_flux,/nan)]
    left_glow_flux = where(trace_profile le .05 and ytrace lt wmax, nglow_flux_left)
    left_cutoff = ytrace[max(left_glow_flux,/nan)]

    ;Remove flux everywhere left and right of these points
    trace_profile[0:left_cutoff] = 0d0
    trace_profile[right_cutoff:*] = 0d0

    if total(finite(trace_profile)) eq 0. then $
      message, ' Could not properly build a trace profile ! (The spectrum is all NaNs, make sure the trace profile is OK)'

    ;Use a gaussian fit to precisely determine the center position of the profile
    yfit_gauss = gaussfit2(ytrace, trace_profile, gauss_parameters,ESTIMATES=[1d0,double(wmax),double(abs(right_cutoff-left_cutoff))*3d0/20d0],NTERMS=3L)
    trace_max_pos = gauss_parameters[1L]

    ;Create and save a figure of the trace profile compared to a Gaussian fit (this is useful for debugging)
    if ~file_test(trace_profiles_dir) then file_mkdir, trace_profiles_dir
    trace_profiles_file = trace_profiles_dir+'trace_profile_order'+strtrim(ord_ids[good_orders[i]],2L)+'_'+object_name+'_'+file_basename(data_file,file_ext(data_file))+'.png'
    ytrace_precise = dindgen(1d3)/double(1d3-1)*double(height)
    gauss_funct, ytrace_precise, gauss_parameters, gauss_fit_precise

    ;Make a plot of the trace profile
    if generate_individual_orders_trace_profile_figures then begin
      pa = plot(buffer=1,ytrace,trace_profile, thick=2, name='Trace Profile', yticklen=.015, xticklen=.015, xthick=2, ythick=2, margin=[.12,.12,.02,.02],xtitle='Pixel',ytitle='Relative Flux')
      pb = plot(ytrace_precise,gauss_fit_precise, thick=3, name='Gaussian Fit', /overplot, color='red')
      pb.order, /send_to_back
      pl = legend(target=[pa,pb],position=[.95,.95],horizontal_align=1,vertical_align=1,font_size=10,transparency=50,shadow=0,thick=2)
      pa.save, trace_profiles_file
      pa.close
    endif

    ;Smooth the data horizontally to align it better
    smoothed_data_order = horizontal_median(straight_data_order,20)

    ;Do a first rough extraction from the straightened trace
    trace_profile_straight_2d = make_array(nx,value=1d0,/double)#trace_profile
    spectrum_initial_estimate = total(smoothed_data_order*trace_profile_straight_2d,2L,/nan)/total(trace_profile_straight_2d,2L,/nan)
    if total(finite(spectrum_initial_estimate)) eq 0. then $
      message, ' Could not properly build an initial spectrum estimate ! (The spectrum is all NaNs, make sure the trace profile is OK)'

    ;Create an array that only covers the center of the order
    gfinite_estim = where(finite(spectrum_initial_estimate),nfinite_estim)
    if nfinite_estim-2*edge_npix_avoid lt 500 then $
      message, ' Removing order edges leaves too few data points for the code to run properly !'
    order_center_positions = gfinite_estim[edge_npix_avoid:-edge_npix_avoid-1]

    ;Normalize the intial spectrum
    spectrum_initial_estimate /= weighted_median(spectrum_initial_estimate[order_center_positions],medval=0.8)

    ;Create a vector of lags for the cross-correlation
    if ~keyword_set(lag) then begin
      lag = dindgen(sizelag)+1d0
      lag = [-reverse(lag),0.,lag]
    endif

    ;Loop on pixel position to refine the trace position
    trace_pos_deviations = dblarr(nx)+!values.d_nan
    for l=0L, nx-1L do begin

      ;Skip pixels where the intial spectrum has a flux of less than 50% of the max
      if ~finite(spectrum_initial_estimate[l]) or spectrum_initial_estimate[l] le 0.5 then continue

      ;Take a slice of the data at that pixel position and normalize it
      data_slice = smoothed_data_order[l,*]
      data_slice /= weighted_median(data_slice,medval=quartile_fraction_for_norm)
      data_slice >= 0d0

      ;Do the cross-correlation
      ccs = fc_correlate(data_slice,trace_profile,lag,/INTEGERS)

      ;Give more importance to cross-correlation values near the center of the array
      ccs *= exp(-(dindgen(n_elements(ccs))-double(height)/2d0)^2/(2.*sizelag^2)*3)

      ;Fit a second-order polynomial near the center of the cross-correlation function
      ; to precisely determine its central value
      void = max(ccs,wmax_ccs,/nan)
      if wmax_ccs gt n_elements(lag)-window_lag_poly_fit-1L or wmax_ccs lt window_lag_poly_fit then continue
      ;message,'Fix illegal subscript range CCS'
      fit_ccs_par = poly_fit(lag[wmax_ccs-window_lag_poly_fit:wmax_ccs+window_lag_poly_fit],ccs[wmax_ccs-window_lag_poly_fit:wmax_ccs+window_lag_poly_fit],2)

      ;Analytically determine the maximum value of this second-order polynomial.
      ; This is the refined position of the trace at this column
      trace_pos_deviations[l] = 0.5d0*fit_ccs_par[1]/fit_ccs_par[2]

    endfor

    ;Drop the last ~5 points on both ends
    if keyword_set(refine_poly_drop_pixels) then begin
      g_finite_devs = where(finite(trace_pos_deviations), ng_finite_devs)
      if ng_finite_devs gt 2L*refine_poly_drop_pixels then begin
        trace_pos_deviations[g_finite_devs[0L:refine_poly_drop_pixels-1L]] = !values.d_nan
        trace_pos_deviations[g_finite_devs[-refine_poly_drop_pixels:*]] = !values.d_nan
      endif
    endif

    ;Create a refined list of trace positions
    y_position_order_i_refined = y_position_order_i_estimate+trace_pos_deviations

    ;Refine the coefficients with a polynomial fit
    gfin_pos_refined = where(finite(y_position_order_i_refined), ngfin_pos_refined)
    if ngfin_pos_refined le 10 then begin
      message, ' The refined trace position has too few finite data points !', /continue
      continue
    endif
    refined_coefficients = reform(poly_fit(xarr[gfin_pos_refined],y_position_order_i_refined[gfin_pos_refined],ndegree_poly_fit-1L))

    ;Shift the first coefficient so that the profile matches the trace exactly
    ; this shift is non-zero only when the trace is not exactly at the center of the order
    refined_coefficients[0L] += trace_max_pos-double(height)/2d0

    ;Recreate an array of refined trace positions (with no NaNs)
    y_position_refined = poly(xarr,refined_coefficients)

    ;Create a curved 2D profile
    profile_2d = dblarr(nx,ny)+!values.d_nan
    left_edge_ypos = poly(xarr,ord_left_coeffs[*,good_orders[i]])
    int_left_edge_ypos = floor(left_edge_ypos)
    frac_left_edge_ypos = left_edge_ypos-double(int_left_edge_ypos)
    for l=0L, nx-1L do begin
      ;Skip this column if the order is partially out of frame
      if left_edge_ypos[l] lt 0d0 then continue
      if left_edge_ypos[l]+double(height-1L) gt double(ny-1L) then continue
      profile_2d[l,int_left_edge_ypos[l]:(int_left_edge_ypos[l]+height-1L)] = interpol2(trace_profile,subyarr,subyarr+(trace_max_pos-(y_position_refined[l]-double(int_left_edge_ypos[l]))),badvalue=0d0,/repairnans)
    endfor

    ;Measure the fraction of pixels with finite data in each row
    frac_finite_pix_per_row = total(finite(smoothed_data_order),1L)
    frac_finite_pix_per_row /= max(frac_finite_pix_per_row,/nan)

    ;Determine what rows should be used to measure the sky:
    ; At least 80% of pixels must be finite compared to the best row
    ; and then the "Nrows_sky" pixels with the least trace profile flux are chosen
    g_enough_finite = where(frac_finite_pix_per_row ge .8, ng_enough_finite)
    sort_indices = sort(trace_profile[g_enough_finite])
    sky_pos = g_enough_finite[sort_indices[0L:nrows_sky-1L]]

    ;Sort the sky position indices
    sky_pos = sky_pos[sort(sky_pos)]

    ;Measure the sky
    ; It is normal that a bit of the science trace bleeds into the sky
    ;  when the target is bright; the slit is not long enough to avoid this
    ;  at least the relative sky flux vs wavelength should be relatively accurate
    ; A consequence of this is that the measured S/N will be slightly under-estimated
    sky = median(smoothed_data_order[*,sky_pos],dim=2L)
    sky_2d = sky#make_array(ny,value=1d0,/double)

    ;Do a raw extraction
    raw_sp = total((im_order-sky_2d)*profile_2d,2,/nan)/total(profile_2d,2,/nan)

    ;Filter out bad pixels
    raw_sp_bpix = correct_bad_pixels_posterior(raw_sp)

    ;Replace NaNs in the raw extraction with a median smooth
    for rep=0L, 1L do begin
      bad = where(~finite(raw_sp_bpix), nbad)
      if nbad ne 0L then raw_sp_bpix[bad] = (median(raw_sp_bpix,3))[bad]
    endfor

    ;Recreate a 2D spectrum of the ata that contains less bad pixels
    raw_sp_bpix_2d = (raw_sp_bpix#make_array(ny,value=1d0,/double))

    ;Create a very smoothed version of the data
    smooth_spectrum = horizontal_median(raw_sp_bpix_2d,75)

    ;Alpha factor is the factor to express the maximum flux on the trace instead of the total flux within the trace,
    ; assuming that the profile is a good representation of the data
    alpha_factor = total(trace_profile,/nan)/total(trace_profile^2,/nan)

    ;Do some bad pixel recognizing by smoothing the data in both directions and looking at outliers
    deviations = (abs(im_order-horizontal_median(im_order,3))<abs(im_order-vertical_median(im_order,3)))/(smooth_spectrum*alpha_factor+abs(sky_2d))

    ;Identify bad pixels
    ;The .08 threshold should be a free parameter. We also have to make sure that it works well in general situations
    bad = where(abs(deviations) ge bad_pixels_threshold, nbad)

    ;Create a bad pixel mask
    badpix_mask = dblarr(nx,ny)+1d0
    if nbad ne 0L then badpix_mask[bad] = 0d0

    ;Re-do a crude extraction with the bad pixel mask
    raw_sp2 = total((im_order-sky_2d)*profile_2d*badpix_mask,2,/nan)/total(profile_2d*badpix_mask,2,/nan)

    ;Mask pixels with no data
    bad = where(total(finite((im_order-sky_2d)*profile_2d*badpix_mask),2) eq 0, nbad)
    if nbad ne 0L then raw_sp2[bad] = !values.d_nan

    ;Compute the effective read noise in electrons
    eff_read_noise = read_noise + dark_current*integration_times[g_science[sci]]

    ;Do the optimal extraction
    opt_spectrum = general_optimal_extract_pmassey(im_order, sky_2d, profile_2d*badpix_mask, eff_read_noise, gain, ESP=err_opt_spectrum)

    ;Some debugging displays:
    ;flat_spectrum = general_optimal_extract_pmassey(flat_field_corrected, sky_2d*0., profile_2d*badpix_mask, eff_read_noise, gain);, NITER=optimal_bpix_iterations)
    ;corr_spectrum=opt_spectrum/flat_spectrum & plot,opt_spectrum,yrange=[0,5d4] & oplot,corr_spectrum/min(corr_spectrum/opt_spectrum,/nan),col=rvb_hex(100,255,100) & oplot,opt_spectrum & oplot,flat_spectrum/median(flat_spectrum/opt_spectrum),col=255

    ;Do the complete 2D fit of the trace if needed
    ; This part may be more prone to crashing, but it will be better at removing bad pixels,
    ; and will account for wavelength-dependent seeing. It will also produce more refinements
    ; to the Y trace position polynomial.
    ;This option should absolutely not be used if the trace profile is very different from a gaussian.
    if keyword_set(do_trace_2d_fit) then begin
      ;Estimate the parameters that will be used to create a 2D model of the data
      estimated_parameters = double([refined_coefficients, gauss_parameters[2L], dblarr(ndegree_sigma_fit-1L)])

      ;Create parameters limits and start values
      n_parameters = n_elements(estimated_parameters)
      parinfo = replicate({VALUE:!values.d_nan,FIXED:0B,LIMITED:[0B,0B],LIMITS:[0.,0.]},n_parameters)
      parinfo.VALUE = estimated_parameters

      ;Fix the width at the origin (x=0) to a positive value
      parinfo[ndegree_poly_fit].LIMITED = [1,0]
      parinfo[ndegree_poly_fit].LIMITS[0] = 0d0

      ;Create maps of x and y pixel values
      xmap = xarr#make_array(ny,value=1d0,/double)
      ymap = make_array(nx,value=1d0,/double)#yarr

      ;Create a data-sky array and make sure there are no NaNs in it
      data_minus_sky = im_order-sky_2d
      bad = where(~finite(data_minus_sky), nbad)
      if nbad ne 0L then begin
        data_minus_sky[bad] = 0d0
        badpix_mask[bad] = 0d0
      endif

      ;Give more weight in the fitting to the center of the order where there are less light scattering effects
      datafit_weights = badpix_mask*flat_field_corrected*sqrt(profile_2d)
      bad = where(~finite(datafit_weights), nbad)
      if nbad ne 0L then datafit_weights[bad] = 0d0

      ;Data to be passed to the fitting function (with the _extra=extra IDL mechanism, mpfit2dfun takes care of this)
      opt_spectrum_2d = opt_spectrum#make_array(ny,value=1d0,/double)
      functargs = {spectrum:opt_spectrum_2d,npoly_ypos:ndegree_poly_fit,data:data_minus_sky,weights:datafit_weights}

      ;Do the 2D fit on the trace, our forward model function is called fit_2d_curved_trace.pro
      print, '  performing 2D fit of the trace...'
      best_fit_parameters = mpfit2dfun('fit_2d_curved_trace',xmap,ymap,data_minus_sky,0., WEIGHTS=datafit_weights, /QUIET, PARINFO=parinfo, FUNCTARGS=functargs, STATUS=status, ERRMSG=errmsg)
      if status le 0 then $
        message, 'The 2D trace fit has failed !'

      ;Use the new fit to flag additional bad pixels
      new_deviations = abs((data_minus_sky-fit_2d_curved_trace(xmap, ymap, best_fit_parameters, _EXTRA=functargs, /nobpixpatch))*badpix_mask/weighted_median(opt_spectrum,medval=quartile_fraction_for_norm))
      bad = where(new_deviations ge second_bad_pixels_threshold, nbad)
      if nbad ne 0L then $
        badpix_mask[bad] = 0.0

      ;Recreate a 2D profile with the 2D trace model that allows for a linear
      ;  dependence on seeing vs spatial pixel position
      if keyword_set(model_refinement_curvature) then begin
        refined_again_coefficients = best_fit_parameters[0:ndegree_poly_fit-1L]
        sigma_coefficients = best_fit_parameters[ndegree_poly_fit:*]
        sigma = poly(xmap, sigma_coefficients)
        gauss_y = poly(xmap, refined_again_coefficients)
        profile_2d = exp(-(ymap-gauss_y)^2/(2d0*sigma^2))
      endif

      ;Store the coefficients in the full-orders array of coefficients
      full_model_polynomials[*,i] = best_fit_parameters

      ;Re-do the optimal extraction with the refined parameters
      opt_spectrum = general_optimal_extract_pmassey(im_order, sky_2d, profile_2d*badpix_mask, eff_read_noise, gain, ESP=err_opt_spectrum)

      ;Apply a final bad pixel masking step
      ;opt_spectrum = correct_bad_pixels_posterior(opt_spectrum)

      ;Refine the 2D model even more for output. This can be useful for debugging.
      opt_spectrum_2d = opt_spectrum#make_array(ny,value=1d0,/double)
      functargs.spectrum = opt_spectrum_2d
      updated_2d_model = fit_2d_curved_trace(xmap, ymap, best_fit_parameters, _EXTRA=functargs, /nobpixpatch)
      wbadpix = where(badpix_mask eq 0, nwbadpix)
      if nwbadpix ne 0L then updated_2d_model[wbadpix] = !values.d_nan

      ;Save the refined 2D model in IDL save format
      models_file = models_dir+'models_order'+strtrim(ord_ids[good_orders[i]],2L)+'_'+object_name+'_'+file_basename(data_file,file_ext(data_file))+'.sav'
      save, best_fit_parameters, functargs, updated_2d_model, badpix_mask, file=models_file, /compress

      ;Put the model at the correct order position in the full model array
      full_2d_model[g_order] = updated_2d_model[g_order]

    endif

    ;Normalize the raw and optimal spectra to relative fluxes
    medval = weighted_median(opt_spectrum,medval=quartile_fraction_for_norm)
    rel_norm = median((raw_sp2/opt_spectrum))
    opt_spectrum_norm = opt_spectrum/medval
    raw_spectrum_norm = raw_sp2/(medval*rel_norm)

    ;Save the intermediate reduction results for the full-orders diagnostic figures
    extraction_data_file = extraction_idl_dir+'extraction_idl_data'+strtrim(ord_ids[good_orders[i]],2L)+'_'+object_name+'_'+file_basename(data_file,file_ext(data_file))+'.sav'
    save, opt_spectrum_norm, raw_spectrum_norm, trace_profile, ytrace_precise, gauss_fit_precise, /compress, file=extraction_data_file

    ;Store these results in arrays
    full_orders_spectra[*,i] = opt_spectrum_norm
    full_orders_raw_spectra[*,i] = raw_spectrum_norm
    *full_orders_trace_profile[i] = trace_profile
    *full_orders_trace_xfit[i] = ytrace_precise
    *full_orders_trace_yfit[i] = gauss_fit_precise

    ;Create a single-order figure with the data
    if generate_individual_orders_spectra_figures then begin
      ;Plot the optimal extraction
      gfinite = where(finite(opt_spectrum),ngfinite)
      a = plot(buffer=1,opt_spectrum_norm,xtitle='Pixels (Order '+strtrim(ord_ids[good_orders[i]],2)+')',ytitle='Relative flux',thick=2,color='black',margin=[.07,.12,.02,.02],xthick=2,ythick=2,yticklen=.015,xticklen=.015,dimensions=[1000,512],xrange=[gfinite[0L]-3L,gfinite[-1L]+3L],name='Optimal')
      ;Plot the raw extraction with a bad pixel mask
      b = plot(/overplot,raw_spectrum_norm,yrange=a.yrange,color=whiten_color('red',40),thick=2,name='Regular')
      ;Plot the raw extraction without a bad pixel mask (hidden as it causes too much clutter)
      ;c = plot(/overplot,raw_sp/(medval*rel_norm),yrange=a.yrange,color=whiten_color('red',40),thick=2)
      ;Arrange the display order
      b.order,/send_to_back
      ;c.order,/send_to_back
      ;Plot a legend
      l = legend(target=[a,b],position=[.97,.95],horizontal_align=1,vertical_align=1,font_size=10,transparency=50,shadow=0,thick=2)
      ;Increase yrange by 5% each side
      yrange = [min(opt_spectrum_norm,/nan),max(opt_spectrum_norm,/nan)]
      yrange += [-1,1]*(max(yrange)-min(yrange))*.05
      ;Limit max yrange to 1.2
      yrange = [yrange[0L],yrange[1L]<1.2]
      a.yrange = yrange
      ;Save the figure
      a.save, previews_dir+file_basename(data_file,file_ext(data_file,/FITS))+'_OBJ_'+object_name+'_ORD_'+strtrim(ord_ids[good_orders[i]],2)+'_preview.png'
      ;Close the figure
      a.close
    endif

    ;Normalize the data for output
    out_data = opt_spectrum_norm
    out_errors = err_opt_spectrum/medval

    ;Idenfity bad pixels to output bad pixel mask
    data_mask = lonarr(n_elements(opt_spectrum_norm))+1
    bad = where(~finite(out_data), nbad)
    if nbad ne 0L then begin
      out_data[bad] = 0d0
      out_errors[bad] = 0d0
      data_mask[bad] = 0
    endif

    ;Add useful information in the header for output
    sxaddpar, header, 'FIT_2D', do_trace_2d_fit, 'Whether a full 2D forward model was fitted'
    sxaddpar, header, 'DARK_SUB', do_dark_subtraction, 'Whether darks were subtracted'
    sxaddpar, header, 'BLZ_CORR', correct_blaze_function_in_flatfield, 'Whether the blaze fcn was corrected in flats'
    sxaddpar, header, 'FRG_CORR', correct_fringing_in_flatfield, 'Whether the fringing was corrected in flats'
    sxaddpar, header, 'MOD_FRG', model_fringing, 'Whether the fringing modelled to better remove it'
    sxaddpar, header, 'HNPX_AVD', edge_npix_avoid, 'Number of additional avoided edge pixels'
    sxaddpar, header, 'HNPX_MSK', edge_npix_mask, 'Number of masked horizontal edge pixels'
    sxaddpar, header, 'BPIX_TD', bad_pixels_threshold, 'The fractional threshold for bad pixels'
    if do_trace_2d_fit eq 1 then $
      sxaddpar, header, 'BPIX_TD2', second_bad_pixels_threshold, 'The 2nd fractional threshold for bad pixels'
    sxaddpar, header, 'NORM_QRT', quartile_fraction_for_norm, 'Quartile used for normalization'
    sxaddpar, header, 'SIZELAG', sizelag, 'Pixel size of lag used to refine y trace pos'
    sxaddpar, header, 'LAGNFIT', window_lag_poly_fit, 'Number of pixels fitted in cross-correlation'
    sxaddpar, header, 'POLYDROP', refine_poly_drop_pixels, 'Number of edge pix to drop in ytrace refinement'
    sxaddpar, header, 'TRCPOLY', ndegree_poly_fit, 'Degree of polynomial to fit to the trace'
    sxaddpar, header, 'PROFMASK', mask_trace_edges, 'Number of edge pixels to mask on trace profile'
    sxaddpar, header, 'SKYNROWS', nrows_sky, 'Number of rows used to measure sky'
    sxaddpar, header, 'DETGAIN', gain, 'Detector gain in e/ADU'
    sxaddpar, header, 'RDNOISE', read_noise, 'read noise in electrons'
    sxaddpar, header, 'DRKCURR', dark_current, 'dark current in electrons per second'
    sxaddpar, header, 'DATE_RED', curcompdate(), 'Date of reduction YYMMDD'
    sxaddpar, header, 'CODEVER', code_version, 'Version of data reduction code'
    sxaddpar, header, 'RMDTPTRN', remove_detector_patterns_from_data, 'Remove detector patterns with block median'
    sxaddpar, header, 'NOFLATS', override_flats, '/override_flats flag to skip use of flat fields'
    sxaddpar, header, 'REFCURVE', model_refinement_curvature, 'Use 2D model to refine trace curvature'

    ;Write the data to a file
    openw, lun, output_file, /get_lun
    for kk=0L, n_elements(header)-1L do $
      printf, lun, ';'+header[kk]
    printf, lun, ';RELATIVE FLUX AS A FUNCTION OF PIXEL, TAB DELIMITER, AND VALUE OF BAD PIXEL MASK (1=good, 0=bad):'
    printf, lun, transpose(strtrim(out_data,2L)+string(9b)+strtrim(out_errors,2L)+string(9b)+strtrim(data_mask,2L))
    close, lun
    free_lun, lun
  endfor

  ;Once the loop over orders is finished, save the full 2D models file and the data - model residuals as well (if needed)
  if keyword_set(do_trace_2d_fit) then begin
    if generate_2d_model_fits_images then begin
      full_model_file = models_dir+'full_model_'+object_name+'_'+file_basename(data_file,file_ext(data_file))+'.fits.gz'
      writefits, full_model_file, full_2d_model, /silent, /compress
      full_model_polynomials_file = models_dir+'full_polynomials_'+object_name+'_'+file_basename(data_file,file_ext(data_file))+'.txt'
      openw, lun2, full_model_polynomials_file, /get_lun
      printf, lun2, '; '+strtrim(ndegree_poly_fit,2)+' coefficients for Y trace pos + '+strtrim(ndegree_sigma_fit,2)+' coefficients for seeing versus wavelength (columns)'
      printf, lun2, '; For all '+strtrim(n_orders,2)+' orders (rows)'
      printf, lun2, '; Trace Y position (pixel) on full '+strtrim(nx,2)+'x'+strtrim(ny,2)+' detector versus X (pixel), coefficients a_0 through a_'+strtrim(ndegree_poly_fit-1,2)+' (Y = SUM_i a_i*x^i)'
      printf, lun2, '; Trace width (Gaussian sigma in pixels) versus X (pixel), coefficients b_0 through b_'+strtrim(ndegree_sigma_fit-1,2)+' (SIGMA = SUM_i b_i*x^i)'
      printf, lun2, full_model_polynomials
      close, lun2
      free_lun, lun2
    endif
    if generate_residuals_fits_images then begin
      full_residuals_file = residuals_dir+'full_residuals'+'_'+object_name+'_'+file_basename(data_file,file_ext(data_file))+'.fits.gz'
      writefits, full_residuals_file, data_im-full_2d_model, /silent, /compress
    endif
  endif

  ;Create the full-orders trace profile figure
  if generate_full_orders_trace_profile_figures then begin
    pa = !NULL
    yshift = 0d0
    ydel = 0.3d0
    for i=0L, n_orders-1L do begin
      if *full_orders_trace_profile[i] eq !NULL then continue
      if *full_orders_trace_yfit[i] eq !NULL then continue
      if *full_orders_trace_xfit[i] eq !NULL then continue
      if pa eq !NULL then begin
        pa = plot(buffer=1,dindgen(n_elements(*full_orders_trace_profile[i])), *full_orders_trace_profile[i], thick=2, name='Trace Profile', yticklen=.015, xticklen=.015, xthick=2, ythick=2, margin=[.12,.12,.02,.02],xtitle='Pixel',ytitle='Relative Flux + Shift')
        pb = plot(*full_orders_trace_xfit[i],*full_orders_trace_yfit[i], thick=3, name='Gaussian Fit', /overplot, color='red')
        pb.order, /send_to_back
        yshift += ydel
        continue
      endif
      pbi = plot(*full_orders_trace_xfit[i],*full_orders_trace_yfit[i]+yshift, thick=3, /overplot, color='red')
      pai = plot(/overplot,dindgen(n_elements(*full_orders_trace_profile[i])), *full_orders_trace_profile[i]+yshift, thick=2)
      yshift += ydel
    endfor
    pa.yrange = [0., double(n_orders-1L)*ydel+1.2]
    pl = legend(target=[pa,pb],position=[.95,.95],horizontal_align=1,vertical_align=1,font_size=10,transparency=50,shadow=0,thick=2)
    trace_profiles_file = trace_profiles_dir+'full_orders_trace_profile_'+object_name+'_'+file_basename(data_file,file_ext(data_file))+'.png'
    pa.save, trace_profiles_file
    pa.close
  endif

  ;Make a full-orders plot of the extraction
  if generate_full_orders_spectra_figures then begin
    nx_win = 3L
    ny_win = ceil(double(n_orders)/double(nx_win))
    window_id = window(dimensions=[1000L*nx_win,512L*n_orders/nx_win]/2L,buffer=1)
    for i=0L, n_orders-1L do begin
      layout = [nx_win,ny_win,i+1]
      gfinite = where(finite(full_orders_spectra[*,i]),ngfinite)
      pai = plot(full_orders_spectra[*,i],xtitle='Order '+strtrim(i,2),thick=2,xrange=[gfinite[0L]-3L,gfinite[-1L]+3L],$
        color='black',margin=[.07,.12,.02,.02],xthick=2,ythick=2,yticklen=.015,xticklen=.015,layout=layout,/current)
      pbi = plot(/overplot,full_orders_raw_spectra[*,i],color='red',thick=2)
      pbi.order,/send_to_back
      ;Increase yrange by 5% each side
      yrange = [min(full_orders_spectra[*,i],/nan),max(full_orders_spectra[*,i],/nan)]
      yrange += [-1,1]*(max(yrange)-min(yrange))*.05
      ;Limit max yrange to 1.2
      yrange = [yrange[0L],yrange[1L]<1.2]
      pai.yrange = yrange
    endfor
  endif

  full_orders_spectra_file = previews_dir+'full_orders_spectra_preview_'+file_basename(data_file,file_ext(data_file,/FITS))+'_OBJ_'+object_name+'.png'
  window_id.save, full_orders_spectra_file
  pai.close
  
End