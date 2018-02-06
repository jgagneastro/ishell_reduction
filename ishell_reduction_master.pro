;DONE: ishell_reduction_master
;"simple block filter" ishell_reduction_master,model_fringing=0,remove_detector_patterns_from_data=1
;"noflats": ishell_reduction_master,/override_flats,model_fringing=0,correct_fringing_in_flatfield=0
;>>"nofringing_corr_": ishell_reduction_master,model_fringing=0,correct_fringing_in_flatfield=0
Pro ishell_reduction_master, data_path, output_dir_root, DEBUG_TRACE_ORDER=debug_trace_order, DO_DARK_SUBTRACTION=do_dark_subtraction, $
  CORRECT_FRINGING_IN_FLATFIELD=correct_fringing_in_flatfield, MODEL_FRINGING=model_fringing, $
  REMOVE_DETECTOR_PATTERNS_FROM_DATA=remove_detector_patterns_from_data, OVERRIDE_FLATS=override_flats
  ;Code version history
  ; - The code started at 0.9, between 0.9 and 1.0 only minor bugs were fixed and more diagnostic outputs were added
  ; Version 1.0: First stable version (J. Gagne)
  ; Version 1.1: Added outputs of polynomial coefficients for Ytrace pos and seeing versus X pixel in ASCII format (J. Gagne), March 3, 2017
  ; Version 1.2: Added N_orders keyword for compatibility with K band, January 26, 2018
  ; Version 1.3: More fixes for mixed bands in a single night, added an option to remove column-dependent detector patterns from the data, January 31, 2018
  ; Version 1.4: Added option to model fringing to separate it from the column-dependent detector patterns. February 4, 2018
  
  ;Code version for headers
  code_version = 1.4
  
  ;List of subroutines
  forward_function readfits,ishell_trace_orders,ishell_flat_fringing,interpol2,weighted_median,horizontal_median,$
    vertical_median,sxpar,whiten_color, sxpar_mul, rtrim, readcol, printuarr, match_unsorted, $
    strtrim_multiple, nan_str, mpfit2dfun, fit_2d_curved_trace, general_optimal_extract_pmassey, $
    correct_bad_pixels_posterior, sxaddpar, strkill
  
  ;/// Adjustable parameters ///
  
  ;The path of the night which will be reduced
  if ~keyword_set(data_path) then $
    ;data_path = '/Users/gagne/Documents/Data_Repository/RAW/IRTF/iShell/20171023UT/'
    ;data_path = '/Users/gagne/Documents/Data_Repository/RAW/IRTF/iShell/20161016UT_Vega_night1/'
    ;data_path = '/Users/gagne/Documents/Data_Repository/RAW/IRTF/iShell/20161107UT_Vega_night2/'
    data_path = '/Users/gagne/Documents/Data_Repository/RAW/IRTF/iShell/20171023UT/'
  
  ;The path where the data reduction products will be stored
  if ~keyword_set(output_dir_root) then $
    output_dir_root = '/Users/gagne/Documents/Data_Repository/RAW/IRTF/iShell/redux/'
  
  ;Make sure that paths end with a path separator
  if strmid(data_path,0,1) ne path_sep() then $
    data_path += path_sep()
  if strmid(output_dir_root,0,1) ne path_sep() then $
    output_dir_root += path_sep()
  
  ;Whether or not to do debugging for trace order detection
  if debug_trace_order eq !NULL then $
    debug_trace_orders = 0
  
  ;Whether or not darks should be subtracted
  if do_dark_subtraction eq !NULL then $
    do_dark_subtraction = 0
  
  ;Whether fringing should be modelled to separate it better from column-dependent detector patterns
  if model_fringing eq !NULL then $
    model_fringing = 1
  
  ;Whether or not column-dependent detector patterns should be preserved in the flat fields and therefore removed from the data
  ;This option is obsolete because it does not appropriately eliminate fringing from the flat fields.
  ;Instead use model_fringing = 1
  if remove_detector_patterns_from_data eq !NULL then $
    remove_detector_patterns_from_data = 0
  
  ;Whether fringing should be corrected in the flat field
  if correct_fringing_in_flatfield eq !NULL then $
    correct_fringing_in_flatfield = 1
  
  ;Avoid N pixels on each side of the detector when normalizing
  edge_npix_avoid = 200L
  
  ;Mask N pixels on each side of the detector horizontally (this is where the Blaze function is at its lowest)
  edge_npix_mask = 200L 
  
  ;Perform a complete 2D fitting of the trace with a gaussian profile whose width
  ; vary linearly with wavelength (should produce better results, but might crash more often) 
  do_trace_2d_fit = 1
  
  ;Which figures should be output
  generate_full_orders_spectra_figures = 1
  generate_individual_orders_spectra_figures = 1
  generate_full_orders_trace_profile_figures = 1
  generate_individual_orders_trace_profile_figures = 1
  generate_residuals_fits_images = 1
  generate_2d_model_fits_images = 1
  
  ;Header key for object name 
  object_header_key = 'OBJECT'
  
  ;Header key for integration time 
  exptime_header_key = 'ITIME'
  
  ;Header key for comments
  comment_header_key = 'COMMENT'
  
  ;Header key for gas cell position
  gascell_header_key = 'GASCELL'
  
  ;Header key for filter
  filter_header_key = 'XDTILT'
  
  ;Header key for slit
  slit_header_key = 'SLIT'
  
  ;Header key for TCS object name
  tcs_obj_header_key = 'TCS_OBJ'
  
  ;Threshold for bad pixel detection from deviations in terms of % of the flux
  bad_pixels_threshold = .08
  
  ;Threshold for bad pixel detection (in terms of % of the flux) on second filter
  ;This is only relevant if do_trace_2d_fit = 1
  second_bad_pixels_threshold = .5
  
  ;Should the blaze function and lamp spectra be corrected in the flat fields ?
  correct_blaze_function_in_flatfield = 0
  
  ;The fractional quartile value of the pixel to use for normalization 
  quartile_fraction_for_norm = .99
  
  ;The total amount of pixels to be y-shifted when the trace position is refined  
  sizelag = 20L
  
  ;Size of the window (in pixels) around the max of the cross-correlation function that
  ; should be used to fit a 2nd order polynomial
  window_lag_poly_fit = 5L
  
  ;The number of data points to be dropped on both ends of the order for polynomial refinement
  refine_poly_drop_pixels = 5L
  
  ;The number of degrees for the polynomial fit to the trace position
  ndegree_poly_fit = 3L
  
  ;The number of degrees for the polynomial fit to the seeing versus wavelength (only used if do_trace_2d_fit = 1)
  ndegree_sigma_fit = 2L
  
  ;Number of pixels to mask on each edge of the trace profile
  mask_trace_edges = 3L
  
  ;Number of rows used to measure the sky
  nrows_sky = 8L
  
  ;Gain of the detector in electrons / ADU
  ;data from http://irtfweb.ifa.hawaii.edu/~ishell/iSHELL_observing_manual.pdf
  gain = 1.8d0
  
  ;Read noise in electrons, assuming NREAD = 32 with the 30-stripes slow readout mode
  ;data from http://irtfweb.ifa.hawaii.edu/~ishell/iSHELL_observing_manual.pdf
  read_noise = 8d0
  
  ;Median measured dark current in electrons per second
  ;http://irtfweb.ifa.hawaii.edu/~ishell/iSHELL_observing_manual.pdf
  dark_current = 0.05
  
  ;/// End of: adjustable parameters ///
  
  ;Add the night ID after the output directory
  date_id = file_basename(data_path)
  output_dir = output_dir_root+date_id+path_sep()
  
  ;Recursively create output directory if it does not exist
  if ~file_test(output_dir) then begin
    dirs = strsplit(output_dir,path_sep(),/extract)
    ndir = n_elements(dirs)
    current_path = path_sep()
    for j=0L, ndir-1L do begin
      if file_test(current_path+dirs[j]+path_sep()) then begin
        current_path += dirs[j]+path_sep()
        continue
      endif else begin
        file_mkdir, current_path+dirs[j]+path_sep()
        current_path += dirs[j]+path_sep()
      endelse
    endfor
  endif
  
  ;Identify all fits data files (they must not be in subdirectories)
  fits_data = file_search(data_path+'*.fits')
  ndata = n_elements(fits_data)
  
  ;Read file header information
  object_names = sxpar_mul(fits_data,object_header_key)
  integration_times = sxpar_mul(fits_data,exptime_header_key)
  data_comments = sxpar_mul(fits_data,comment_header_key)
  tcs_obj = sxpar_mul(fits_data,tcs_obj_header_key)
  gascell_position = sxpar_mul(fits_data,gascell_header_key)
  data_filters = sxpar_mul(fits_data,filter_header_key)
  data_slits = sxpar_mul(fits_data,slit_header_key)
  
  ;Clean up header information
  strkill, object_names, ['=','''',',','|']
  strkill, data_comments, ['=','''',',']
  strkill, tcs_obj, ['=','''',',','|']
  strkill, data_slits, '|'
  strkill, data_filters, '|'
  
  ;Remove trailing spaces
  object_names = strtrim(object_names,2)
  data_comments = strtrim(data_comments,2)
  tcs_obj = strtrim(tcs_obj,2)
  data_slits = strtrim(data_slits,2)
  data_filters = strtrim(data_filters,2)
  
  ;Check whether the log file exists
  logfile = output_dir+'logfile_'+file_basename(data_path)+'.txt'
  
  if file_test(logfile) then begin
    print, ' Found an existing log file, using it...'
  endif else begin
    
    ;By default, darks or ThAr data won't be reduced
    do_reduce = replicate('Yes',ndata)
    bad = where(strpos(strlowcase(object_names),'dark') ne -1 or strpos(strlowcase(object_names),'thar') ne -1, nbad)
    if nbad ne 0L then do_reduce[bad] = 'No'
    
    ;Create an crude log file and put some header information in it 
    printuarr, logfile, file_basename(fits_data), object_names, rtrim(integration_times,3), data_comments, tcs_obj, $
      gascell_position, data_filters, data_slits, do_reduce, title=[';File','Object','ExpTime','Comment','TCS Object','GasCell','Filter','Slit','Reduce'], $
      /justify, symbol='|'
    
    print, ' A log file was just created at :'
    print, '   '+logfile
    print, ' Please review header information, make sure that object names are consistent and correct, (...)'
    print, '   chose No instead of Yes to skip reducing any file, and make sure that all flats (...)'
    print, '   are named QTH. Only flats with "Recude = Yes" will be used. (...)'
    print, '   Flats can have an object name in their "TCS Object" header position, in which case they will (...)'
    print, '   be used with the appropriate object. When ready, launch this code again.'
    return
  endelse

  ;Read the log file back
  readcol, logfile, filenames_log, object_names, integration_times, data_comments, tcs_obj, $
    gascell_position, data_filters, data_slits, do_reduce, comment=';', delimiter='|', $
    format='A,A,F,A,A,A,A,A,A', /silent
  
  ;Remove trailing spaces
  strtrim_multiple, filenames_log, object_names, data_comments, gascell_position, data_filters, data_slits, do_reduce
  
  ;Make sure that the data files did not change
  if ~array_equal(filenames_log,file_basename(fits_data)) then begin
    file_move, logfile, logfile+'.old'
    print, ' Mismatches were found in file names ! The log file has been moved away. Please run the code again.'
    return
  endif
  
  ;Identify all darks
  if do_dark_subtraction eq 1 then begin
    g_darks = where(strpos(strlowcase(object_names),'dark') ne -1 and strlowcase(data_slits) eq 'mirror' and strlowcase(do_reduce) eq 'yes', ng_darks)
    if ng_darks eq 0 then begin
      message, ' No darks were found in the data ! Skipping the dark subtraction step.', /continue
      do_dark_subtraction = 0
    endif
  endif
  
  if do_dark_subtraction eq 1 then begin
  
    ;Identify all types of dark IDs (i.e., their exposure times)
    darks_ids = rtrim(integration_times[g_darks],3)
    
    ;Create a set of unique dark IDs
    darks_ids_uniq = darks_ids
    darks_ids_uniq = darks_ids_uniq[sort(darks_ids_uniq)]
    darks_ids_uniq = darks_ids_uniq[uniq(darks_ids_uniq)]
    ndarks_uniq = n_elements(darks_ids_uniq)
    
    ;If combined darks already exist, restore them
    darks_dir = output_dir+'darks'+path_sep()
    if ~file_test(darks_dir) then file_mkdir, darks_dir
    comb_darks_file = darks_dir+'combined_darks_'+date_id+'.sav'
    if file_test(comb_darks_file) then begin
      restore, comb_darks_file;, darks_uniq_cube, nx, ny
    endif else begin
  
      print, ' Median-combining dark exposures...'
  
      ;Create combined dark fields for each data ID
      for f=0L, ndarks_uniq-1L do begin
        ;Identify all darks with the appropriate ID
        g_darks_f = where(darks_ids eq darks_ids_uniq[f], ng_darks_f)
        darks_f = fits_data[g_darks[g_darks_f]]
        ndarks_f = n_elements(darks_f)
        ;Reset the darks cube
        darks_f_cube = !NULL
        ;Read the required darks
        for subf=0L, ndarks_f-1L do begin
          dark_im = double(readfits(darks_f[subf],/silent))
  
          ;Determine array size and create a cube to store darks
          if ~keyword_set(nx) then begin
            nx = (size(dark_im))[1]
            ny = (size(dark_im))[2]
            darks_uniq_cube = dblarr(nx,ny,ndarks_uniq)+!values.d_nan
          endif
          if ~keyword_set(darks_f_cube) then $
            darks_f_cube = dblarr(nx,ny,ndarks_f)+!values.d_nan
          ;Store dark in cube
          darks_f_cube[*,*,subf] = dark_im
        endfor
  
        ;Median-combine the dark and store it in cube
        darks_uniq_cube[*,*,f] = median(darks_f_cube,dim=3)
        
        ;Save the dark as a fits file
        writefits, darks_dir+'dark_medcomb_'+date_id+'_TEXP'+darks_IDs_uniq[f]+'s.fits', darks_uniq_cube[*,*,f]
      endfor
      save, darks_uniq_cube, nx, ny, file=comb_darks_file, /compress
    endelse
  endif
  
  ;If combined flats already exist, restore them
  nx = !NULL
  flats_dir = output_dir+'flats'+path_sep()
  if ~file_test(flats_dir) then file_mkdir, flats_dir
  comb_flats_file = flats_dir+'combined_flats_'+date_id+'.sav'
  if file_test(comb_flats_file) then begin
    restore, comb_flats_file;, flats_uniq_cube, nx, ny, flat_ids_uniq, nflat_uniq
  endif else begin
    
    print, ' Median-combining flat field exposures...'
    
    ;Identify all flat fields
    g_flats = where(object_names eq 'QTH' and strlowcase(gascell_position) eq 'out' and strlowcase(do_reduce) eq 'yes', ng_flats)
    if ng_flats eq 0 then $
      message, ' No flats were found in the data !'

    ;Determine all different flat field identifications from their comments
    flat_ids = tcs_obj[g_flats]+'|'+data_filters[g_flats]+'|'+data_slits[g_flats]
    bad = where(flat_ids eq '', nbad)
    if nbad ne 0L then flat_ids[bad] = 'BLANK'

    ;Create a set of unique flat IDs
    flat_ids_uniq = flat_ids
    flat_ids_uniq = flat_ids_uniq[sort(flat_ids_uniq)]
    flat_ids_uniq = flat_ids_uniq[uniq(flat_ids_uniq)]
    nflat_uniq = n_elements(flat_ids_uniq)
    
    ;Create combined flat fields for each data ID
    for f=0L, nflat_uniq-1L do begin
      ;Identify all flats with the appropriate ID
      g_flat_f = where(flat_ids eq flat_ids_uniq[f], ng_flat_f)
      flats_f = fits_data[g_flats[g_flat_f]]
      nflats_f = n_elements(flats_f)
      ;Reset the flat cube
      flat_f_cube = !NULL
      ;Read the required flats
      for subf=0L, nflats_f-1L do begin
        flat_im = double(readfits(flats_f[subf],/silent))
        
        ;Determine array size and create a cube to store flats
        if ~keyword_set(nx) then begin
          nx = (size(flat_im))[1]
          ny = (size(flat_im))[2]
          flats_uniq_cube = dblarr(nx,ny,nflat_uniq)+!values.d_nan
        endif
        if ~keyword_set(flat_f_cube) then $
          flat_f_cube = dblarr(nx,ny,nflats_f)+!values.d_nan
        ;Normalize flat
        flat_im /= weighted_median(flat_im,medval=quartile_fraction_for_norm)
        ;Store flat in cube
        flat_f_cube[*,*,subf] = flat_im
      endfor
      
      ;Median-combine the flat field and store it in cube
      flats_uniq_cube[*,*,f] = median(flat_f_cube,dim=3)
    endfor
    save, flats_uniq_cube, nx, ny, flat_ids_uniq, nflat_uniq, file=comb_flats_file, /compress
  endelse
  
  ;Build all orders masks
  max_n_orders = 32L;This should always be updated when adding stuff to the case statement below
  orders_mask_cube = dblarr(nx,ny,nflat_uniq)+!values.d_nan
  for f=0L, nflat_uniq-1L do begin
    ;Verify whether this was already done and saved to disk
    order_mask_dir = output_dir+'order_masks'+path_sep()
    if ~file_test(order_mask_dir) then file_mkdir, order_mask_dir
    order_mask_file = order_mask_dir+'order_mask_'+date_id+'_ID_'+flat_ids_uniq[f]+'.sav'
    if file_test(order_mask_file) then begin
      restore, order_mask_file;, orders_mask_f, orders_structure_f, n_orders, min_order_spacing
    endif else begin
      
      ;Determine the number of orders and the minimum order spacing expected for each filter
      current_filter = (strsplit(flat_ids_uniq[f],'|',/extract))[1]
      case strlowcase(current_filter) of
        ;KS-band data has 29 orders and 15 pixels minimum between each order
        'kgas': begin
          n_orders = 29L
          min_order_spacing = 15
        end
        ;K2-band data has 32 orders and 15 pixels minimum between each order
        'k2': begin
          n_orders = 32L
          min_order_spacing = 15
        end
        else: message, 'There are no options set for filter ID '+strtrim(current_filter,2)
      endcase
      
      print, 'Building orders mask for flat ID ['+strtrim(f+1L,2L)+'/'+strtrim(nflat_uniq,2L)+']: '+flat_ids_uniq[f]+'...'
      orders_mask_f = ishell_trace_orders(flats_uniq_cube[*,*,f],orders_structure=orders_structure_f,N_ORDERS=n_orders, MIN_ORDER_SPACING=min_order_spacing,debug=debug_trace_orders)
      save, file=order_mask_file, orders_mask_f, orders_structure_f, n_orders, min_order_spacing, /compress
    endelse
    ;Store data in cubes
    orders_mask_cube[*,*,f] = orders_mask_f
    if ~keyword_set(orders_structure_cube) then begin
      orders_structure_cube = orders_structure_f[0L]
      nan_str, orders_structure_cube
      orders_structure_cube = replicate(orders_structure_cube,max_n_orders,nflat_uniq)
    endif
    orders_structure_cube[0:n_orders-1L,f] = orders_structure_f
  endfor
  
  ;Correct flat fields for fringing
  flats_uniq_cube_corrected = flats_uniq_cube+!values.d_nan
  for f=0L, nflat_uniq-1L do begin
    flat_field_file = flats_dir+'corr_flat_field_'+date_id+'_ID_'+flat_ids_uniq[f]+'.fits'
    lumcorr_flat_field_file = flats_dir+'lumcorr_flat_field_'+date_id+'_ID_'+flat_ids_uniq[f]+'.fits'
    fringe_flat_field_file = flats_dir+'fringe_flat_field_'+date_id+'_ID_'+flat_ids_uniq[f]+'.fits'
    fringe_1d_file = flats_dir+'fringing_1d'+path_sep()+'fringe_1d_'+date_id+'_ID_'+flat_ids_uniq[f]
    fringe_1d_fits_file = flats_dir+'fringe_1d_'+date_id+'_ID_'+flat_ids_uniq[f]+'.fits'
    if file_test(flat_field_file) then begin
      flat_corrected = readfits(flat_field_file,/silent)
    endif else begin
      print, ' Correcting flat field for fringing, flat ID ['+strtrim(f+1L,2L)+'/'+strtrim(nflat_uniq,2L)+']: '+flat_ids_uniq[f]+'...'
      flat_corrected = ishell_flat_fringing(flats_uniq_cube[*,*,f], orders_structure_cube[*,f], orders_mask_cube[*,*,f], $
        CORRECT_BLAZE_FUNCTION=correct_blaze_function_in_flatfield, LUMCORR_FLAT=lumcorr_flat, FRINGING_FLAT=fringing_flat, $
        CORRECT_FRINGING=correct_fringing_in_flatfield, FRINGING_SOLUTION_1D=fringing_solution_1d,fringe_nsmooth=fringe_nsmooth, $
        MODEL_FRINGING=model_fringing)
      
      ;Try to fit a model to the fringing
      ;min_fit_index = 800L
      ;max_fit_index = 1200L
      
      ;=====================================
      ;START OF: NEW METHOD -- All of this should be moved inside of "ishell_flat_fringing.pro" 
      ;=====================================
      
;      ;Extract fringing without horizontal smoothing so that we can bring out the detector patterns
;      flat_corrected = ishell_flat_fringing(flats_uniq_cube[*,*,f], orders_structure_cube[*,f], orders_mask_cube[*,*,f], $
;        CORRECT_BLAZE_FUNCTION=correct_blaze_function_in_flatfield, LUMCORR_FLAT=lumcorr_flat, FRINGING_FLAT=fringing_flat, $
;          CORRECT_FRINGING=correct_fringing_in_flatfield, FRINGING_SOLUTION_1D=fringing_solution_1d,fringe_nsmooth=fringe_nsmooth, $
;          MODEL_FRINGING=model_fringing)
      
;      nrowi = (size(fringing_solution_1d))[2]
;      model_pars = dblarr(4L,nrowi)+!values.d_nan
;      models = dblarr(nx,nrowi)+!values.d_nan
;      nfit = 40
;      for rowi=0L, nrowi-1L do begin & $
;        min_fit_index = 400L & $
;        max_fit_index = 1500L & $
;        ;Data to be fitted
;        fit_y = fringing_solution_1d[min_fit_index:max_fit_index,rowi] & $
;        fit_x = dindgen(n_elements(fit_y))+min_fit_index & $
;        ;Estimate initial fitting paramters
;        fit_par_estim = [weighted_median(abs(fit_y-1),medval=.9),$;Amplitude
;          36.5d0,$;Period
;          0d0,$;Phase
;          -0.0015d0] & $;Period slope
;        fit_par_scatter = fit_par_estim/1d2 & $
;        fit_par_scatter[1] = 5. & $
;        fit_par_scatter[2] = !dpi/10. & $
;        fit_par_scatter[3] = 1d-4 & $
;        fit_par_noise = fit_par_estim#make_array(nfit,value=1d0,/double) + (fit_par_scatter#make_array(nfit,value=1d0,/double))*randomn(seed,n_elements(fit_par_estim),nfit) & $
;        fit_par_noise[*,0] = fit_par_estim & $
;        redchi2s = dblarr(nfit)+!values.d_nan & $
;        fitpars = dblarr(n_elements(fit_par_estim),nfit)+!values.d_nan & $
;        for fiti=0L, nfit-1L do begin & $
;          fit_pari = mpfitfun('ishell_fringing_1d_model',fit_x,fit_y,1d0,fit_par_noise[*,fiti],YFIT=yfit,status=status,err=err,/nan,/quiet) & $
;          redchi2s[fiti] = total((fit_y-yfit)^2,/nan)/double(total(finite(fit_y-yfit))) & $
;          fitpars[*,fiti] = fit_pari & $
;        endfor & $
;        ;Idenfity min chi2 and keep it
;        void = min(redchi2s,wmin) & $
;        fit_par = fitpars[*,wmin] & $
;        ;fit_par = mpfitfun('ishell_fringing_1d_model',fit_x,fit_y,1d0,fit_par_estim,YFIT=yfit,status=status,err=err,/nan,/quiet) & $
;        wset, 0 & $
;        wait,0.1 & $
;        plot,fit_x,fit_y, xtitle=strtrim(rowi,2) & oplot, fit_x, ishell_fringing_1d_model(fit_x,fit_par), col=255 & $
;        wait,0.1 & $
;        model_pars[*,rowi] = fit_par & $
;        models[*,rowi] = ishell_fringing_1d_model(dindgen(nx),fit_par) & $
;      endfor
;      gneg = where(reform(model_pars[0,*]) lt 0, ngneg)
;      if ngneg ne 0L then model_pars[0,gneg] *= -1
;      if ngneg ne 0L then model_pars[2,gneg] += !dpi/2d0
;      model_pars[2,*] = ((model_pars[2,*]+2*!dpi) mod (!dpi*2d0))
;      
;      nsmooth_fringe_illum = 100L
;      fringe_illum_function = smooth(median(max(fringing_solution_1d,dim=2,/nan),nsmooth_fringe_illum),nsmooth_fringe_illum)
;      ;Normalize
;      fringe_illum_function = (fringe_illum_function-1.)/max((fringe_illum_function-1.),/nan)
;      fringe_illum_function[0:nsmooth_fringe_illum] = 0.
;      fringe_illum_function[-nsmooth_fringe_illum:*] = 0.
;      model_with_edges = (fringe_illum_function#make_array((size(fringing_solution_1d))[2],value=1d0,/double))*(models-1.)+1.
;      
;      ;Bring out detector patterns
;      det_patterns = median(fringing_solution_1d/model_with_edges,dim=2)
;      plot,det_patterns-1,yrange=[-.03,.04] & for i=0L, 32L do oplot, i*64+[0,0], [-10,10], col=255 & oplot,det_patterns-1
;      
;      det_patterns_2d = (det_patterns#make_array((size(fringing_solution_1d))[2],value=1d0,/double))
;      
;      ;Smooth the fringing solution w/o detector patterns
;      fringing_no_det_patterns = fringing_solution_1d/det_patterns_2d
;      fringe_nsmooth = 3L
;      for sri=0L, (size(fringing_solution_1d))[2]-1L do $
;        fringing_no_det_patterns[*,sri] = median(fringing_no_det_patterns[*,sri],fringe_nsmooth)
;      stop
      
      ;=====================================
      ;END OF: NEW METHOD
      ;=====================================
      ;flatgas=mrdfits('/Users/gagne/Documents/Data_Repository/RAW/IRTF/iShell/20161107UT_Vega_night2/icm.2016B107.161107.flat.00132.a.fits',0)
      ;flat=mrdfits('/Users/gagne/Documents/Data_Repository/RAW/IRTF/iShell/20161016UT_Vega_night1/icm.2016B107.161016.flat.00042.a.fits',0)
      ;flat_diff = double(flatgas)/double(flat)
      ;bad = where(double(flatgas)/median(double(flatgas)) le 1. or double(flat)/median(double(flat)) le 0)
      ;flat_diff[bad]=0.
      ;flat_corrected = ishell_flat_fringing(flat_diff, orders_structure_cube[*,f], orders_mask_cube[*,*,f], CORRECT_BLAZE_FUNCTION=correct_blaze_function_in_flatfield, LUMCORR_FLAT=lumcorr_flat, FRINGING_FLAT=fringing_flat, CORRECT_FRINGING=correct_fringing_in_flatfield, FRINGING_SOLUTION_1D=fringing_solution_1d,fringe_nsmooth=1L)
      ;CTRL+C then i=0 then run lines 76-99
      ;vf,straight_flat_order
      ;plot,median(straight_flat_order,dim=2)/2.,yrange=[1.9,2.1]-1.
      
      
      ;Any pattern present in all orders is likely related to the detector and should be left in the flat fields.
      if remove_detector_patterns_from_data eq 1 then begin
        
        detector_patterns_block_size = 64L
        
        frac_part = double(nx)/double(detector_patterns_block_size)-long(double(nx)/double(detector_patterns_block_size))
        if frac_part gt 1e-3 then $
          message, ' The detector size is not a multiple of the detector patterns block size !' 
        
        ;Compute the median of detector patterns in moving blocks
        ndpi = long(nx/detector_patterns_block_size)
        median_blocks = dblarr(ndpi)
        median_blocks_image = finite(fringing_solution_1d)*0d0
        for dpi=0L, ndpi-1L do begin & $
          min_hindex = dpi*detector_patterns_block_size & $
          max_hindex = (dpi+1L)*detector_patterns_block_size-1L & $
          median_blocks[dpi] = median(fringing_solution_1d[min_hindex:max_hindex,*]) & $
          median_blocks_image[min_hindex:max_hindex,*] = median_blocks[dpi] & $
        endfor
        
        detector_patterns = median(fringing_solution_1d,dim=2)
        
        ;Remove detector patterns from flat and fringing solutions
        fringing_solution_1d /= (detector_patterns#make_array((size(fringing_solution_1d))[2],value=1d0,/double))
        detector_patterns_ny = detector_patterns#make_array(ny,value=1d0,/double)
        fringing_flat /= detector_patterns_ny
        
        ;Add back detector patterns in flats
        flat_corrected *= detector_patterns_ny
        lumcorr_flat *= detector_patterns_ny
      endif
      
      ;Output 1D fringing solutions as text files
      if ~file_test(flats_dir+'fringing_1d'+path_sep()) then $
        file_mkdir, flats_dir+'fringing_1d'+path_sep()
      for no=0L, (size(fringing_solution_1d))[2]-1L do $
        printuarr, fringe_1d_file+'_ORDER'+strtrim(no+1,2L)+'.txt', fringing_solution_1d[*,no]
      writefits, fringe_1d_fits_file, fringing_solution_1d
      writefits, flat_field_file, flat_corrected
      writefits, lumcorr_flat_field_file, lumcorr_flat
      writefits, fringe_flat_field_file, fringing_flat
    endelse
    
    ;If flat field correction is to be overrided
    if keyword_set(override_flats) then $
      flat_corrected[*] = 1d0
    
    flats_uniq_cube_corrected[*,*,f] = flat_corrected
  endfor
  
  ;Identify all science targets
  g_science = where(object_names ne 'QTH' and strlowcase(object_names) ne 'flat' and strlowcase(do_reduce) eq 'yes' and strpos(strlowcase(object_names),'dark') eq -1, ng_science)
  
  ;Create directories for data output
  previews_dir = output_dir+'previews'+path_sep()
  if ~file_test(previews_dir) then file_mkdir, previews_dir
  spectra_dir = output_dir+'spectra'+path_sep()
  if ~file_test(spectra_dir) then file_mkdir, spectra_dir
  
  ;Determine the IDs that flat files would have for all science exposures.
  ;Basically we want to use the flats intended for that object, with the same filter and the same slit
  science_flat_ids = tcs_obj[g_science]+'|'+data_filters[g_science]+'|'+data_slits[g_science] 
  
  ;Loop on science targets to perform the data reduction
  for sci=0L, ng_science-1L do begin
    print, ' Extracting Exposure ['+strtrim(sci+1L,2L)+'/'+strtrim(ng_science,2L)+']: '+object_names[g_science[sci]]+' !'
    
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
    orders_structure = orders_structure_cube[*,g_flat_id[0L]]
    orders_mask = orders_mask_cube[*,*,g_flat_id[0L]]
    
    ;If all orders of this exposure are already present, skip this exposure altogether
    object_name = object_names[g_science[sci]]
    output_files = spectra_dir+file_basename(data_file,'.fits')+'_OBJ_'+object_name+'_ORD_'+strtrim(orders_structure.order_id,2)+'_spectrum.txt'
    if min(file_test(output_files)) eq 1 then begin
      print, '   Skipping extraction because it already exists...'
      continue
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
    
    ;Determine the number of orders
    good_orders = where(orders_structure.order_id ge 0,n_orders)
    
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
        extraction_data_file = extraction_idl_dir+'extraction_idl_data'+strtrim(orders_structure[i].order_id,2L)+'_'+object_name+'_'+file_basename(data_file,file_ext(data_file))+'.sav'
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
          models_file = models_dir+'models_order'+strtrim(orders_structure[i].order_id,2L)+'_'+object_name+'_'+file_basename(data_file,file_ext(data_file))+'.sav'
          if file_test(models_file) then begin
            ;Make sure the save file contains the updated 2d model
            updated_2d_model = !NULL
            restore, models_file;, best_fit_parameters, functargs, updated_2d_model, badpix_mask
            if keyword_set(best_fit_parameters) then $
              full_model_polynomials[*,i] = best_fit_parameters
            if keyword_set(updated_2d_model) then begin
              ;Identify order positions and patch the 2D model at the order position
              g_order = where(orders_mask eq orders_structure[i].ORDER_ID, ng_order)
              full_2d_model[g_order] = updated_2d_model[g_order]
            endif
          endif
        endif
        continue
      endif
      
      ;Create a list of y positions for that orders trace
      y_position_order_i_estimate = poly(xarr,orders_structure[i].mid_coeffs)

      ;Define the height of the order
      height = ceil(orders_structure[i].height)

      ;Create an image where only the relevant order is seen
      im_order = data_im
      g_order = where(orders_mask eq orders_structure[i].ORDER_ID, ng_order, complement=bad_order, ncomplement=nbad_order)
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
        trace_profile[0L:mask_trace_edges-1L] = 0d0
        trace_profile[-mask_trace_edges:*] = 0d0
      endif
      
      ;Re-scale the trace profile
      trace_profile -= weighted_median(trace_profile,medval=.05)
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
      trace_profiles_dir = output_dir+'trace_profiles'+path_sep()
      if ~file_test(trace_profiles_dir) then file_mkdir, trace_profiles_dir
      trace_profiles_file = trace_profiles_dir+'trace_profile_order'+strtrim(orders_structure[i].order_id,2L)+'_'+object_name+'_'+file_basename(data_file,file_ext(data_file))+'.png'
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
      left_edge_ypos = poly(xarr,orders_structure[i].left_coeffs)
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
        
        ;Recreate a 2D profile
        message, ' Do not trust these best fit parameters for the trace shape, just use the model to flag bad pixels. This is what caused continuum to crap out sometimes', /continue
        refined_again_coefficients = best_fit_parameters[0:ndegree_poly_fit-1L]
        sigma_coefficients = best_fit_parameters[ndegree_poly_fit:*]
        sigma = poly(xmap, sigma_coefficients)
        gauss_y = poly(xmap, refined_again_coefficients)
        profile_2d = exp(-(ymap-gauss_y)^2/(2d0*sigma^2))
        
        ;Store the coefficients in the full-orders array of coefficients
        full_model_polynomials[*,i] = best_fit_parameters
        
        ;Re-do the optimal extraction with the refined parameters
        opt_spectrum = general_optimal_extract_pmassey(im_order, sky_2d, profile_2d*badpix_mask, eff_read_noise, gain, ESP=err_opt_spectrum)
        
        ;Refine the 2D model even more for output. This can be useful for debugging.
        opt_spectrum_2d = opt_spectrum#make_array(ny,value=1d0,/double)
        functargs.spectrum = opt_spectrum_2d
        updated_2d_model = fit_2d_curved_trace(xmap, ymap, best_fit_parameters, _EXTRA=functargs, /nobpixpatch)
        wbadpix = where(badpix_mask eq 0, nwbadpix)
        if nwbadpix ne 0L then updated_2d_model[wbadpix] = !values.d_nan
        
        ;Save the refined 2D model in IDL save format
        models_file = models_dir+'models_order'+strtrim(orders_structure[i].order_id,2L)+'_'+object_name+'_'+file_basename(data_file,file_ext(data_file))+'.sav'
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
      extraction_data_file = extraction_idl_dir+'extraction_idl_data'+strtrim(orders_structure[i].order_id,2L)+'_'+object_name+'_'+file_basename(data_file,file_ext(data_file))+'.sav'
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
        a = plot(buffer=1,opt_spectrum_norm,xtitle='Pixels (Order '+strtrim(orders_structure[i].order_id,2)+')',ytitle='Relative flux',thick=2,color='black',margin=[.07,.12,.02,.02],xthick=2,ythick=2,yticklen=.015,xticklen=.015,dimensions=[1000,512],xrange=[gfinite[0L]-3L,gfinite[-1L]+3L],name='Optimal')
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
        yrange = a.yrange
        yrange += [-1,1]*(max(yrange)-min(yrange))*.05
        ;Limit max yrange to 1.2
        yrange = [yrange[0L],yrange[1L]<1.2]
        a.yrange = yrange
        ;Save the figure
        a.save, previews_dir+file_basename(data_file,'.fits')+'_OBJ_'+object_name+'_ORD_'+strtrim(orders_structure[i].order_id,2)+'_preview.png'
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
      sxaddpar, header, 'NOFLATS', override_flats, '/override_flats flag to skip flat field correction'
      
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
        full_model_file = models_dir+'full_model_'+object_name+'_'+file_basename(data_file,file_ext(data_file))+'.fits'
        writefits, full_model_file, full_2d_model, /silent
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
        full_residuals_file = residuals_dir+'full_residuals'+'_'+object_name+'_'+file_basename(data_file,file_ext(data_file))+'.fits'
        writefits, full_residuals_file, data_im-full_2d_model, /silent
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
        yrange = pai.yrange
        yrange += [-1,1]*(max(yrange)-min(yrange))*.05
        ;Limit max yrange to 1.2
        yrange = [yrange[0L],yrange[1L]<1.2]
        pai.yrange = yrange
      endfor
    endif
    
    full_orders_spectra_file = previews_dir+'full_orders_spectra_preview_'+file_basename(data_file,'.fits')+'_OBJ_'+object_name+'.png'
    window_id.save, full_orders_spectra_file
    pai.close
    
  endfor
  
  print, ' All Done !'
  
End