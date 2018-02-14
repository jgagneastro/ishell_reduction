Pro ishell_reduction_master, data_path, output_dir_root, DEBUG_TRACE_ORDERS=debug_trace_orders, DO_DARK_SUBTRACTION=do_dark_subtraction, $
  CORRECT_FRINGING_IN_FLATFIELD=correct_fringing_in_flatfield, MODEL_FRINGING=model_fringing, $
  REMOVE_DETECTOR_PATTERNS_FROM_DATA=remove_detector_patterns_from_data, OVERRIDE_FLATS=override_flats, $
  MODEL_REFINEMENT_CURVATURE=model_refinement_curvature, CORRECT_BLAZE_FUNCTION_IN_FLATFIELD=correct_blaze_function_in_flatfield, $
  NTHREADS=nthreads
  ;Code version history
  ; The code started at 0.9, between 0.9 and 1.0 only minor bugs were fixed and more diagnostic outputs were added
  ; Version 1.0: First stable version (J. Gagne)
  ; Version 1.1: Added outputs of polynomial coefficients for Ytrace pos and seeing versus X pixel in ASCII format (J. Gagne), March 3, 2017
  ; Version 1.2: Added N_orders keyword for compatibility with K band, January 26, 2018
  ; Version 1.3: More fixes for mixed bands in a single night, added an option to remove column-dependent detector patterns from the data, January 31, 2018
  ; Version 1.4: Added option to model fringing to separate it from the column-dependent detector patterns. February 4, 2018
  ; Version 1.5: Added support for J2 band, .gz data, fixed problems with continuum in the optimal spectrum and improved bad pixel rejection. Feb. 6, 2018
  ; Version 1.6: Fixed a problem with fringing residuals left in flat fields. Feb 8, 2018
  ; Version 1.7: Added support for split_for with NTHREADS keyword
  ; 
  ;Planned modifications:
  ; - Use A star to derive Blaze function: reduce Vega with lumcorr to do that
  ; - The file name for detector patterns preview must have flat ID in it
  ; - Try optimally extracting the non-lumcorr flat field and see if the resulting Blaze function works better.
  
  ;Code version for headers
  code_version = 1.7
  
  ;List of subroutines
  forward_function readfits, ishell_trace_orders, ishell_flat_fringing, interpol2, weighted_median, $
    horizontal_median, vertical_median, sxpar, whiten_color, sxpar_mul, rtrim, readcol, printuarr, $
    match_unsorted, strtrim_multiple, nan_str, mpfit2dfun, fit_2d_curved_trace, $
    general_optimal_extract_pmassey, correct_bad_pixels_posterior, sxaddpar, strkill, $
    ishell_extract_exposure
  
  ;/// Adjustable parameters ///
  
  ;The path of the night which will be reduced
  if ~keyword_set(data_path) then $
    data_path = '/Users/gagne/Documents/Data_Repository/RAW/IRTF/iShell/20171023UT/'
    ;data_path = '/Users/gagne/Documents/Data_Repository/RAW/IRTF/iShell/20161016UT_Vega_night1/'
    ;data_path = '/Users/gagne/Documents/Data_Repository/RAW/IRTF/iShell/20161107UT_Vega_night2/'
    ;data_path = '/Users/gagne/Documents/Data_Repository/RAW/IRTF/iShell/20161107UT_Vega_night2/'
    ;data_path = '/Users/gagne/Documents/Data_Repository/RAW/IRTF/iShell/20161016UT_Vega_No_Gas/'
  
  ;The number of CPU threads to use when extracting exposures
  if ~keyword_set(nthreads) then $
    nthreads = 1L
  
  ;The path where the data reduction products will be stored
  if ~keyword_set(output_dir_root) then $
    output_dir_root = '/Users/gagne/Documents/Data_Repository/RAW/IRTF/iShell/redux/'
  
  ;Add trailing path separators if needed
  if strmid(data_path,0,1,/reverse_offset) ne path_sep() then $
    data_path += path_sep()
  if strmid(output_dir_root,0,1,/reverse_offset) ne path_sep() then $
    output_dir_root += path_sep()
  
  ;Make sure that paths end with a path separator
  if strmid(data_path,0,1) ne path_sep() then $
    data_path += path_sep()
  if strmid(output_dir_root,0,1) ne path_sep() then $
    output_dir_root += path_sep()
  
  ;Whether or not to do debugging for trace order detection
  if debug_trace_orders eq !NULL then $
    debug_trace_orders = 0
  
  ;Whether or not to avoid flat fields correction
  if ~keyword_set(override_flats) then $
    override_flats = 0
  
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
  
  ;Whether to use 2D modelling of the trace to refine trace curvature.
  ; This also allows a linear seeing dependence over spatial pixel position
  ; However this was found to cause bad model curvatures when the PSF is weird
  ; enough on low SNR data, hence we set it off by default. 
  if model_refinement_curvature eq !NULL then $
    model_refinement_curvature = 0L
  
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
  if correct_blaze_function_in_flatfield eq !NULL then $
    correct_blaze_function_in_flatfield = 1
  
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
  
  ;List of available filters with their relevant settings
  filters_names =               ['KGAS', 'K2', 'J2']
  filters_n_orders =            [29L,     32L, 38L]
  filters_min_order_spacing =   [15L,     15L, 12L]
  filters_unequal_order_edges = [1,       1,   1]
  ;min_order_spacing is the minimal possible space in pixels between any two orders
  ;"unequal_order_edges" must be set to 1 if there is 1 more left edges than right edges
  ;  in the order detection plot (use debug_trace_orders to view that plot) *and* you
  ;  want to extract even the last order to the right
  
  max_n_orders = max(filters_n_orders,/nan)
  
  ;/// End of: adjustable parameters ///
  
  ;Add the night ID after the output directory
  date_id = file_basename(data_path)
  output_dir = output_dir_root+date_id+path_sep()
  
  ;Define directory for trace profiles
  trace_profiles_dir = output_dir+'trace_profiles'+path_sep()
  
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
  fits_data = file_search(data_path+'*.fits*')
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
    bad = where(strpos(strlowcase(object_names),'dark') ne -1 or strpos(strlowcase(object_names),'thar') ne -1 or (strpos(strlowcase(object_names),'qth') ne -1 and strlowcase(gascell_position) eq 'in'), nbad)
    if nbad ne 0L then do_reduce[bad] = 'No'
    
    ;Create an crude log file and put some header information in it 
    printuarr, logfile, file_basename(fits_data), object_names, rtrim(integration_times,3), data_comments, tcs_obj, $
      gascell_position, data_filters, data_slits, do_reduce, title=[';File','Object','ExpTime','Comment','TCS Object','GasCell','Filter','Slit','Reduce'], $
      /justify, symbol='|', /new
    
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
  strtrim_multiple, filenames_log, object_names, data_comments, tcs_obj, gascell_position, data_filters, data_slits, do_reduce
  
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
  
  darks_ids_uniq = ''
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
        writefits, darks_dir+'dark_medcomb_'+date_id+'_TEXP'+darks_IDs_uniq[f]+'s.fits.gz', darks_uniq_cube[*,*,f], /compress
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
      ;If there is only one flat, do not median combine
      if (size(flat_f_cube))[0] eq 2 then $
        flats_uniq_cube[*,*,f] = flat_f_cube else $
        flats_uniq_cube[*,*,f] = median(flat_f_cube,dim=3)
    endfor
    save, flats_uniq_cube, nx, ny, flat_ids_uniq, nflat_uniq, file=comb_flats_file, /compress
  endelse
  
  ;Build all orders masks
  orders_mask_cube = dblarr(nx,ny,nflat_uniq)+!values.d_nan
  for f=0L, nflat_uniq-1L do begin
    ;Verify whether this was already done and saved to disk
    order_mask_dir = output_dir+'order_masks'+path_sep()
    if ~file_test(order_mask_dir) then file_mkdir, order_mask_dir
    order_mask_file = order_mask_dir+'order_mask_'+date_id+'_ID_'+flat_ids_uniq[f]+'.sav'
    if file_test(order_mask_file) then begin
      restore, order_mask_file;, orders_mask_f, orders_structure_f, n_orders, min_order_spacing
    endif else begin
      
      ;Determine the relevant set-up for the current filter
      current_filter = (strsplit(flat_ids_uniq[f],'|',/extract))[1]
      g_filter = (where(strlowcase(filters_names) eq strlowcase(current_filter), ng_filter))[0]
      if ng_filter ne 1 then $
        message, ' Could not find 1 unique match to filter "'+current_filter+'" in the list of known filters (variable filters_names) !'
      
      n_orders = filters_n_orders[g_filter]
      min_order_spacing = filters_min_order_spacing[g_filter]
      unequal_order_edges = filters_unequal_order_edges[g_filter]
      
      print, 'Building orders mask for flat ID ['+strtrim(f+1L,2L)+'/'+strtrim(nflat_uniq,2L)+']: '+flat_ids_uniq[f]+'...'
      orders_mask_f = ishell_trace_orders(flats_uniq_cube[*,*,f],orders_structure=orders_structure_f,N_ORDERS=n_orders, $
        MIN_ORDER_SPACING=min_order_spacing,debug=debug_trace_orders,unequal_order_edges=unequal_order_edges)
      save, file=order_mask_file, orders_mask_f, orders_structure_f, n_orders, min_order_spacing, /compress
      
      ;Save order mask
      writefits, order_mask_dir+'order_mask_'+date_id+'_ID_'+flat_ids_uniq[f]+'.fits.gz', orders_mask_f, /compress
      
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
    flat_field_file = flats_dir+'corr_flat_field_'+date_id+'_ID_'+flat_ids_uniq[f]+'.fits.gz'
    lumcorr_flat_field_file = flats_dir+'lumcorr_flat_field_'+date_id+'_ID_'+flat_ids_uniq[f]+'.fits.gz'
    fringe_flat_field_file = flats_dir+'fringe_flat_field_'+date_id+'_ID_'+flat_ids_uniq[f]+'.fits.gz'
    fringe_1d_file = flats_dir+'fringing_1d'+path_sep()+'fringe_1d_'+date_id+'_ID_'+flat_ids_uniq[f]
    fringe_1d_fits_file = flats_dir+'fringe_1d_'+date_id+'_ID_'+flat_ids_uniq[f]+'.fits.gz'
    if file_test(flat_field_file) then begin
      flat_corrected = readfits(flat_field_file,/silent)
    endif else begin
      if keyword_set(override_flats) then begin
        flat_corrected = flats_uniq_cube[*,*,f]
      endif else begin
        print, ' Correcting flat field for fringing, flat ID ['+strtrim(f+1L,2L)+'/'+strtrim(nflat_uniq,2L)+']: '+flat_ids_uniq[f]+'...'
        flat_corrected = ishell_flat_fringing(flats_uniq_cube[*,*,f], orders_structure_cube[*,f], orders_mask_cube[*,*,f], $
          CORRECT_BLAZE_FUNCTION=correct_blaze_function_in_flatfield, LUMCORR_FLAT=lumcorr_flat, FRINGING_FLAT=fringing_flat, $
          CORRECT_FRINGING=correct_fringing_in_flatfield, FRINGING_SOLUTION_1D=fringing_solution_1d,fringe_nsmooth=fringe_nsmooth, $
          MODEL_FRINGING=model_fringing,DETECTOR_PATTERNS=detector_patterns,MODELS=models, FLAT_ILLUMINATION=flat_illumination)
      
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
          printuarr, fringe_1d_file+'_ORDER'+strtrim(no+1,2L)+'.txt', fringing_solution_1d[*,no], /new
        writefits, fringe_1d_fits_file, fringing_solution_1d, /compress
        writefits, flat_field_file, flat_corrected, /compress
        writefits, lumcorr_flat_field_file, lumcorr_flat, /compress
        writefits, fringe_flat_field_file, fringing_flat, /compress
        
        ;Output illumination function
        save, flat_illumination, file=flats_dir+'flat_illumination_'+date_id+'_ID_'+flat_ids_uniq[f]+'.sav', /compress
        
        ;Output model fringing and detector patterns model_fringing = 1
        if model_fringing eq 1 then begin
          writefits, flats_dir+'fringe_model_2d_'+date_id+'_ID_'+flat_ids_uniq[f]+'.fits.gz', models, /compress
          printuarr, flats_dir+'detector_patterns_'+date_id+'_ID_'+flat_ids_uniq[f]+'.txt', detector_patterns, /new
          
          ;Output 1D fringing models as text files
          if ~file_test(flats_dir+'fringing_1d'+path_sep()) then $
            file_mkdir, flats_dir+'fringing_1d'+path_sep()
          fringe_models_file = flats_dir+'fringing_1d'+path_sep()+'model_fringe_1d_'+date_id+'_ID_'+flat_ids_uniq[f]
          for no=0L, (size(fringing_solution_1d))[2]-1L do $
            printuarr, fringe_models_file+'_ORDER'+strtrim(no+1,2L)+'.txt', models[*,no], /new
          
          ;Make a figure with detector patterns
          plpat1 = plot(detector_patterns, xtitle='Pixel position',xticklen=.015,yticklen=.015,$
            margin=[.14,.12,.03,.025],font_size=16,thick=2,xthick=2,ythick=2,$
            xrange=[-1L,nx+1L],ytitle='Persistent flux bias',/buffer)
          for no=0L, 32L do $
            plpati = plot(/overplot,no*64+[0,0], [-10,10], color='red',linestyle='--',yrange=plpat1.yrange)
          plpat1.order,/bring_to_front
          plpat1.save, flats_dir+'preview_detector_patterns.png'
          plpat1.close
        endif
      endelse
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
  
  ;Explicit orders structure for compatibility with split_for
  order_ids = orders_structure_cube.order_id
  order_heights = orders_structure_cube.height
  order_left_coeffs = orders_structure_cube.left_coeffs
  order_mid_coeffs = orders_structure_cube.mid_coeffs
  
  ;Do a normal loop to allow output in the console when NTHREADS < 2
  if nthreads lt 2L then begin
    
    ;Loop on science targets to perform the data reduction
    for sci=0L, ng_science-1L do begin
      print, ' Extracting Exposure ['+strtrim(sci+1L,2L)+'/'+strtrim(ng_science,2L)+']: '+object_names[g_science[sci]]+' !'
      
      ishell_extract_exposure, fits_data, tcs_obj, object_names, integration_times, $; Information from observing log
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
      
    endfor
  
  ;Otherwise use split_for
  endif else begin
    
    
    ;The list of variables that need to be passed to the CPU threads
    ; sci does not go in varnames explicitly
    varnames = ['fits_data','tcs_obj','object_names','integration_times',$
      'darks_ids_uniq','flat_ids_uniq','n_orders','g_science','ng_science','science_flat_ids','nx','ny',$
      'edge_npix_mask','edge_npix_avoid','sizelag','window_lag_poly_fit','ndegree_sigma_fit','ndegree_poly_fit','nrows_sky','quartile_fraction_for_norm',$
      'bad_pixels_threshold','second_bad_pixels_threshold','read_noise','dark_current','gain',$
      'generate_full_orders_spectra_figures','generate_residuals_fits_images','generate_full_orders_trace_profile_figures',$
      'spectra_dir','output_dir','trace_profiles_dir','previews_dir',$
      'correct_blaze_function_in_flatfield','correct_fringing_in_flatfield','model_fringing','code_version','remove_detector_patterns_from_data','override_flats',$
      'flats_uniq_cube_corrected','orders_mask_cube',$
      'order_ids','order_heights','order_left_coeffs','order_mid_coeffs',$
      'do_dark_subtraction','do_trace_2d_fit','mask_trace_edges','refine_poly_drop_pixels',$
      'model_refinement_curvature','generate_individual_orders_trace_profile_figures',$
      'generate_individual_orders_spectra_figures']
    
    ;The command that will be launched on each CPU thread (do not include comments or line breaks !)
    commands = ['ishell_extract_exposure, fits_data, tcs_obj, object_names, integration_times, '+$
      'darks_ids_uniq, flat_ids_uniq, n_orders, g_science, ng_science, sci, science_flat_ids, nx, ny, '+$
      'edge_npix_mask, edge_npix_avoid, sizelag, window_lag_poly_fit, ndegree_sigma_fit, ndegree_poly_fit, nrows_sky, quartile_fraction_for_norm, '+$
      'bad_pixels_threshold, second_bad_pixels_threshold, read_noise, dark_current, gain, '+$
      'generate_full_orders_spectra_figures, generate_residuals_fits_images, generate_full_orders_trace_profile_figures, '+$
      'spectra_dir, output_dir, trace_profiles_dir, previews_dir, '+$
      'correct_blaze_function_in_flatfield,correct_fringing_in_flatfield,model_fringing,code_version,remove_detector_patterns_from_data,override_flats, '+$
      'flats_uniq_cube_corrected, orders_mask_cube, '+$
      'order_ids, order_heights, order_left_coeffs, order_mid_coeffs, '+$
      'DO_DARK_SUBTRACTION=do_dark_subtraction, DO_TRACE_2D_FIT=do_trace_2d_fit, MASK_TRACE_EDGES=mask_trace_edges, REFINE_POLY_DROP_PIXELS=refine_poly_drop_pixels, '+$
      'MODEL_REFINEMENT_CURVATURE=model_refinement_curvature, GENERATE_INDIVIDUAL_ORDERS_TRACE_PROFILE_FIGURES=generate_individual_orders_trace_profile_figures, '+$
      'GENERATE_INDIVIDUAL_ORDERS_SPECTRA_FIGURES=generate_individual_orders_spectra_figures']
    
      split_for, 0L, ng_science-1L, $
        commands=commands, NSPLIT=nthreads, ctvariable_name='sci',varnames=varnames
      
  endelse
  
  print, ' All Done !'
  
End
;"normal run": ishell_reduction_master
;"simple block filter" ishell_reduction_master,model_fringing=0,remove_detector_patterns_from_data=1
;"noflats": ishell_reduction_master,/override_flats,model_fringing=0,correct_fringing_in_flatfield=0
;"nofringing_corr_": ishell_reduction_master,model_fringing=0,correct_fringing_in_flatfield=0