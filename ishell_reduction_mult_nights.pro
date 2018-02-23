Pro ishell_reduction_mult_nights, root_data_dir, output_dir_root, DEBUG_TRACE_ORDERS=debug_trace_orders, DO_DARK_SUBTRACTION=do_dark_subtraction, CORRECT_FRINGING_IN_FLATFIELD=correct_fringing_in_flatfield, MODEL_FRINGING=model_fringing, REMOVE_DETECTOR_PATTERNS_FROM_DATA=remove_detector_patterns_from_data, OVERRIDE_FLATS=override_flats, MODEL_REFINEMENT_CURVATURE=model_refinement_curvature, CORRECT_BLAZE_FUNCTION_IN_FLATFIELD=correct_blaze_function_in_flatfield, NTHREADS=nthreads

	if ~keyword_set(root_data_dir) then message, ' ERROR: MUST INCLUDE ROOT_DATA_DIR AS FIRST ARGUMENT.'
	if ~keyword_set(output_dir_root) then message, ' ERROR: MUST INCLUDE OUTPUT_DIR_ROOT AS SECOND ARGUMENT.'

	; Only grab the directories, just in case
	nights = file_search(root_data_dir, COUNT=n_nights, /FULLY_QUALIFY_PATH, /TEST_DIRECTORY)

	; Simultaneously reduce these nights instead of parallelizing over individual spectra
	; NOTE: If the number of nights < the number of available threads,
	; one should instead manually reduce these nights, and not use this!
	;command = 'ishell_reduction_master, nights[i], output_dir_root, DEBUG_TRACE_ORDERS=debug_trace_orders, DO_DARK_SUBTRACTION=do_dark_subtraction, CORRECT_FRINGING_IN_FLATFIELD=correct_fringing_in_flatfield, MODEL_FRINGING=model_fringing, REMOVE_DETECTOR_PATTERNS_FROM_DATA=remove_detector_patterns_from_data, OVERRIDE_FLATS=override_flats, MODEL_REFINEMENT_CURVATURE=model_refinement_curvature, CORRECT_BLAZE_FUNCTION_IN_FLATFIELD=correct_blaze_function_in_flatfield, NTHREADS=1'
	;varnames = ['nights', 'output_dir_root']
	;split_for, 0, n_nights-1, commands=command, varnames=varnames

	print, nights
	stop

	for i=0, n_nights-1 do $
		ishell_reduction_master, nights[i], output_dir_root, DEBUG_TRACE_ORDERS=debug_trace_orders, DO_DARK_SUBTRACTION=do_dark_subtraction, CORRECT_FRINGING_IN_FLATFIELD=correct_fringing_in_flatfield, MODEL_FRINGING=model_fringing, REMOVE_DETECTOR_PATTERNS_FROM_DATA=remove_detector_patterns_from_data, OVERRIDE_FLATS=override_flats, MODEL_REFINEMENT_CURVATURE=model_refinement_curvature, CORRECT_BLAZE_FUNCTION_IN_FLATFIELD=correct_blaze_function_in_flatfield, NTHREADS=nthreads

end