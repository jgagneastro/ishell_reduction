Pro ishell_reduction_both_fringing_modes
  ;This code will run "ishell_reduction_master" twice, with and without fringing correction in the flat fields
  
  ;/// Adjustable parameters ///
  
  ;The path of the night which will be reduced
  if ~keyword_set(data_path) then $
    data_path = '/Volumes/bryson/iSHELL/Data/Raw/20171023UT/'

  ;The path where the data reduction products will be stored
  if ~keyword_set(output_dir_root) then $
    output_dir_root = '/Volumes/bryson/iSHELL/redux/'

  ;Whether or not to do debugging for trace order detection
  if debug_trace_order eq !NULL then $
    debug_trace_orders = 1

  ;Whether or not darks should be subtracted
  if do_dark_subtraction eq !NULL then $
    do_dark_subtraction = 0L
  
  ;/// End of: Adjustable parameters ///
  
  ;If the log file does not exist, create it
  logfile = file_search(output_dir_root+'*.log')
  if logfile[0] eq '' then begin
    ishell_reduction_master
  endif
  
  ;If it exists, run the reduction code twice with appropriate output directories
  if logfile[0] ne '' then begin
    ishell_reduction_master, data_path, output_dir_root+'fringing_corrected'+path_sep(), DEBUG_TRACE_ORDER=debug_trace_order, DO_DARK_SUBTRACTION=do_dark_subtraction, CORRECT_FRINGING_IN_FLATFIELD=1
    ishell_reduction_master, data_path, output_dir_root+'fringing_not_corrected'+path_sep(), DEBUG_TRACE_ORDER=debug_trace_order, DO_DARK_SUBTRACTION=do_dark_subtraction, CORRECT_FRINGING_IN_FLATFIELD=0
  endif
  
End