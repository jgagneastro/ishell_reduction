Pro ishell_astar_debug
  dir = '/Users/gagne/Documents/Data_Repository/RAW/IRTF/iShell/20161107UT_Vega_night2/'
  redux0 = '/Users/gagne/Documents/Data_Repository/RAW/IRTF/iShell/redux/'
  reduxdir = redux0+file_basename(dir)+path_sep()
  if file_test(reduxdir) then stop
  
  ;Create log
  ishell_reduction_master, dir, redux0
  
  ;Run
  ishell_reduction_master, dir, redux0
  
  ;Rename output
  file_move,reduxdir,file_dirname(reduxdir)+path_sep()+'modelling_fringing_'+file_basename(reduxdir)
  
  ishell_reduction_master, dir, redux0,model_fringing=0,remove_detector_patterns_from_data=1
  ishell_reduction_master, dir, redux0,model_fringing=0,remove_detector_patterns_from_data=1
  file_move,reduxdir,file_dirname(reduxdir)+path_sep()+'simple_block_filter_'+file_basename(reduxdir)
  
  ishell_reduction_master, dir, redux0,/override_flats,model_fringing=0,correct_fringing_in_flatfield=0
  ishell_reduction_master, dir, redux0,/override_flats,model_fringing=0,correct_fringing_in_flatfield=0
  file_move,reduxdir,file_dirname(reduxdir)+path_sep()+'noflats_'+file_basename(reduxdir)
  
  ishell_reduction_master, dir, redux0,model_fringing=0,correct_fringing_in_flatfield=0
  ishell_reduction_master, dir, redux0,model_fringing=0,correct_fringing_in_flatfield=0
  file_move,reduxdir,file_dirname(reduxdir)+path_sep()+'nofringing_corr_'+file_basename(reduxdir)
  
  print, ' Done !'
  stop
  
End