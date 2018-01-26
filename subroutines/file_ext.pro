Function file_ext, file
  ;Outputs the extention of a file
  ;WRITTEN, Jonathan Gagne, April 25 2014.
  forward_function strmid_arr, strpos_arr
  ff = file_basename(file)
  return, strmid_arr(ff,strpos_arr(ff,'.',/reverse_search))
End