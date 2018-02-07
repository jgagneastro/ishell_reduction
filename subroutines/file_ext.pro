Function file_ext, file, FITS=fits
  ;Outputs the extention of a file
  ;WRITTEN, Jonathan Gagne, April 25 2014.
  ;Added FITS keyword to output .fits.gz not just .gz, Jonathan Gagne, Feb 6, 2018
  forward_function strmid_arr, strpos_arr
  ff = file_basename(file)
  ext = strmid_arr(ff,strpos_arr(ff,'.',/reverse_search))
  if keyword_set(fits) then begin
    bad = where(ext eq '.gz', nbad)
    if nbad ne 0L then $
      ext[bad] = '.fits.gz'
  endif
  return, ext
End