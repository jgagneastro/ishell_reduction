Pro nan_str, str
  ;This function was written by J. Gagne to fill up an existing structure with NaN values of the appropriate data type.
  forward_function create_nan
  if n_elements(str) gt 1 then message, ' La structure ne peut pas etre un array ! Utilisez le resultat avec replicate() pour avoir un array.'
  tags = tag_names(str)
  ntags = n_elements(tags)
  for i=0, ntags-1 do begin
    void = execute('nel = n_elements(str.'+tags[i]+')')
    void = execute('issubstr = size(str.'+tags[i]+',/type) eq 8L')
    if issubstr then begin
      s0 = str.(i)
      nan_str, s0
      str.(i) = temporary(s0)
      continue
    endif
    if nel eq 1 then $
      void = execute('str.'+tags[i]+' = create_nan(size(str.'+tags[i]+',/type))') else $
      void = execute('str.'+tags[i]+' = replicate(create_nan(size(str.'+tags[i]+',/type)),nel)')
  endfor
End