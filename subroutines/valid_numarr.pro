Function valid_numarr, string, value, INTEGER=integer
  ;Comme valid_num, mais fonctionne pour les arrays
  
  compile_opt hidden
  on_error, 2
  forward_function valid_num
  
  if ~isa(string,/array) then return, valid_num(string, value, INTEGER=integer)
  
  ns = n_elements(string)
  ret = intarr(ns)
  if keyword_set(integer) then value = intarr(ns) else value = dblarr(ns)
  for i=0, ns-1 do begin
    ret[i] = valid_num(string[i], valuei, INTEGER=integer)
    if keyword_set(valuei) then $
      value[i] = temporary(valuei) else $
      if ~keyword_set(integer) then value[i] = !values.f_nan
  endfor
  return, ret
  
End