Function strinvert, str
;  Cette fonction inverse un string ou un array de strings
  compile_opt hidden
  forward_function reverse
  on_error, 2
  
  ;Mode array
  if isa(str,/array) then begin
    result = strarr(n_elements(str))
    for i=0, n_elements(str)-1 do $
      result[i] = strinvert(str[i])
    return, result
  endif
  
  return, strjoin(reverse(string(transpose(byte(str)))))
End