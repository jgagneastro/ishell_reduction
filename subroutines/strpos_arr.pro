Function strpos_arr, strin, expin, _extra=ex
;  Comme strpos mais avec "expin" en array
  compile_opt hidden
  on_error, 2
  
  ;Pour éviter d'altérer str et exp
  str = strin
  exp = expin
  
  ;Cas similaire a strpos
  if ~isa(exp,/array) then return, strpos(str,exp, _extra=ex)
  
  ;Cas avec deux arrays
  sortie = intarr(n_elements(exp))
  if isa(str,/array) and isa(exp,/array) then begin
    if n_elements(str) ne n_elements(exp) then $
      message, ' Erreur. EXP et STR doivent etre de la meme taille !'
    for i=0, n_elements(str)-1 do $
      sortie[i] = strpos_arr(str[i],exp[i], _extra=ex)
  endif
  
  ;Cas avec seulement EXP en array
  if ~isa(str,/array) and isa(exp,/array) then begin
    str = strarr(n_elements(exp))+str[0]
    for i=0, n_elements(str)-1 do $
      sortie[i] = strpos_arr(str[i],exp[i], _extra=ex)
  endif
  
  return, sortie
End