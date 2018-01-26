Function kill_value, vec_in, tokill, toreplace_in, NAN=nan
  ;Remplace toutes les valeurs "tokill" par "toreplace_in" dans un vecteur donne.
  ; Par defaut, toreplace_in est la mediane du vecteur ou un string vide.
  ; /nan permet de remplacer des valeurs NaN. Ex: print, kill_value(vec,/nan) : Remplace les NaN par la mediane. 
  ; toreplace_in = 'zero' -> 0
  compile_opt hidden
  vec = vec_in
  if keyword_set(toreplace_in) then toreplace = toreplace_in
  if ~keyword_set(toreplace_in) then toreplace = size(vec,/tname) eq 'STRING' ? '' : median([vec])
  if keyword_set(toreplace_in) then begin
    if size(toreplace_in,/tname) eq 'STRING' and size(vec,/tname) ne 'STRING' then $
      if toreplace_in eq 'zero' then toreplace = 0
  endif
  kill_ind = keyword_set(nan) ? where(~finite([vec]),nkill) : where([vec] eq tokill,nkill)
  if nkill eq n_elements(vec) and size(vec,/tname) ne 'STRING' then return, -1
  if nkill gt 0 then vec[kill_ind] = toreplace
  return, vec
End