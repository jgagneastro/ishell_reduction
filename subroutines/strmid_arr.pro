function strmid_arr, str, dep, len, REVERSE_OFFSET=reverse_offset
  ;Cette fonction est semblable a strmid, mais elle permet de fonctionner avec "dep" et "len" comme des arrays de la meme longueur que "str"
  
  ;Evite de toujours afficher 'Compiled module'
  ;compile_opt hidden  
  on_error, 2
  
  store = strarr(n_elements(str))
  if ~keyword_set(len) then len = strlen(str)
  if total(len) eq 0 then return, str
  if not isa(dep,/array) and not isa(len,/array) then begin
    if len lt 0 then return, strmid(str, dep) else $
      return, strmid(str, dep, len)
  endif
  
  ;Maximum de memoire
  maxm = 32000l
  if n_elements(str)-1 gt maxm then begin
    limit = maxm
    doagain = 1
  endif else begin
    limit = n_elements(str)-1
    doagain = 0
  endelse
  
  for i=0, limit do begin
    if isa(dep,/array) then a = dep[i] else a = dep
    if isa(len,/array) then b = len[i] else b = len
    if b lt 0 then b = strlen(str[i])
    store[i] = strmid(str[i],a,b, REVERSE_OFFSET=reverse_offset)
  endfor
  
  ;-- Les lignes suivantes sont utiles dans les cas d'arrays TRES gros
  if doagain then begin
    limit2 = n_elements(str)-1-maxm
    if limit2 gt maxm then begin
      limit2 = maxm
      doagain2 = 1
    endif else begin
      doagain2 = 0
    endelse
    
    for i=0, limit2 do begin
      if isa(dep,/array) then a = dep[i+maxm] else a = dep
      if isa(len,/array) then b = len[i+maxm] else b = len
      if b lt 0 then b = strlen(str[i+maxm])
      store[i+maxm] = strmid(str[i+maxm],a,b, REVERSE_OFFSET=reverse_offset)
    endfor
  endif else doagain2 = 0
  
  if doagain2 then begin
    limit3 = n_elements(str)-1-2*maxm
    for i=0, limit3 do begin
      if isa(dep,/array) then a = dep[i+2*maxm] else a = dep
      if isa(len,/array) then b = len[i+2*maxm] else b = len
      if b lt 0 then b = strlen(str[i+2*maxm])
      store[i+2*maxm] = strmid(str[i+2*maxm],a,b, REVERSE_OFFSET=reverse_offset)
    endfor
  endif

  return, store
end