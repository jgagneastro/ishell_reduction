Function rtrim, xx, i, SIGNIFICATIF = significatif, VIRGULE=virgule
  ;  Cette fonction arrondit le chiffre x a i decimales,
  ;  puis retourne ceci en string.
  ;  SIGNIFICATIF = x : Pour retourner un chiffre a x chiffres significatifs
  
  compile_opt hidden ;Evite de toujours afficher 'Compiled module'
  forward_function round, strepjo, sign, kill_value, string_replace
  on_error, 2
  
  ;Evite d'alterer l'entree
  x = xx
  
  if n_params() eq 1 then i = 0
  
  ;Chiffres significatifs
  if keyword_set(significatif) then begin
    if significatif lt 0 then $
      message, 'ROUNDEFFSTRING : Le keyword SIGNIFICATIF doit etre positif !'
    ordre = -1*floor(alog10(abs(kill_value(xx,0,1))))
    i = ordre + (significatif-1)
  endif
  
  ;Version vecteur
  if isa(x,/array) then begin
    xout = strarr(n_elements(x))
    for j=0, n_elements(x)-1 do $
      xout[j] = isa(i,/array) ? rtrim(x[j],i[j], SIGNIFICATIF = significatif, VIRGULE = virgule) : rtrim(x[j],i, SIGNIFICATIF = significatif, VIRGULE = virgule)
    return, xout
  endif
  
  ;Pour les NaN
  if ~finite(x) then return, 'NaN'
  
  ;On verifie pour les puissances
  stn = strtrim(x,2)
  if strpos(stn,'e') ne -1 then begin
    power = fix(strmid(stn,strlen(stn)-3))
    if power lt i-1 then begin
      return, '0.'+strjoin((strarr(i)+'0'))
    endif
    ;message, 'ROUNDEFFSTRING : N''est pas prevu pour tenir en compte des nombres en notation scientifique.'
  endif
  
  ;On ajuste le string de sortie
  ssgn = sign(x)
  x = abs(x)
  s = round(x*10.d0^(i))/10.d0^(i)
  nstring = strlen(strtrim(fix(s),2))
  if i eq 0 then add = 0 else add = 1
  addzer = (-(i+1) > 0)
  ss = strtrim(s,2)
  rstring = ss
  number = nstring+add+i
  for j=0, n_elements(s)-1 do rstring[j] = strmid(ss[j],0,number[j])
  ;for j=0, n_elements(s)-1 do rstring(j) = strmid(ss(j),0,number(j))
  rstring += strepjo('0',addzer)
  
  if keyword_set(virgule) then $
    string_replace, rstring, '.', ','
    
  ;On renvoie le string
  return, ssgn eq 1 ? rstring : '-'+rstring
End