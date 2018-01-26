Function addzero, numin, nzeros, SPACE=space, TRIM=trim, REVERSE=reverse
;+
; NAME:
;       ADDZERO
;       
; PURPOSE:
;       Adds a number of trailing zeros (or spaces) before a number while transforming it to a string.
;       
; CALLING SEQUENCE:
;       string = ADDZERO( numin, nzeros[, /SPACE] )
;
; INPUTS:
;       NUMIN = Scalar or Array of integer numbers.
;       NZEROS = Scalar indicating the number of trailing zeros to add.
;       TRIM = Trims trailing zeros at the end of decimal numbers (see trim() function)
;
; OPTIONAL INPUT KEYWORD:
;       /SPACE - If this keyword is set then trailing spaces are added instead of trailing zeros.
;
; OUTPUTS:
;       STRING - A string consisting of a number with some trailing zeros.
;
; NOTE:
; 
;       You can avoid using this routine if you only want to directy print a number with some trailing zeros. In
;       this case, use the following command :
;  
;       IDL> print, numin, format='(I0'+strtrim(nzeros,2)+')')
;
; RESTRICTIONS:
;       (1) If you give a floating-point number to this function, it will convert it into an integer.
;
; PROCEDURES USED:
;       SIGN(), STRREP(), FRAC(), TRIM()
;
; MODIFICATION HISTORY:
;       WRITTEN, Jonathan Gagne, September, 3 2011
;       Added support for fractional numbers, Jonathan Gagne, June, 22 2012
;-
  
  compile_opt hidden
  forward_function sign, strrep, frac, trim
  on_error, 2
  
  if ~keyword_set(nzeros) then nzeros = 2
  nzeros >= 0
  
  ;Array version
  if n_elements(numin) ne 1 then begin
    e_fin = strarr(n_elements(numin))
    for i=0, n_elements(numin)-1 do $
      e_fin[i] = addzero(numin[i], nzeros, SPACE = space, TRIM = trim, REVERSE = reverse)
    return, e_fin
  endif
  if keyword_set(space) then add = ' ' else add = '0'
  x = long(numin)
  if keyword_set(trim) then $
    xfrac = trim(abs(numin)) else $
    xfrac = strtrim(abs(numin),2)
  neg = replicate('',n_elements(x))
  wneg = where(sign(ceil(x)) eq -1,nc)
  if nc ne 0 then begin
    neg[wneg] = '-'
  endif
  x = strtrim(abs(x),2)
  if keyword_set(reverse) then $
    x = xfrac+strrep(add,nzeros-strlen(x)) else $
    x = strrep(add,nzeros-strlen(x))+xfrac
  return, (neg+x)[0]
end