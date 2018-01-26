Pro printarr, arr1in, arr2in, arr3in, arr4in, arr5in, arr6in, arr7in, arr8in, arr9in, arr10in, arr11in, arr12in, $
              arr13in, arr14in, arr15in, arr16in, arr17in, arr18in, arr19in, arr20in, arr21in, arr22in, arr23in, $
              arr24in, arr25in, arr26in, arr27in, arr28in, arr29in, arr30in, arr31in, arr32in, arr33in, arr34in, $
              arr35in, arr36in, arr37in, arr38in, arr39in, arr40in, arr41in, arr42in, arr43in, arr44in, arr45in, $
              arr46in, arr47in, arr48in, arr49in, arr50in, arr51in, arr52in, arr53in, arr54in, arr55in, arr56in, $
              arr57in, arr58in, arr59in, arr60in, arr61in, arr62in, arr63in, arr64in, arr65in, arr66in, arr67in, $
              arr68in, arr69in, arr70in, arr71in, arr72in, arr73in, arr74in, arr75in, arr76in, arr77in, arr78in, $
              arr79in, arr80in, SPACE=space, JUSTIFY=justify, $
              NOCOMMA=nocomma, SYMBOL=symbol, INDICES=indices, TITLE=title, TITSPACE=titspace, $
              CUTDEC=cutdec, LUN=lun, TEXT=text, COPY=copy, IND2=ind2, REPLACENAN=replacenan
;  Sert a imprimer un array element par element
;  Toujours mettre le plus long array en premier
;  Le keyword Space sert a mettre un espace avant chaque print
;  Il est aussi possible de debuter par exemple par un ! avec SPACE='!'
;  CUTDEC : imprime roundeffstring(arri,cutdec)
  
  ;Evite de toujours afficher 'Compiled module'
  compile_opt hidden
  on_error, 2
  forward_function strepjo, valid_numarr, roundeffstring
  
  ;On determine le plus long array
  npar0 = n_params()
  npar = 0L
  outstring = []
  for i=0, npar0-1 do void = execute('if keyword_set(arr'+strtrim(i+1,2)+'in) then begin & arr'+strtrim(i+1,2)+' = arr'+strtrim(i+1,2)+'in & npar += 1 & endif')
  
  if keyword_set(ind2) then begin
    for i=0, npar-1 do void = execute('if keyword_set(arr'+strtrim(i+1,2)+'in) then begin & arr'+strtrim(i+1,2)+' = arr'+strtrim(i+1,2)+'in[ind2] & endif')
  endif
  
  if ~keyword_set(arr1) then return
  ;if keyword_set(title) then $
  ;  if n_elements(title) ne npar then message, ' Title doit avec la meme longueur que le nombre de vecteurs imprimes !'
  sz = intarr(npar)+keyword_set(title)+keyword_set(titspace)
  for z=1, npar do begin
    R = execute('sz[z-1] += n_elements(arr'+strtrim(z,2)+')')
  endfor
  if ~keyword_set(symbol) then symbol = ' , '
  if keyword_set(title) and keyword_set(indices) then indices = [0, indices+1]
  
  if keyword_set(cutdec) and n_elements(cutdec) eq 1 then cutdec = replicate(cutdec, npar)
  if keyword_set(cutdec) then $
    for i=0, npar-1 do begin
      void = execute('arri = arr'+strtrim(i+1,2))
      if min(valid_numarr(arri)) eq 0 then continue
      arri = strtrim(arri, 2)
      arri = roundeffstring(arri, cutdec[i])
      void = execute('arr'+strtrim(i+1,2)+' = arri')
    endfor
  
  maxlen = 0
  for i=0, npar-1L do begin
    void = execute('arri = arr'+strtrim(i+1,2))
    if keyword_set(titspace) then sp = ','''' ' else sp = ''
    if keyword_set(title) then void = execute('arr'+strtrim(i+1,2)+' = [title[i]'+sp+',strtrim(arr'+strtrim(i+1,2)+',2)]')
    narri = n_elements(arri)+keyword_set(title)+keyword_set(titspace)
    if narri gt maxlen then maxlen = narri
  endfor
  if ~keyword_set(indices) then indices = indgen(maxlen)
  maxsz = max(sz)
  
  text = []
  for ii=0, maxsz-1L do begin
    if ii lt n_elements(indices) then $
      i = indices[ii] else i = ii
    if n_elements(indices) ne 0 and ii ge n_elements(indices) then return ;Ligne ajoutee pour tenter de regler le probleme qui affiche tout en double quand on fournit "indices"
    void = where(indices eq i, ni)
    if ni lt 1 then continue
    if keyword_set(space) then begin
      if strtrim(space,2) eq '1' then str = ' ' else str = space
    endif else str = ''
    for j=1, npar do begin
      virg = ''
      if ~keyword_set(nocomma) then $
        if j eq 1 then virg = '' else virg = symbol
      if keyword_set(justify) then begin
        R = execute('if keyword_set(arr'+strtrim(j,2)+') then len = max(strlen(strtrim(arr'+strtrim(j,2)+',2)))')
        R = execute('if keyword_set(arr'+strtrim(j,2)+') and i lt n_elements(arr'+strtrim(j,2)+') then add = len-strlen(strtrim(arr'+strtrim(j,2)+'[i],2))')
        R = execute('if keyword_set(arr'+strtrim(j,2)+') and i ge n_elements(arr'+strtrim(j,2)+') then add = len-1')
      endif else add = 0
      R = execute('if keyword_set(arr'+strtrim(j,2)+') and i lt n_elements(arr'+strtrim(j,2)+') then str += virg+strtrim(arr'+strtrim(j,2)+'[i],2)+strepjo('' '',add) else if keyword_set(arr'+strtrim(j,2)+') then str+= virg+''-''+strepjo('' '',add)')
    endfor
    if keyword_set(replacenan) then $
      string_replace, str, 'NaN', replacenan
    if keyword_set(lun) then printf, lun, str else $
      print, str
    text = [text, str]
  endfor
  if keyword_set(copy) then copy_clipboard, text
End