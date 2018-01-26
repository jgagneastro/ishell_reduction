Function pathextension, path, sep=sep, remain=remain, SILENT=silent
;  Cette fonction prend en entree un path et retourne seulement l'extention du fichier au bout de 
;  celui-ci.
;  Sep sert a definir le path_separation
;  Remain : Retourne seulement le nom du fichier, sans extension
  compile_opt hidden
  forward_function strinvert, pathfile
  
  ;Version array
  if isa(path,/array) then begin
    retour = strarr(n_elements(path))
    for i=0, n_elements(path)-1 do $
      retour[i] = pathextension(path[i],sep=sep,remain=remain)
    return, retour
  endif
  
  file = pathfile(path,sep=sep,/silent)
  invfile = strinvert(file)
  
  dot = strpos(invfile,'.')
  if dot eq 0 then begin
    invfile = strmid(invfile,1,strlen(invfile)-1)
    dot = strpos(invfile,'.')
    if keyword_set(remain) then return, ''
  endif
  if dot eq -1 then begin & if ~keyword_set(silent) then print, 'pathextension : ERREUR !!!' & return, '' & endif
  invfile = keyword_set(remain) ? strmid(invfile,dot+1) : strmid(invfile,0,dot+1)
  ext = strinvert(invfile)
  return, ext
End