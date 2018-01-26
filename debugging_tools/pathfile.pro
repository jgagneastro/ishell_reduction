Function pathfile, pathin, sep=sep, SILENT=silent, REMAIN=remain, NOEXTENSION=noextension
;  Cette fonction prend en entree un path et retourne seulement le fichier au bout de celui-ci
;  /Sep sert a definir le path_separation
;  /Remain retourne tout le chemin sauf le fichier
;  /Noexp enleve l'extension a la fin du fichier
  
  compile_opt hidden
  forward_function pathextension, strinvert, addslash
  on_error, 2
  
  ;Pour eviter de modifier la variable dans le programme-pere
  if ~isa(pathin,'string') then message, 'Le chemin entre doit etre un string !'
  path = pathin
  
  ;Version array
  if isa(path,/array) then begin
    retour = strarr(n_elements(path))
    for i=0, n_elements(path)-1 do $
      retour[i] = pathfile(path[i], sep=sep, SILENT=silent, REMAIN=remain, NOEXTENSION=noextension)
    return, retour
  endif
  
  if keyword_set(remain) then $
    if pathfile(path) eq path then return, ''
  
  if keyword_set(noextension) then begin
    ext = pathextension(path)
    path = strmid(path,0,strlen(path)-strlen(ext))
  endif
  
  if keyword_set(sep) then pathsep = sep else pathsep = path_sep()
  invpath = strinvert(path)
  dash = strpos(invpath,pathsep)
  if dash eq 0 then begin
    invpath = strmid(invpath,1,strlen(invpath)-1)
    dash = strpos(invpath,pathsep)
  endif
  ;Si le nom donné n'est pas un path, mais juste un fichier
  if dash eq -1 and strtrim(path,2) ne '' then begin & return, path & endif ;if ~keyword_set(silent) then print, 'pathfile : ERREUR !!! : '+strtrim(path,2) & 
  invfile = strmid(invpath,0,dash)
  invpath_remain = strmid(invpath,dash)
  file = strinvert(invfile)
  path_remain = strinvert(invpath_remain)
  if keyword_set(remain) then begin
    addslash, path_remain
    return, path_remain
  endif
  return, file
End