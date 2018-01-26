function strepjo, str, n
  ;C'est strep.pro, mais pour eviter le conflit avec celui de John Hopkins
  compile_opt hidden ;Evite de toujours afficher 'Compiled module'
  on_error, 2
  
  if n[0] le 0 and ~isa(n,/array) then return, ''
  
  ;Version array
  if isa(n,/array) then begin
    retour = strarr(n_elements(n))
    for i=0,n_elements(n)-1 do $
      retour[i] = strepjo(str,n[i])
    return, retour
  endif
  
  ;Cette fonction définit un string composé de la séquence str répétée n fois
  return, strjoin(strarr(n)+str)
end