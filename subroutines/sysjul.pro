Function sysjul
  ;Donne le moment présent (system time) en jour julien.
  compile_opt hidden ;Evite de toujours afficher 'Compiled module'
  forward_function jdcnv
  months = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
  mon = where(months eq strmid(systime(),4,3))
  if mon[0] eq -1 then return, '-1'
  jdcnv, fix(strmid(systime(),20,4)), mon[0]+1, fix(strmid(systime(),8,2)), float(strmid(systime(),11,2))+float(strmid(systime(),14,2))/60.d0+float(strmid(systime(),17,2))/3600.d0, sysjul
  return, double(sysjul)
End