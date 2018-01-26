Function correct_bad_pixels_posterior, flux, NSIG=nsig, NSHIFTS=nshifts
  
  forward_function robust_sigma, create_nan
  
  ;Parameters
  if ~keyword_set(nshifts) then $
    nshifts = [1L,2L,3L,1L]
  if ~keyword_set(nsig) then $
    nsig = [12.,15.,20.,15.]
  nfilters = n_elements(nshifts)
  
  ;Create a flux array that won't modify input
  fl = flux
  
  for i=0L, nfilters-1L do begin
    ;Create an array of deviations
    devleft = (abs(fl-shift(fl,nshifts[i])))
    devright = (abs(fl-shift(fl,-nshifts[i])))
    dev = devleft<devright
    
    ;Fix NaN problems
    bad = where(finite(devleft) and ~finite(devright), nbad)
    if nbad ne 0L then dev[bad] = devleft[bad]
    bad = where(finite(devright) and ~finite(devleft), nbad)
    if nbad ne 0L then dev[bad] = devright[bad]
    
    ;Compute robust stddev
    sig = robust_sigma(dev,/nan)
    
    ;Identify outliers
    bad = where(dev/sig gt nsig[i], nbad)
    
    ;Mask bad pixels
    if nbad ne 0L then fl[bad] = create_nan(size(fl,/type))
;    stop
;    g = where(finite(flux))
;    xrange = [min(g),max(g)]+[-20,20]
;    ;xrange = [600,800]
;    xrange = [700,750]
;    plot,flux,xrange=xrange,/ps,/xsty & oplot, fl, col=255,/ps
  endfor
  return, fl
  
End