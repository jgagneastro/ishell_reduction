Function fit_2d_curved_trace, x, y, A, _EXTRA=extra, NOBPIXPATCH=nobpixpatch
  
  ;Recreate the y position of the trace as a 2D map
  ytrace_coeff = A[0:extra.npoly_ypos-1L]
  ytrace_pos = poly(x,ytrace_coeff)
  
  ;Recreate the width of the gaussian as a 2D map
  sigma_coeff = A[extra.npoly_ypos:*]
  sigma = poly(x,sigma_coeff)
  
  ;Recreate the gaussian profile (normalized)
  gaussian_profile = exp(-(y-ytrace_pos)^2/(2d0*sigma^2))
  
  ;Put the spectrum on top of it
  model = gaussian_profile*extra.spectrum
  
  ;Mask the 100 worst pixels
  if ~keyword_set(nobpixpatch) then begin
    chi2 = (extra.data-model)^2*extra.weights
    gweights = where(extra.weights ne 0L, ngweights)
    chi2s = chi2[sort(chi2)]
    bad = where(chi2 gt chi2s[-100L], nbad)
    if nbad ne 0L then $
      model[bad] = extra.data[bad]
  endif
  
  ;Return the model
  return, model
  
End