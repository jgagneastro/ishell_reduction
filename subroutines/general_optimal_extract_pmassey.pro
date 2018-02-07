Function general_optimal_extract_PMassey, trace, sky, profile, readnoise, gain, NITER=niter, ESP=esp, NOSKYSUBTRACTION=noskysubtraction

  forward_function robust_sigma
  
  ;Determine the stddev in the sky (in photons)
  err_sky = robust_sigma(((sky-shift(sky,1,0))[1:*,1:*]+(sky-shift(sky,0,1))[1:*,1:*])/2.)*gain
  
  ;Loop on spatial pixels
  npix = (size(trace))[1]
  xi = lindgen((size(trace))[2])
  spext = findgen(npix)+!values.f_nan
  espext = findgen(npix)+!values.f_nan
  for i=0L, npix-1L do begin
    
    ;Take a slice of the extraction profile
    profilei = reform(profile[i,*])
    
    ;Normalize profile by its maximum
    ; (this ensures that flux * profile re-generates the data exactly)
    profilei /= max(profilei,/nan)
    
    ;Take a slice of the data (in photon units)
    tracei = reform(trace[i,*])*gain
    
    ;Take a slice of the sky (in photon units)
    skyi = reform(sky[i,*])*gain
    
    if keyword_set(noskysubtraction) then $
      skyi *= 0d0
    
    ;Determine the stellar flux (trace-sky)
    if keyword_set(noskysubtraction) then $
      datai = tracei else $
      datai = tracei-skyi
    
    ;Skip if all data is negative
    if total(datai,/nan) eq 0 then begin
      spext[i] = !values.d_nan
      espext[i] = !values.d_nan
      continue
    endif
    
    ;Prepare masking negative data
    wneg_data = where(datai lt 0., nneg_data)
    
    ;Determine variance
    vari = readnoise^2+datai+skyi+err_sky^2
    
    ;Determine extraction weights based on the variance
    weightsi = profilei^2/vari
    
    ;Mask negative data
    if nneg_data ne 0L then $
      weightsi[wneg_data] = 0.
    
    bad = where(~finite(datai*weightsi),nbad)
    if nbad ne 0L then begin
      datai[bad] = !values.f_nan
      weightsi[bad] = 0d0
    endif
    
    ;Normalize weights
    weightsi /= total(weightsi,/nan)
    
    ;Prepare to normalize the data in a way that "flux * profile"
    ; re-generates exactly the observed data
    norm = total(weightsi*profilei,/nan)
    
    ;Extract the data
    spexti = total(datai*weightsi,/nan)/gain/norm
    espexti = sqrt(total(vari,/nan))/gain/norm
    
    ;Store the data
    spext[i] = spexti
    espext[i] = espexti
    
  endfor
  
  ESP = espext
  return, spext
  
End