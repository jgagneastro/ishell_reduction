Function whiten_color, rgb_code, intensity
  if size(rgb_code,/type) eq 7L then begin
    string_replace, rgb_code, ' ', '_'
    void = execute('rgb_code = !color.'+rgb_code[0L])
    if void eq 0 then message, 'An execution statement has failed !'
  endif
  return, byte(round(double(rgb_code)+(1.-double(rgb_code)/255)*double(intensity)/100d0*255d0))
End