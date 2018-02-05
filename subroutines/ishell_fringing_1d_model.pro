Function ishell_fringing_1D_model, x, A
  ; Fits a 1D fringing solution for iShell
  ; A[0] : Amplitude
  ; A[1] : Period
  ; A[2] : Phase
  ; A[4] : Period Slope
  
  model = sin(x*2d0*!dpi/(A[1]+A[3]*x)+A[2])*A[0]+1.
  
  return, model
End