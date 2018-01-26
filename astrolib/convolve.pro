function convolve, image, psf, FT_PSF=psf_FT, FT_IMAGE=imFT, NO_FT=noft, $
                        CORRELATE=correlate, AUTO_CORRELATION=auto
;+
; NAME:
;       CONVOLVE
; PURPOSE:
;       Convolution of an image with a Point Spread Function (PSF)
; EXPLANATION:
;       The default is to compute the convolution using a product of 
;       Fourier transforms (for speed).
;
; CALLING SEQUENCE:
;
;       imconv = convolve( image1, psf, FT_PSF = psf_FT )
;  or:
;       correl = convolve( image1, image2, /CORREL )
;  or:
;       correl = convolve( image, /AUTO )
;
; INPUTS:
;       image = 2-D array (matrix) to be convolved with psf
;       psf = the Point Spread Function, (size < or = to size of image).
;
; OPTIONAL INPUT KEYWORDS:
;
;       FT_PSF = passes out/in the Fourier transform of the PSF,
;               (so that it can be re-used the next time function is called).
;       FT_IMAGE = passes out/in the Fourier transform of image.
;
;       /CORRELATE uses the conjugate of the Fourier transform of PSF,
;               to compute the cross-correlation of image and PSF,
;               (equivalent to IDL function convol() with NO rotation of PSF)
;
;       /AUTO_CORR computes the auto-correlation function of image using FFT.
;
;       /NO_FT overrides the use of FFT, using IDL function convol() instead.
;               (then PSF is rotated by 180 degrees to give same result)
; METHOD:
;       When using FFT, PSF is centered & expanded to size of image.
; HISTORY:
;       written, Frank Varosi, NASA/GSFC 1992.
;       Appropriate precision type for result depending on input image
;                               Markus Hundertmark February 2006
;       Fix the bug causing the recomputation of FFT(psf) and/or FFT(image)
;                               Sergey Koposov     December 2006
;-
        compile_opt idl2
        sp = size( psf_FT,/str )  &  sif = size( imFT, /str )
        sim = size( image )  &  sc = sim/2  &  npix = N_elements( image )

        if (sim[0] NE 2) OR keyword_set( noft ) then begin
                if keyword_set( auto ) then begin
                        message,"auto-correlation only for images with FFT",/INF
                        return, image
                  endif else if keyword_set( correlate ) then $
                                return, convol( image, psf ) $
                        else    return, convol( image, rotate( psf, 2 ) )
           endif

        if (sif.N_dimensions NE 2) OR ((sif.type NE 6) AND (sif.type NE 9)) OR $
           (sif.dimensions[0] NE sim[1]) OR (sif.dimensions[1] NE sim[2]) then imFT = FFT( image,-1 )

        if keyword_set( auto ) then $
         return, shift( npix*real_part(FFT( imFT*conj( imFT ),1 )), sc[1],sc[2] )

        if (sp.N_dimensions NE 2) OR ((sp.type NE 6) AND (sp.type NE 9)) OR $
           (sp.dimensions[0] NE sim[1]) OR (sp.dimensions[1] NE sim[2]) then begin
                sp = size( psf )
                if (sp[0] NE 2) then begin
                        message,"must supply PSF matrix (2nd arg.)",/INFO
                        return, image
                   endif
                Loc = ( sc - sp/2 ) > 0         ;center PSF in new array,
                s = (sp/2 - sc) > 0        ;handle all cases: smaller or bigger
                L = (s + sim-1) < (sp-1)
                psf_FT = conj(image)*0 ;initialise with correct size+type according 
                ;to logic of conj and set values to 0 (type of psf_FT is conserved)  
                psf_FT[ Loc[1], Loc[2] ] = psf[ s[1]:L[1], s[2]:L[2] ]
                psf_FT = FFT( psf_FT, -1, /OVERWRITE )
           endif

        if keyword_set( correlate ) then $
                conv = npix * real_part(FFT( imFT * conj( psf_FT ), 1 ))  $
          else  conv = npix * real_part(FFT( imFT * psf_FT, 1 )) 

        sc = sc + (sim MOD 2)   ;shift correction for odd size images.

return, shift( conv, sc[1], sc[2] )
end
