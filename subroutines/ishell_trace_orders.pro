Function ishell_trace_orders, flat_input, DEBUG=debug, ORDERS_STRUCTURE=orders_structure, $
  N_ORDERS=n_orders, MIN_ORDER_SPACING=min_order_spacing, UNEQUAL_ORDER_EDGES=unequal_order_edges
  ;Takes a flat image as input and traces order positions
  ;IDEA (for future): instead of using a central median to trace the initial order positions, use a Z-shaped line to catch all orders at once
  ;Code version history
  ; Version 1.0: First stable version (J. Gagne), March 3, 2017
  ; Version 1.1: Added N_orders keyword for compatibility with K band, January 26, 2018
  
  forward_function weighted_median,robust_sigma,horizontal_median,rvb_hex,remove,poly_fit
  
  ;/// PARAMETERS ///
  ;Total number of orders - default is 29, which is appropriate for K band w/ the gas cell
  if ~keyword_set(n_orders) then $
    n_orders = 29L
  
  if ~keyword_set(unequal_order_edges) then $
    unequal_order_edges = 0L
  
  ;Minimum pixel space between orders (default is 30 pixels for KS band)
  if ~keyword_set(min_order_spacing) then $
    min_order_spacing = 15L
  
  ;Median smooth box size used to detect bad pixels
  bpix_mask_nsmooth = 3L
  
  ;Nsigma used to detect bad pixels
  bpix_mask_nsigma = 1d0
  
;  ;Kernel size for the convolution that will bring out the order edges
;  detection_kernel_size = 13d0
  
  ;Vertical shift for the edge detection
  nvshift = 3L
  
  ;Horizontal smooth to be applied on the order edge detection map
  nhsmooth = 25L 
  
  ;Number of rows/columns to be masked at the edges of the flat field
  nmask_bottom_rows = 4L
  nmask_top_rows = 24L
  nmask_right_cols = 4L
  nmask_left_cols = 4L
  
  ;The height of the running window that traces edges
  ; It must contain enough data points to fit a gaussian,
  ; but not too wide so that it includes more than a single edge 
  running_window_height = 10L
  
  ;Number of data points to skip at both end of each edge for polynomial fitting
  nskip_edge_borders = 3L
  
  ;Degree of polynomial fitting
  ndegree_fit = 3L
  
  ;/// END OF: PARAMETERS ///
  
  ;Create a normalized flat image
  flat = flat_input/weighted_median(flat_input,medval=.9)
  
  ;Do a horizontal normalization of the flat image
  nx = (size(flat_input))[1]
  ny = (size(flat_input))[2]
  column_med_flux = dblarr(nx)
  for i=0L, nx-1L do $
    column_med_flux[i] = weighted_median(reform(flat[i,*]),medval=.96)
  
  column_med_flux_smooth = median(column_med_flux,nhsmooth)>0.05
  column_med_flux_smooth_2d = column_med_flux_smooth#make_array(ny,value=1d0,/double)

  ;Correct the flux
  flat /= column_med_flux_smooth_2d
  
  ;Identify order data and calculate standard deviation
  g_orders_data = where(flat ge 0.1, ng_orders_data)
  sig = robust_sigma(flat[g_orders_data])
  
  ;Create a temporary bad pixel mask
  dev = flat-median(flat,bpix_mask_nsmooth)
  bad = where(abs(dev) ge (sig*bpix_mask_nsigma), nbad)
  flat_masked = flat
  if nbad ne 0L then flat_masked[bad] = !values.d_nan
  
  ;Use median smoothing to remove bad pixels
  bad = where(~finite(flat_masked), nbad)
  if nbad ne 0L then flat_masked[bad] = (median(flat_masked,bpix_mask_nsmooth))[bad]
  
  ;Remove any remaining bad pixels with a second-pass
  bad = where(~finite(flat_masked), nbad)
  if nbad ne 0L then flat_masked[bad] = (median(flat_masked,bpix_mask_nsmooth))[bad]
  
  ;Remove any other remaining bad pixels with a final median
  bad = where(~finite(flat_masked), nbad)
  if nbad ne 0L then flat_masked[bad] = median(flat[g_orders_data])
  
  ;Bring out the order edges
  detection_image = flat_masked - shift(flat_masked,0,-nvshift)
  
;  ;Create a kernel for the convolution that will bring out order edges
;  detection_kernel = [-(findgen(detection_kernel_size)+1), reverse(findgen(detection_kernel_size)+1)] / (2.*detection_kernel_size) 
;  
;  ;Do the convolution
;  detection_image = convol(flat_masked, detection_kernel, /center, /edge_truncate)
  
  ;Smooth the detection image horizontally
  detection_image_smooth = horizontal_median(detection_image,nhsmooth)
  
  ;Mask the flat edges
  if keyword_set(nmask_bottom_rows) then $
    detection_image_smooth[*,0:nmask_bottom_rows-1L] = !values.d_nan
  if keyword_set(nmask_top_rows) then $
    detection_image_smooth[*,-nmask_top_rows:*] = !values.d_nan
  if keyword_set(nmask_left_cols) then $
    detection_image_smooth[0:nmask_left_cols-1L,*] = !values.d_nan
  if keyword_set(nmask_right_cols) then $
    detection_image_smooth[-nmask_right_cols:*,*] = !values.d_nan
  
  ;Pick the detection pattern at the center of the image
  detection_image_center = median(detection_image_smooth[(nx/2L-nhsmooth/2L):(nx/2L+nhsmooth/2L),*],dim=1)
  
  ;Detect order edges in the center of the image first
  ;Loop on orders to select the most solid ones first
  left_positions = lonarr(n_orders+unequal_order_edges)-1L
  left_detection_image_center = -detection_image_center
  
  ;If number of left edges = right edges, loop goes over 1 extra order
  for i=0L, n_orders-1L+unequal_order_edges do begin & $
    maxval = max(left_detection_image_center,/nan,wmax) & $
    left_positions[i] = wmax & $
    ;Mask pixels around the detected max
    left_detection_image_center[(wmax-min_order_spacing)>0L:$
      (wmax+min_order_spacing)<(ny-1L)] = !values.d_nan & $
  endfor
  
  ;Sort the left positions
  left_positions = left_positions[sort(left_positions)]
  
  ;If the number of order edges is not odd, add an artificial left edge at the end
  if unequal_order_edges eq 0 then $
    left_positions = [left_positions,nx-1]
  
  ;Find the right positions by taking the minimum value between each left position
  right_positions = lonarr(n_orders)-1L
  right_detection_image_center = detection_image_center
  for i=0L, n_orders-1L do begin & $
    maxval = max(detection_image_center[left_positions[i]:left_positions[i+1L]],/nan,wmax) & $
    right_positions[i] = wmax+left_positions[i] & $
  endfor
  
  ;If the number of order edges is not odd, remove the artifical left edge at the end
  if unequal_order_edges eq 0 then $
    remove, n_elements(left_positions)-1, left_positions
  
  ;Make some debugging plots
  if keyword_set(debug) then begin
    plot,lindgen(ny),detection_image_center
    oplot,left_positions,interpol(detection_image_center,lindgen(ny),left_positions),col=255,/ps
    oplot,right_positions,interpol(detection_image_center,lindgen(ny),right_positions),col=rvb_hex(100,255,100),/ps
    stop
  endif
  
  ;Drop edges in excess
  if n_elements(left_positions) gt n_elements(right_positions) then $
    remove,n_elements(left_positions)-1,left_positions
  if n_elements(right_positions) gt n_elements(left_positions) then $
    remove,0L,left_positions
  
  ;Move along the edges horizontally to collect maximum positions and trace them
  negative_detection_image_smooth = -detection_image_smooth
  order_left_y_positions = dindgen(nx,n_orders)+!values.d_nan
  order_left_y_positions[nx/2L,*] = double(left_positions)
  for i=0L, n_orders-1L do begin
    last_position = double(left_positions[i])
    for j=nx/2L+1L,nx-1L do begin
      y0 = round((last_position-double(running_window_height)/2d0))>0L
      y1 = round((last_position+double(running_window_height)/2d0))<(ny-1L)
      if ((ny-1L)-y0) lt double(running_window_height)/2d0 then continue
      if y1 lt double(running_window_height)/2d0 then continue
      vertical_slice = negative_detection_image_smooth[j,y0:y1]
      if total(~finite(vertical_slice)) gt 0 then continue
      gauss_par = !NULL
      estimates = [max(vertical_slice,/nan),double(running_window_height)/2d0,double(nvshift)/2d0*0.8553]
      yfit = gaussfit2(dindgen(n_elements(vertical_slice)), vertical_slice, gauss_par,NTERMS=3,ESTIMATES=estimates)
      order_left_y_positions[j,i] = gauss_par[1]+y0
      last_position = order_left_y_positions[j,i]
    endfor
    last_position = double(left_positions[i])
    for j=nx/2L-1L,0L,-1L do begin
      y0 = round((last_position-double(running_window_height)/2d0))>0L
      y1 = round((last_position+double(running_window_height)/2d0))<(ny-1L)
      if ((ny-1L)-y0) lt double(running_window_height)/2d0 then continue
      if y1 lt double(running_window_height)/2d0 then continue
      vertical_slice = negative_detection_image_smooth[j,y0:y1]
      if total(~finite(vertical_slice)) gt 0 then continue
      gauss_par = !NULL
      estimates = [max(vertical_slice,/nan),double(running_window_height)/2d0,double(nvshift)/2d0*0.8553]
      yfit = gaussfit2(dindgen(n_elements(vertical_slice)), vertical_slice, gauss_par,NTERMS=3,ESTIMATES=estimates)
      order_left_y_positions[j,i] = gauss_par[1]+y0
      last_position = order_left_y_positions[j,i]
    endfor
  endfor
  
  ;Repeat for right edges
  order_right_y_positions = dindgen(nx,n_orders)+!values.d_nan
  order_right_y_positions[nx/2L,*] = double(right_positions)
  for i=0L, n_orders-1L do begin
    last_position = double(right_positions[i])
    for j=nx/2L+1L,nx-1L do begin
      y0 = round((last_position-double(running_window_height)/2d0))>0L
      y1 = round((last_position+double(running_window_height)/2d0))<(ny-1L)
      if ((ny-1L)-y0) lt double(running_window_height)/2d0 then continue
      if y1 lt double(running_window_height)/2d0 then continue
      vertical_slice = detection_image_smooth[j,y0:y1]
      if total(~finite(vertical_slice)) gt 0 then continue
      gauss_par = !NULL
      estimates = [max(vertical_slice,/nan),double(running_window_height)/2d0,double(nvshift)/2d0*0.8553]
      yfit = gaussfit2(dindgen(n_elements(vertical_slice)), vertical_slice, gauss_par,NTERMS=3,ESTIMATES=estimates)
      order_right_y_positions[j,i] = gauss_par[1]+y0
      last_position = order_right_y_positions[j,i]
    endfor
    last_position = double(right_positions[i])
    for j=nx/2L-1L,0L,-1L do begin
      y0 = round((last_position-double(running_window_height)/2d0))>0L
      y1 = round((last_position+double(running_window_height)/2d0))<(ny-1L)
      if ((ny-1L)-y0) lt double(running_window_height)/2d0 then continue
      if y1 lt double(running_window_height)/2d0 then continue
      vertical_slice = detection_image_smooth[j,y0:y1]
      if total(~finite(vertical_slice)) gt 0 then continue
      gauss_par = !NULL
      estimates = [max(vertical_slice,/nan),double(running_window_height)/2d0,double(nvshift)/2d0*0.8553]
      yfit = gaussfit2(dindgen(n_elements(vertical_slice)), vertical_slice, gauss_par,NTERMS=3,ESTIMATES=estimates)
      order_right_y_positions[j,i] = gauss_par[1]+y0
      last_position = order_right_y_positions[j,i]
    endfor
  endfor
  
  ;Shift the right edges depending on the vertical shift used to create the edge detection image
  order_right_y_positions += double(nvshift)
  
  if keyword_set(nskip_edge_borders) then begin
    for i=0L, n_orders-1L do begin
      gfin = where(finite(order_left_y_positions[*,i]), ngfin)
      if ngfin ge nskip_edge_borders*2L then begin
        order_left_y_positions[gfin[0L]:(gfin[0L]+nskip_edge_borders-1L),*] = !values.d_nan
        order_left_y_positions[gfin[-1L]-nskip_edge_borders-1L:gfin[-1L],*] = !values.d_nan
      endif
      gfin = where(finite(order_right_y_positions[*,i]), ngfin)
      if ngfin ge nskip_edge_borders*2L then begin
        order_right_y_positions[gfin[0L]:(gfin[0L]+nskip_edge_borders-1L),*] = !values.d_nan
        order_right_y_positions[gfin[-1L]-nskip_edge_borders-1L:gfin[-1L],*] = !values.d_nan
      endif
    endfor
  endif
  
  ;Fit polynomial orders through the edges
  xarr = dindgen(nx)
  left_coeffs = dblarr(n_orders,ndegree_fit)+!values.d_nan
  left_chi2s = dblarr(n_orders)
  left_medianchi2s = dblarr(n_orders)
  right_coeffs = dblarr(n_orders,ndegree_fit)+!values.d_nan
  right_chi2s = dblarr(n_orders)
  right_medianchi2s = dblarr(n_orders)
  mid_coeffs = dblarr(n_orders,ndegree_fit)+!values.d_nan
  mid_chi2s = dblarr(n_orders)
  mid_medianchi2s = dblarr(n_orders)
  for i=0L, n_orders-1L do begin
    ;Fit the left edges with a polynomial
    gfin = where(finite(order_left_y_positions[*,i]), ngfin)
    coeffs = reform(poly_fit(xarr[gfin],order_left_y_positions[gfin,i],ndegree_fit-1L,status=status))
    left_coeffs[i,*] = coeffs
    left_chi2s[i] = total((poly(xarr[gfin],coeffs)-order_left_y_positions[gfin,i])^2,/nan)/double(ngfin)
    left_medianchi2s[i] = median((poly(xarr[gfin],coeffs)-order_left_y_positions[gfin,i])^2)
    
    ;Fit the right edges with a polynomial
    gfin = where(finite(order_right_y_positions[*,i]), ngfin)
    coeffs = reform(poly_fit(xarr[gfin],order_right_y_positions[gfin,i],ndegree_fit-1L,status=status))
    right_coeffs[i,*] = coeffs
    right_chi2s[i] = total((poly(xarr[gfin],coeffs)-order_right_y_positions[gfin,i])^2,/nan)/double(ngfin)
    right_medianchi2s[i] = median((poly(xarr[gfin],coeffs)-order_right_y_positions[gfin,i])^2)
    gfin = where(finite(order_left_y_positions[*,i]) and finite(order_right_y_positions[*,i]), ngfin)
    
    ;Fit the middle position of edges with a polynomial
    order_mid_pos = (order_left_y_positions[gfin,i]+order_right_y_positions[gfin,i])/2d0
    coeffs = reform(poly_fit(xarr[gfin],order_mid_pos,ndegree_fit-1L,status=status))
    mid_coeffs[i,*] = coeffs
    mid_chi2s[i] = total((poly(xarr[gfin],coeffs)-order_mid_pos)^2,/nan)/double(ngfin)
    mid_medianchi2s[i] = median((poly(xarr[gfin],coeffs)-order_mid_pos)^2)
    
    if keyword_set(debug) then begin
      alldata = [order_left_y_positions[gfin,i],order_right_y_positions[gfin,i]]
      yrange = [min(alldata,/nan),max(alldata,/nan)]
      yrange += [-1,1]*(yrange[1]-yrange[0])*.05
      plot, xarr[gfin],order_left_y_positions[gfin,i],xtitle='Pixel position',ytitle='Edge positions',title='Edges for Order #'+strtrim(i+1,2)+'/'+strtrim(n_orders,2)+' (Left, Right and Middle)',thick=2,linestyle=2,yrange=yrange,/ystyle,charsize=1.3
      oplot, xarr[gfin],order_right_y_positions[gfin,i],thick=2,linestyle=2
      oplot, xarr[gfin],order_mid_pos,thick=2,linestyle=0
      
      oplot, xarr, poly(xarr, left_coeffs[i,*]), col=255, linestyle=0
      oplot, xarr, poly(xarr, right_coeffs[i,*]), col=255, linestyle=0
      oplot, xarr, poly(xarr, mid_coeffs[i,*]), col=255, linestyle=0
      stop
    endif
  endfor
  
  ;Create an order mask image
  order_mask = fltarr(nx,ny)+!values.f_nan
  yarr = dindgen(ny)
  x2d = xarr#make_array(ny,value=1d0,/double)
  y2d = make_array(nx,value=1d0,/double)#yarr
  for i=0L, n_orders-1L do begin
    left_2d = poly(x2d,reform(left_coeffs[i,*]))
    right_2d = poly(x2d,reform(right_coeffs[i,*]))
    g_order = where(y2d ge left_2d and y2d le right_2d, ng_order)
    if ng_order ne 0L then $
      order_mask[g_order] = float(i)
  endfor
  
  structure_out = {order_id:-1L,height:!values.f_nan,mid_coeffs:dblarr(ndegree_fit),mid_chi2:!values.d_nan,mid_medianchi2:!values.d_nan,$
    left_coeffs:dblarr(ndegree_fit),left_chi2:!values.d_nan,left_medianchi2:!values.d_nan,order_left_y_positions:dblarr(nx)+!values.d_nan,$
    right_coeffs:dblarr(ndegree_fit),right_chi2:!values.d_nan,right_medianchi2:!values.d_nan,order_right_y_positions:dblarr(nx)+!values.d_nan}
  structure_out = replicate(structure_out,n_orders)
  
  order_heights = median(order_right_y_positions-order_left_y_positions,dim=1)
  for i=0L, n_orders-1L do begin
    structure_out[i].order_id = i
    structure_out[i].height = order_heights[i]
    structure_out[i].mid_coeffs = reform(mid_coeffs[i,*])
    structure_out[i].mid_chi2 = mid_chi2s[i]
    structure_out[i].mid_medianchi2 = mid_medianchi2s[i]
    structure_out[i].left_coeffs = reform(left_coeffs[i,*])
    structure_out[i].left_chi2 = left_chi2s[i]
    structure_out[i].left_medianchi2 = left_medianchi2s[i]
    structure_out[i].right_coeffs = reform(right_coeffs[i,*])
    structure_out[i].right_chi2 = right_chi2s[i]
    structure_out[i].right_medianchi2 = right_medianchi2s[i]
    structure_out[i].order_left_y_positions = order_left_y_positions[*,i]
    structure_out[i].order_right_y_positions = order_right_y_positions[*,i]
  endfor
  
  orders_structure = structure_out
  return, order_mask
  
;  stop
;  plot,order_left_y_positions[*,0],yrange=[0,2048],xrange=[0,2048],/xsty,/ysty
;  for i=1L, n_orders-1L do oplot,order_left_y_positions[*,i]
;  for i=0L, n_orders-1L do oplot,order_right_y_positions[*,i],col=255
;  for i=0L, n_orders-1L do oplot,poly(x2d,reform(mid_coeffs[i,*])),col=rvb_hex(100,255,100)
;  stop
End
