Function horizontal_median, image, width, EVEN=even
  nx = (size(image))[1]
  ny = (size(image))[2]
  output_image = dblarr(nx,ny)+!values.d_nan
  for i=0L, ny-1L do $
    output_image[*,i] = median(image[*,i], width, EVEN=even)
  return, output_image
End