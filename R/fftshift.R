#' one line helper function to do quadrant flip for 2D Fourier transforms
#' shifts the low spatial frequencies to the centre of the spectrum
#' this is included in most programming languages, but is absent from the base R Fourier transform tools
#' @export
fftshift <- function(im) {im * (-1)^(row(im) + col(im))}
