#' one line helper function to do quadrant flip
#' @export
fftshift <- function(im) {im * (-1)^(row(im) + col(im))}
