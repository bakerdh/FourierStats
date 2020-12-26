fftshift <- function(im) {im * (-1)^(row(im) + col(im))}  
