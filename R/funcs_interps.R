#' Interpolating the values between two vectors.
#'
#' @description A wrapper of \emph{diff_1D} function,
#' returning the mean values of concentration during each
#' time step of the diffusion.
#'
#' @usage vectInterp(vect1, vect2, dist)
#'
#' @param vect1 First vector variable.
#' @param vect2 Second vector variable.
#' @param dist Desired interpolation value of distance between the two vectors.
#'
#' @return A vector with interpolated values between the first and the second vectors.
#'
#' @examples
#' library(Rdiffsolver)
#' vect1 <- seq(1,100,length.out = 10)
#' vect2 <- seq(2.5,250,length.out = 10)
#' dist <- 0.12
#' vect_int <- vectInterp(vect1,vect2,dist)

vectInterp <- function(vect1, vect2, dist){

  if(length(vect1) < length(vect2)){

    warning("Lengths of the two vectors are not equal")
    return("Error: interpolation failed")

  } else {

    if(length(vect1) > length(vect2)){

      warning("Lengths of the two vectors are not equal")
      return("Error: interpolation failed")

    } else {

      fun_int <- apply(rbind(vect1,vect2), 2, approxfun, x=c(0,1))

      v_int <- vector(length = length(vect1))

      for (i in 1:length(v_int)) {

        v_int[i] <- fun_int[[i]](dist)

      }

      return(v_int)

    }

  }
}
