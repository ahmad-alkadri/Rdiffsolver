#' 1D diffusion, boundary constant, analytical method
#'
#' @description Function to solve diffusion equation
#' using analytical method, applied with the following conditions:
#' a) the diffusion occurs inside a slab of material unidimensionally,
#' or only along one axis, b) the coefficient of diffusion is constant,
#' and c) the boundary conditions are constant.
#'
#' @usage analyticdiffu(C_i,C_f,D,l,xreq,treq,N=2)
#'
#' @param D coefficient of diffusion, constant
#' @param l half-length of the slab, usually in cm
#' @param xreq requested x point
#' @param treq requested t point
#' @param C_i initial concentration value inside the slab
#' @param C_f dirichlet boundary condition,
#' final concentration coming from the outside of the slab,
#' with one or two elements. If there is only one element,
#' the two sides of the slab (x = -l and x = l) will have
#' the same C_f. If there are
#' two elements, the first element (C_f[1]) will be on
#' x = -l while the second one will be on x = l.
#' @param N the number of series taken from the analytical
#' solution (see Crank(1975)). Summation of the series
#' will be conducted from n=0 to n=N-1.
#'
#' @return a matrix with {length(xreq)} number of columns and
#' {length(treq)} number of rows.
#'
#' @examples
#'
#' #1st Example: single value for a single x and a single t points ----
#' C_i = 0.00 # Initial concentration inside the slab
#' C_f = 1.00	# Final concentration coming from outside
#' D = 10^-7 # Coefficient of diffusion, cm^2/s
#' l = 0.25 # half-thickness of the slab, in cm
#' treq = 1e5 # time in s
#' xreq = 0 # point space
#' N = 2 # number of series, thus sum from n=0 to n=N-1
#' u <- analyticdiffu(C_i,C_f,D,l,xreq,treq,N)
#'
#' #2nd Example: compared with finite difference, explicit method ----
#'
#' C_i = 0.00 # Initial concentration inside the slab
#' C_f = 1.00	# Final concentration coming from outside
#' D = 10^-7 # Coefficient of diffusion, cm^2/s
#' l = 0.25 # half-thickness of the slab, in cm
#' dt = 60 # timestep
#' treq = 1e3 # time in s
#' xreq = seq(-0.25,0.25,length.out = 100)
#' F = 0.5 # Fourier's mesh number
#' T = 432000 # Total calculated time of diffusion in seconds (~5 days)
#'
#' N_analytic <- c(2,5,10,25)
#' u_analytic <- matrix(0,
#'                  nrow = length(N_analytic),
#'                  ncol = length(xreq))

#'## Analytical (Crank) method
#'for (i in 1:length(N_analytic)) {
#'  for (j in 1:length(xreq)) {
#'    u_analytic[i,j] <- analyticdiffu(C_i,C_f,D,l,x=xreq[j],
#'                               t=treq,N=N_analytic[i])
#'  }
#'}
#'
#'## Finite difference, explicit
#'u_diffex <- as.vector(mdfexdiffxtreq(D,dt,l,xreq,treq,T,C_i,C_f,F))
#'
#'## Plot
#'par(mfrow=c(2,2),
#'    oma = c(0, 0, 2, 0))
#'for (i in 1:nrow(u_analytic)) {
#'  plot(xreq,u_diffex,
#'       xlab=bquote(italic(x)~(cm)),
#'       ylab=bquote(italic(c(x*","*t))),
#'       main = bquote(N~"="~.(N_analytic[i])),
#'       ylim = c(min(u_analytic),C_f),
#'       type = "l")
#'  points(xreq,u_analytic[i,])
#'}
#'mtext(bquote(italic(t)~"="~.(treq)~s*","
#'      ~italic(D)~"="~.(D)~(cm^2/s)*","~italic(l)~"="~.(l)~cm),
#'      outer = TRUE,
#'      cex = 1)
#'
#' #3rd Example: contour plotting using plotly ----
#'
#' C_i = 0.00 # Initial concentration inside the slab
#' C_f = 1.00	# Final concentration coming from outside
#' D = 10^-7 # Coefficient of diffusion, cm^2/s
#' l = 0.25 # half-thickness of the slab, in cm
#' N=100
#' treq = seq(0,1e5,length.out = 100)
#' xreq = seq(-0.25,0.25,length.out = 100)
#'
#' ureq <- analyticdiffu(C_i,C_f,D,l,xreq,treq,N)
#' # Using plotly for plotting a contour plot
#' library(plotly)
#'
#' df.list <- list(x = xreq,
#'                 y = treq,
#'                 z = ureq)
#'
#' plot_ly() %>%
#'   add_contour(x = df.list$x, y = df.list$y, z = df.list$z) %>%
#'   layout(title = "Contour plot diffusion",
#'          xaxis = list(title = "Requested x points (cm)"),
#'          yaxis = list(title = "t (s)"))

analyticdiffu <- function(C_i,C_f,D,l,xreq,treq,N=2){

  # Prepare empty matrix with column for X and rows for T
  ureq <- matrix(0,
                 ncol = length(xreq),
                 nrow = length(treq))

  # Calculation
  for(i in 1:length(treq)){

    for(j in 1:length(xreq)){

      x <- xreq[j]
      t <- treq[i]
      u <- matrix(0, ncol = 1, nrow = N)
      for (k in 1:N) {
        n <- k-1
        u[k] <- (-1)^n/(2*n+1)*
          exp(-D*(2*n+1)^2*pi^2*t/(4*l^2))*
          cos((2*n+1)*pi*x/(2*l))
      }
      E <- 1-4/pi*sum(u)
      C_xt <- E*(C_f-C_i)+C_i

      # Put the result in the matrix ureq
      ureq[i,j] <- C_xt
    }

  }

  if(length(xreq) == 1 | length(treq) == 1){

    ureq <- as.vector(ureq)

  }

  return(ureq)

}
