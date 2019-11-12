#' 1D diffusion, boundary constant, analytical method
#'
#' @description Function to solve diffusion equation
#' using analytical method, applied with the following conditions:
#' a) the diffusion occurs inside a slab of material unidimensionally,
#' or only along one axis, b) the coefficient of diffusion is constant,
#' and c) the boundary conditions are constant.
#'
#' @usage crankdiffu(C_i,C_f,D,l,xreq,treq,N=2)
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
#' @return the concentration for said x and t points.
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
#' u <- crankdiffu(C_i,C_f,D,l,xreq,treq,N)
#'
#' #2nd Example: comparison with finite difference, explicit method ----
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
#' N_crank <- c(2,5,10,25)
#' u_crank <- matrix(0,
#'                  nrow = length(N_crank),
#'                  ncol = length(xreq))

#'## Analytical (Crank) method
#'for (i in 1:length(N_crank)) {
#'  for (j in 1:length(xreq)) {
#'    u_crank[i,j] <- crankdiffu(C_i,C_f,D,l,x=xreq[j],
#'                               t=treq,N=N_crank[i])
#'  }
#'}
#'
#'## Finite difference, explicit
#'u_diffex <- as.vector(mdfexdiffxtreq(D,dt,l,xreq,treq,T,C_i,C_f,F))
#'
#'## Plot
#'par(mfrow=c(2,2))
#'for (i in 1:nrow(u_crank)) {
#'  plot(xreq,u_diffex,
#'       xlab=bquote(italic(x)),
#'       ylab=bquote((c-c[i])/(c[f]-c[i])),
#'       main = bquote(t~"="~.(treq)~s*","~N~"="~.(N_crank[i])),
#'       ylim = c(min(u_crank),1),
#'       type = "l")
#'  points(xreq,u_crank[i,],pch = i)
#'}

crankdiffu <- function(C_i,C_f,D,l,xreq,treq,N=2){

  x <- xreq
  t <- treq
  # Prepare empty matrix
  u <- matrix(0, ncol = 1, nrow = N)

  # Calculation
  for (i in 1:N) {
    n <- i-1
    u[i] <- (-1)^n/(2*n+1)*exp(-D*(2*n+1)^2*pi^2*t/(4*l^2))*cos((2*n+1)*pi*x/(2*l))
  }

  E <- 1-4/pi*sum(u)

  C_xt <- E*(C_f-C_i)+C_i

  return(C_xt)

}

