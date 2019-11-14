#' 1D diffusion, boundary constant, FD implicit
#'
#' @description Function to solve diffusion equation
#' using finite difference method, backward euler or
#' implicit scheme, applied in the condition where:
#' a) the diffusion occurs inside a slab of material unidimensionally,
#' or only along one axis, b) the coefficient of diffusion is constant,
#' and c) the boundary conditions are constant.
#'
#' @usage mdfimpdiffu(D,Nx,Nt,l,T,C_i=0,C_f=1)
#'
#' @param D Coefficient of diffusion, constant
#' @param Nx number of points in space (x)
#' @param Nt number of points in time (t)
#' @param l half-length of the slab, usually in cm
#' @param T Total calculated diffusion time, usually in seconds
#' @param C_i Initial concentration value inside the slab
#' @param C_f Dirichlet boundary condition,
#' final concentration coming from the outside of the slab,
#' with one or two elements. If there is only one element,
#' the two sides of the slab (x = -l and x = l) will have
#' the same C_f. If there are
#' two elements, the first element (C_f[1]) will be on
#' x = -l while the second one will be on x = l.
#'
#' @return A matrix with Nt number of row and
#' Nx number of column, profiling the diffusion.
#'
#' @examples
#' C_i = 0.00 # Initial concentration inside the slab
#' C_f = 1.00	# Final concentration coming from outside
#' D = 1e-7 # Coefficient of diffusion, cm^2/s
#' Nt = 100 # points in time
#' Nx = 20 # points in space
#' l = 0.25 # half-thickness of the slab, in cm
#' T = 432000 # Total measured time in seconds (~5 days)
#' u <- mdfimpdiffu(D,Nx,Nt,l,T,C_i,C_f)
#'
#' # Using plotly for plotting a contour plot
#' library(plotly)
#'
#' df.list <- list(x = seq(-l,l,length.out = Nx),
#'                 y = seq(0,T,length.out = Nt),
#'                 z = u)
#'
#' plot_ly() %>%
#'     add_contour(x = df.list$x, y = df.list$y, z = df.list$z) %>%
#'     layout(title = "Contour plot diffusion",
#'            xaxis = list(title = "x"),
#'            yaxis = list(title = "t (s)"))

mdfimpdiffu <- function(D,Nx,Nt,l,T,C_i=0,C_f=1){

  if(length(C_f) > 2){
    stop("ERROR: C_f must consist of maximum 2 elements")
  }

  L = 2*l
  # Mesh point in time
  t = seq(0, T, length.out=Nt)
  # Mesh point in space
  x = seq(0,L,length.out=Nx)
  # Step size for space (dx)
  dx = x[2]-x[1]
  # Step size for space (dt)
  dt = t[2]-t[1]
  # F
  F = D*dt/(dx^2)
  
  # Structures for the linear system
  A = matrix(0, nrow = Nx+1, ncol = Nx+1)

  for(i in 2:Nx){
     A[i,i-1] = -F
     A[i,i+1] = -F
     A[i,i] = 1+2*F
  }

  # Boundary condition
  A[1,1] = 1
  A[Nx+1,Nx+1] = 1

  # Empty matrix for u and u_n
  u = matrix(0, nrow = Nx+1, ncol = Nt+1)
  # Boundary condition on the borders
  u[1,] = C_f[1]
  u[Nx+1,] = C_f[length(C_f)]
  # Initial condition
  u[2:Nx,] = C_i

  # Compute
  for(n in 2:Nt+1){
     u_n = u[,n-1]
     b = matrix(0, nrow = Nx+1, ncol = 1)
     b[2:Nx] = u_n[2:Nx]
     b[1] = C_f
     b[Nx+1] = C_f[length(C_f)]
     ui = solve(A,b)
     u[,n] = ui
  }
  
  u = t(u)

  return(u)

}
