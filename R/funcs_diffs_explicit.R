#' 1D diffusion, boundary constant, FD explicit
#'
#' @description Function to solve diffusion equation
#' using finite difference method, forward euler or
#' explicit scheme, applied in the condition where:
#' a) the diffusion occurs inside a slab of material unidimensionally,
#' or only along one axis, b) the coefficient of diffusion is constant,
#' and c) the boundary conditions are constant.
#'
#' @usage mdfexdiffu(D,dt,l,T,C_i=0,C_f=1,F=0.5)
#'
#' @param D Coefficient of diffusion, constant
#' @param dt different between each time steps
#' @param l half-length of the slab, usually in cm
#' @param T Total calculated diffusion time
#' @param C_i Initial concentration value inside the slab
#' @param C_f Dirichlet boundary condition,
#' final concentration coming from the outside of the slab,
#' with one or two elements. If there is only one element,
#' the two sides of the slab (x = -l and x = l) will have
#' the same C_f. If there are
#' two elements, the first element (C_f[1]) will be on
#' x = -l while the second one will be on x = l.
#' @param F Fourier's mesh number, should be less than or equal to 0.5
#' to make sure that the solution is stable
#'
#' @return A matrix with {round(T/dt,0)} number of row and
#' {round(L/(sqrt(D*dt/F)),0)} number of column,
#' profiling the diffusion on slab.
#'
#' @examples
#' C_i = 0.00       # Initial concentration inside the slab
#' C_f = 1.00			  # Final concentration coming from outside
#' D = 10^-7				# Coefficient of diffusion, cm^2/s
#' dt = 60			  	# difference between each time step
#' l = 0.25				  # half-thickness of the slab, in cm
#' F = 0.5			    # Fourier's mesh number
#' T = 432000   		# Total measured time in seconds (~5 days)
#' matC <- mdfexdiffu(D,dt,l,T,C_i,C_f,F)
#'
#' # Using plotly for plotting a contour plot
#' library(plotly)
#'
#' df.list <- list(x = seq(-l,l,length.out = ncol(matC)),
#'                 y = seq(0,T,length.out = nrow(matC)),
#'                 z = matC)
#'
#' plot_ly() %>%
#'     add_contour(x = df.list$y, y = df.list$x, z = matC) %>%
#'     layout(title = "Contour plot diffusion",
#'            xaxis = list(title = "x"),
#'            yaxis = list(title = "t (s)"))

mdfexdiffu <- function(D,dt,l,T,C_i=0,C_f=1,F=0.5){

  if(length(C_f) > 2){
    stop("ERROR: C_f must consist of maximum 2 elements")
  }

  L = 2*l
  Nt = round(T/dt,0)
  # Mesh point in time
  t = seq(0, Nt*dt, length.out=Nt)
  dx = sqrt(D*dt/F)
  Nx = round(L/dx,0)
  # Mesh point in space
  x = seq(0,L,length.out=Nx)
  # Verify
  dx = x[2]-x[1]
  dt = t[2]-t[1]

  # Matrice vide pour u et u_n
  u = matrix(0,ncol=Nx,nrow=Nt)
  u_n	= seq(0,0,length.out=Nx)

  # Condition initiale
  for(i in 1:(Nx+1)){
    u_n[i] = C_i
  }

  # Compute
  for(n in 1:Nt){
    # Compute u at mesh points
    for(i in 2:Nx){
      u[n,i] = u_n[i]+F*(u_n[i-1]-2*u_n[i]+u_n[i+1])
    }
    # Insert boundary conditions
    u[n,1] = C_f[1] ; u[n,Nx] = C_f[length(C_f)]

    # Update u_n before next step
    u_n = u[n,]
  }

  return(u)

}
