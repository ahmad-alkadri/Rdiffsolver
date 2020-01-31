#' 1D diffusion, boundary constant, FD explicit
#' for specific time (t) point(s)
#'
#' @description Function to solve diffusion equation
#' using finite difference method, forward euler or
#' explicit scheme, applied in the condition where:
#' a) the diffusion occurs inside a slab of material unidimensionally,
#' or only along one axis, b) the coefficient of diffusion is constant,
#' and c) the boundary conditions are constant. It will return
#' the solution values for the requested time point(s).
#'
#' @usage mdfexdifftreq(D,dt,l,treq,T,C_i=0,C_f=1,F=0.5)
#'
#' @param D Coefficient of diffusion, constant
#' @param dt different between each time steps
#' @param l half-length of the slab, usually in cm
#' @param treq requested time value(s)
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
#' @return A matrix with {length(treq)} number of row and
#' {round(L/(sqrt(D*dt/F)),0)} number of column,
#' profiling the diffusion on the slab for the requested
#' time point(s).
#'
#' @examples
#' C_i = 0.00 # Initial concentration inside the slab
#' C_f = 1.00	# Final concentration coming from outside
#' D = 10^-7 # Coefficient of diffusion, cm^2/s
#' dt = 60 # difference between each time step
#' l = 0.25 # half-thickness of the slab, in cm
#' F = 0.5 # Fourier's mesh number
#' T = 432000 # Total measured time in seconds (~5 days)
#' treq = seq(1e2,1e5,length.out = 40) # requested time points
#' ureq <- mdfexdifftreq(D,dt,l,treq,T,C_i,C_f,F)
#'
#' # Using plotly for plotting a contour plot
#' library(plotly)
#'
#' df.list <- list(x = seq(-l,l,length.out = ncol(ureq)),
#'                 y = treq,
#'                 z = ureq)
#'
#' plot_ly() %>%
#'     add_contour(x = df.list$x, y = df.list$y, z = df.list$z) %>%
#'     layout(title = "Contour plot diffusion",
#'            xaxis = list(title = "x"),
#'            yaxis = list(title = "Requested time points (s)"))

mdfexdifftreq <- function(D,dt,l,treq,T,C_i=0,C_f=1,F=0.5){

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

  u_req <- matrix(0,nrow=length(treq),ncol = ncol(u))

  for (j in 1:nrow(u_req)) {
    for (i in 1:ncol(u)) {
      u_req[j,i] <- pracma::interp1(t,u[,i],treq[j])
    }
  }

  return(u_req)

}

#' 1D diffusion, boundary constant, FD explicit
#' for specific space (x) point(s)
#'
#' @description Function to solve diffusion equation
#' using finite difference method, forward euler or
#' explicit scheme, applied in the condition where:
#' a) the diffusion occurs inside a slab of material unidimensionally,
#' or only along one axis, b) the coefficient of diffusion is constant,
#' and c) the boundary conditions are constant. It will return
#' the solution values for the requested space (x) point(s).
#'
#' @usage mdfexdiffxreq(D,dt,l,xreq,T,C_i=0,C_f=1,F=0.5)
#'
#' @param D Coefficient of diffusion, constant
#' @param dt different between each time steps
#' @param l half-length of the slab, usually in cm
#' @param xreq requested x point(s)
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
#' {length(xreq)} number of column,
#' profiling the diffusion in the slab for the requested
#' x points.
#'
#' @examples
#' C_i = 0.00 # Initial concentration inside the slab
#' C_f = 1.00	# Final concentration coming from outside
#' D = 10^-7 # Coefficient of diffusion, cm^2/s
#' dt = 60 # difference between each time step
#' l = 0.25 # half-thickness of the slab, in cm
#' F = 0.5 # Fourier's mesh number
#' T = 432000 # Total measured time in seconds (~5 days)
#' xreq = seq(-0.25,0,length.out = 20) # requested space (x) points
#' ureq <- mdfexdiffxreq(D,dt,l,xreq,T,C_i,C_f,F)
#'
#' # Using plotly for plotting a contour plot
#' library(plotly)
#'
#' df.list <- list(x = xreq,
#'                 y = seq(0,T,length.out = nrow(ureq)),
#'                 z = ureq)
#'
#' plot_ly() %>%
#'     add_contour(x = df.list$x, y = df.list$y, z = df.list$z) %>%
#'     layout(title = "Contour plot diffusion",
#'            xaxis = list(title = "Requested x points (cm)"),
#'            yaxis = list(title = "t (s)"))

mdfexdiffxreq <- function(D,dt,l,xreq,T,C_i=0,C_f=1,F=0.5){

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

  u_req <- matrix(0,nrow = nrow(u), ncol = length(xreq))

  x_ins <- seq(-l,l,length.out = ncol(u))

  for(i in 1:ncol(u_req)){
    for(j in 1:nrow(u_req)){
      u_req[j,i] <- pracma::interp1(x_ins,u[j,],xreq[i])
    }
  }

  return(u_req)

}

#' 1D diffusion, boundary constant, FD explicit
#' for specific space (x) and time (t) point(s)
#'
#' @description Function to solve diffusion equation
#' using finite difference method, forward euler or
#' explicit scheme, applied in the condition where:
#' a) the diffusion occurs inside a slab of material unidimensionally,
#' or only along one axis, b) the coefficient of diffusion is constant,
#' and c) the boundary conditions are constant. It will return
#' the solution values for the requested space (x)
#' and time (t) point(s).
#'
#' @usage mdfexdiffxtreq(D,dt,l,xreq,treq,T,C_i=0,C_f=1,F=0.5)
#'
#' @param D Coefficient of diffusion, constant
#' @param dt different between each time steps
#' @param l half-length of the slab, usually in cm
#' @param xreq requested x point(s)
#' @param treq requested t point(s)
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
#' @return A matrix with {length(treq)} number of row and
#' {length(xreq)} number of column,
#' profiling the diffusion in the slab for the requested
#' x and t points.
#'
#' @examples
#' C_i = 0.00 # Initial concentration inside the slab
#' C_f = 1.00	# Final concentration coming from outside
#' D = 10^-7 # Coefficient of diffusion, cm^2/s
#' dt = 60 # difference between each time step
#' l = 0.25 # half-thickness of the slab, in cm
#' F = 0.5 # Fourier's mesh number
#' T = 432000 # Total measured time in seconds (~5 days)
#' xreq = seq(-0.25,0,length.out = 20) # requested space (x) points
#' treq = seq(1e2,1e5,length.out = 40) # requested time (t) points
#' ureq <- mdfexdiffxtreq(D,dt,l,xreq,treq,T,C_i,C_f,F)
#'
#' # Using plotly for plotting a contour plot
#' library(plotly)
#'
#' df.list <- list(x = xreq,
#'                 y = treq,
#'                 z = ureq)
#'
#' plot_ly() %>%
#'     add_contour(x = df.list$x, y = df.list$y, z = df.list$z) %>%
#'     layout(title = "Contour plot diffusion",
#'            xaxis = list(title = "Requested x points (cm)"),
#'            yaxis = list(title = "t (s)"))

mdfexdiffxtreq <- function(D,dt,l,xreq,treq,T,C_i=0,C_f=1,F=0.5){

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

  u_treq <- matrix(0,nrow=length(treq),ncol = ncol(u))

  for (j in 1:nrow(u_treq)) {
    for (i in 1:ncol(u)) {
      u_treq[j,i] <- pracma::interp1(t,u[,i],treq[j])
    }
  }

  u_xtreq <- matrix(0,nrow = nrow(u_treq), ncol = length(xreq))

  x_ins <- seq(-l,l,length.out = ncol(u_treq))

  for(i in 1:ncol(u_xtreq)){
    for(j in 1:nrow(u_xtreq)){
      u_xtreq[j,i] <- pracma::interp1(x_ins,u_treq[j,],xreq[i])
    }
  }

  return(u_xtreq)

}
