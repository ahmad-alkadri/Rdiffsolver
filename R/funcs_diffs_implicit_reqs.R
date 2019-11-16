#' 1D diffusion, boundary constant, FD implicit
#' for specific time (t) point(s)
#'
#' @description Function to solve diffusion equation
#' using finite difference method, backward euler or
#' implicit scheme, applied in the condition where:
#' a) the diffusion occurs inside a slab of material unidimensionally,
#' or only along one axis, b) the coefficient of diffusion is constant,
#' and c) the boundary conditions are constant. It will return
#' the solution values for the requested time point(s).
#'
#' @usage mdfimpdifftreq(D,Nx,Nt,l,T,C_i=0,C_f=1,treq)
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
#' @param treq requested time value(s)
#'
#' @return A matrix with {length(treq)} number of row and
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
#' treq = 1e4 # requested time in seconds
#' u <- mdfimpdifftreq(D,Nx,Nt,l,T,C_i,C_f,treq)
#'
#' # Using base plot
#' x <- seq(-l,l,length.out = Nx)
#' plot(x,u,type="b", xlab = "x (cm)", ylab = bquote(c["x,t"]))

mdfimpdifftreq <- function(D,Nx,Nt,l,T,C_i=0,C_f=1,treq){

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
  A = matrix(0, nrow = Nx, ncol = Nx)

  for(i in 2:Nx-1){
    A[i,i-1] = -F
    A[i,i+1] = -F
    A[i,i] = 1+2*F
  }

  # Boundary condition
  A[1,1] = 1
  A[1,2] = 0
  A[Nx,Nx] = 1
  A[Nx,Nx-1] = 0

  # Empty matrix for u and u_n
  u = matrix(0, nrow = Nx, ncol = Nt)
  # Boundary condition on the borders
  u[1,] = C_f[1]
  u[Nx,] = C_f[length(C_f)]
  # Initial condition
  u[2:(Nx-1),] = C_i

  # Compute
  for(n in 2:Nt){
    u_n = u[,n-1]
    b = matrix(0, nrow = Nx, ncol = 1)
    b[2:(Nx-1)] = u_n[2:(Nx-1)]
    b[1] = C_f
    b[Nx] = C_f[length(C_f)]
    ui = solve(A,b)
    u[,n] = ui
  }

  u = t(u)

  u_req <- matrix(0,nrow=length(treq),ncol = ncol(u))

  for (j in 1:nrow(u_req)) {
    for (i in 1:ncol(u)) {
      u_req[j,i] <- pracma::interp1(t,u[,i],treq[j])
    }
  }

  return(u_req)
}

#' 1D diffusion, boundary constant, FD implicit
#' for specific space (x) point(s)
#'
#' @description Function to solve diffusion equation
#' using finite difference method, backward euler or
#' implicit scheme, applied in the condition where:
#' a) the diffusion occurs inside a slab of material unidimensionally,
#' or only along one axis, b) the coefficient of diffusion is constant,
#' and c) the boundary conditions are constant. It will return
#' the solution values for the requested space (x) point(s).
#'
#' @usage mdfimpdiffxreq(D,Nx,Nt,l,T,C_i=0,C_f=1,xreq)
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
#' @param xreq requested x point(s)
#'
#' @return A matrix with Nt number of row and
#' {length(xreq)} number of column, profiling the diffusion.
#'
#' @examples
#' C_i = 0.00 # Initial concentration inside the slab
#' C_f = 1.00	# Final concentration coming from outside
#' D = 5e-7 # Coefficient of diffusion, cm^2/s
#' Nt = 100 # points in time
#' Nx = 20 # points in space
#' l = 0.25 # half-thickness of the slab, in cm
#' T = 432000 # Total measured time in seconds (~5 days)
#' xreq = 0 # requested x point
#' u <- mdfimpdiffxreq(D,Nx,Nt,l,T,C_i,C_f,xreq)
#'
#' # Using base plot
#' t <- seq(0,T,length.out = Nt)
#' plot(t, u, type="l",
#'      xlab = "t (seconds)",
#'      ylab = bquote(c["x,t"]))

mdfimpdiffxreq <- function(D,Nx,Nt,l,T,C_i=0,C_f=1,xreq){

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
  A = matrix(0, nrow = Nx, ncol = Nx)

  for(i in 2:Nx-1){
    A[i,i-1] = -F
    A[i,i+1] = -F
    A[i,i] = 1+2*F
  }

  # Boundary condition
  A[1,1] = 1
  A[1,2] = 0
  A[Nx,Nx] = 1
  A[Nx,Nx-1] = 0

  # Empty matrix for u and u_n
  u = matrix(0, nrow = Nx, ncol = Nt)
  # Boundary condition on the borders
  u[1,] = C_f[1]
  u[Nx,] = C_f[length(C_f)]
  # Initial condition
  u[2:(Nx-1),] = C_i

  # Compute
  for(n in 2:Nt){
    u_n = u[,n-1]
    b = matrix(0, nrow = Nx, ncol = 1)
    b[2:(Nx-1)] = u_n[2:(Nx-1)]
    b[1] = C_f
    b[Nx] = C_f[length(C_f)]
    ui = solve(A,b)
    u[,n] = ui
  }

  u = t(u)

  u_req <- matrix(0, nrow = nrow(u),
                  ncol = length(xreq))

  x_ins <- seq(-l,l,length.out = ncol(u))

  for(i in 1:ncol(u_req)){
    for(j in 1:nrow(u_req)){
      u_req[j,i] <- pracma::interp1(x_ins,u[j,],xreq[i])
    }
  }

  return(u_req)
}

#' 1D diffusion, boundary constant, FD implicit
#' for specific space (x) and time (t) point(s)
#'
#' @description Function to solve diffusion equation
#' using finite difference method, backward euler or
#' implicit scheme, applied in the condition where:
#' a) the diffusion occurs inside a slab of material unidimensionally,
#' or only along one axis, b) the coefficient of diffusion is constant,
#' and c) the boundary conditions are constant. It will return
#' the solution values for the requested space (x)
#' and time (t) point(s).
#'
#' @usage mdfimpdiffxtreq(D,Nx,Nt,l,T,C_i=0,C_f=1,xreq,treq)
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
#' @param xreq requested x point(s)
#' @param treq requested t point(s)
#'
#' @return A matrix with Nt number of row and
#' {length(xreq)} number of column, profiling the diffusion.
#'
#' @examples
#' C_i = 0.00 # Initial concentration inside the slab
#' C_f = 1.00	# Final concentration coming from outside
#' D = 5e-7 # Coefficient of diffusion, cm^2/s
#' Nt = 100 # points in time
#' Nx = 20 # points in space
#' l = 0.25 # half-thickness of the slab, in cm
#' T = 432000 # Total measured time in seconds (~5 days)
#' xreq = seq(-0.25,0,length.out = 20) # requested space (x) points
#' treq = seq(1e2,1e5,length.out = 40) # requested time (t) points
#' ureq <- mdfimpdiffxtreq(D,Nx,Nt,l,T,C_i,C_f,xreq,treq)
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

mdfimpdiffxtreq <- function(D,Nx,Nt,l,T,C_i=0,C_f=1,xreq,treq){

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
  A = matrix(0, nrow = Nx, ncol = Nx)

  for(i in 2:Nx-1){
    A[i,i-1] = -F
    A[i,i+1] = -F
    A[i,i] = 1+2*F
  }

  # Boundary condition
  A[1,1] = 1
  A[1,2] = 0
  A[Nx,Nx] = 1
  A[Nx,Nx-1] = 0

  # Empty matrix for u and u_n
  u = matrix(0, nrow = Nx, ncol = Nt)
  # Boundary condition on the borders
  u[1,] = C_f[1]
  u[Nx,] = C_f[length(C_f)]
  # Initial condition
  u[2:(Nx-1),] = C_i

  # Compute
  for(n in 2:Nt){
    u_n = u[,n-1]
    b = matrix(0, nrow = Nx, ncol = 1)
    b[2:(Nx-1)] = u_n[2:(Nx-1)]
    b[1] = C_f
    b[Nx] = C_f[length(C_f)]
    ui = solve(A,b)
    u[,n] = ui
  }

  u = t(u)

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
