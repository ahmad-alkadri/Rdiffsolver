#' 1D diffusion, boundary constant, FD newtonian
#' method, concentration-dependent diffusion coefficient
#'
#' @description Function to solve diffusion equation
#' using finite difference method, backward euler or
#' implicit scheme, applied in the condition where:
#' a) the diffusion occurs inside a slab of material unidimensionally,
#' or only along one axis, b) the coefficient of diffusion
#' is dependent on the concentration (D = f(c)),
#' and c) the boundary conditions are constant.
#'
#' @usage mdfnewtondiffu(Nt,Nx,T,C_f,C_i,l,Dfunstr)
#'
#' @param Nt number of points in time (t)
#' @param Nx number of points in space (x)
#' @param T Total calculated diffusion time
#' @param C_i Initial concentration value inside the slab
#' @param C_f Dirichlet boundary condition,
#' final concentration coming from the outside of the slab,
#' with one or two elements. If there is only one element,
#' the two sides of the slab (x = -l and x = l) will have
#' the same C_f. If there are
#' two elements, the first element (C_f[1]) will be on
#' x = -l while the second one will be on x = l.
#' @param l half-length of the slab, usually in cm
#' @param Dfunstr the coefficient of diffusion written as
#' a function of concentration (c) as string;
#' please check the example for details
#'
#' @return A matrix with Nx number of row and
#' Nt number of column, profiling the diffusion.
#'
#' @examples
#' C_i = 0.00 # Initial concentration inside the slab
#' C_f = c(0.10,0.07)	# Final concentration coming from outside
#' Nt = 30 # points in time
#' Nx = 20 # points in space
#' l = 0.25 # half-thickness of the slab, in cm
#' T = 432000 # Total measured time in seconds (~5 days)
#'
#' # Coefficient of diffusion as a linear function of c
#' penta = 2.21e-7 #cm^2/s
#' bro = 1.91e-7 #cm^2/s
#' Dfunstr = paste("function(c){",penta,"*c+",bro,"}",sep="")
#' print(Dfunstr)
#'
#' # Calculation
#' U <- mdfnewtondiffu(Nt,Nx,T,C_f,C_i,l,Dfunstr)
#'
#' # Visualization
#' x_s = seq(-l,l,length.out=Nx)
#' y_s = seq(0,T,length.out=Nt)
#' persp(x=x_s, y=y_s, z=U,
#'       theta = -35, phi = 25,
#'       ticktype = "detailed",
#'       xlab = "x (mm)",
#'       ylab = "t (seconds)",
#'       zlab = expression(c[x,t]))

mdfnewtondiffu <- function(Nt,Nx,T,C_f,C_i,l,Dfunstr){

  if(length(C_f) > 2){
    stop("ERROR: C_f must consist of maximum 2 elements")
  }

  a = 0
  b = 2*l

  # FIRST STEPS

  # Left boundary condition
  alpha = C_f[1]
  # Right boundary condition
  beta = C_f[length(C_f)]
  tau=T/(Nt-1)
  h=(b-a)/(Nx-1)
  x = matrix(0,Nx,1)
  u0 = matrix(0,Nx,1)
  for(i in 1:Nx){
    x[i] = a+(i-1)*h
    u0[i] = C_i
  }
  t = matrix(0,Nt,1)
  for(n in 1:Nt){
    t[n] = (n-1)*tau
  }
  u_1 = matrix(0,Nx,1)
  u = matrix(0,Nx,1)
  uNxext = matrix(0,Nx,1)
  G = matrix(0,Nx,1)
  L = matrix(0,Nx,Nx)
  L[1,1] = 1
  L[Nx,Nx] = 1

  # EMPTY MATRIX FOR RESULTS
  U = matrix(0,Nx,Nt)
  U[,1] = u0

  # CALCULATIONS
  for(n in 2:Nt){
    u=U[,(n-1)]
    u_1 = U[,(n-1)]
    eps = 1
    while(eps > 0.0001){
      G[1] = u[1]-alpha
      G[Nx] = u[Nx]-beta
      for(i in 2:(Nx-1)){
        dfun = eval(parse(text = Dfunstr))
        k = dfun(c=u[i])
        Dk = pracma::fderiv(dfun,
                            u[i], n=1)
        D2k = pracma::fderiv(dfun,
                             u[i], n=2)
        v = (u[i+1]-u[i-1])/(2*h)
        A = 1/tau
        phi = A*(u[i]-u_1[i])-Dk*v*v
        f = phi/k
        q=(-f*Dk+A-D2k*v*v)/k
        p=-2*Dk*v/k
        G[i]=u[i+1]-2*u[i]+u[i-1]-h*h*f;
        L[i,i-1]=1+0.5*h*p;
        L[i,i]=-2-h*h*q;
        L[i,i+1]=1-0.5*h*p;
      }
      uNxext = u-solve(L,G)
      eps = sqrt(h*t(uNxext-u) %*% (uNxext-u))
      u = uNxext
    }
    U[,n] = u
  }
  return(U)
}
