#' One-dimensional diffusion, constant boundary conditions, without messages
#'
#' @description Function to solve diffusion equation using numerical implicit finite differences method,
#' applied in the condition where:
#' a) the diffusion occurs inside a slab of material unidimensionally, or only along one axis,
#' b) the coefficient of diffusion is constant, and
#' c) the boundary conditions are constant.
#' The only difference between this function and the diff.1D() is that this one does not return any messages with notice about the boundary conditions at the end.
#'
#' @usage diff_1D_nm(Lx, Tt, nt, nx, D, C_ini, C_lim)
#'
#' @param Lx Total length of slab in x direction, usually in mm.
#' @param Tt Total diffusion time, usually in second.
#' @param nt Number of time discretization.
#' @param nx Number of dimension x discretization
#' @param D Coefficient of diffusion. Must be constant value.
#' @param C_ini Initial concentration value inside the slab.
#' @param C_lim Boundary condition, written as a vector with one or two elements.
#' If there is only one element, the other side of the slab will be presented as having Neumann boundary condition with flux = 0.
#' If there are two elements, the first element is the dirichlet concentration on the left side, while the second element is the dirichlet concentration on the right side.
#'
#' @return A matrix with {nx+1} number of row and {nt} number of column, profiling the diffusion on slab with Lx length along the time Tt.
#'
#' @examples
#' Lx <- 5 #Length of slab, in mm
#' Tt <- 20 #Total measured diffusion time, in seconds
#' nt <- 100 #Number of time discretization, wherein the Tt will be divided into 100 equal parts (0, 0.2, 0.4, ..., 20)
#' nx <- 50 #Number of dimension discretization
#' D <- 0.05 #Coefficient of diffusion in mm^2/s.
#' C_ini <- 0.05 #Initial concentration inside the slab.
#' C_lim <- c(0.10,0.25) #Dirichlet boundary concentration diffusing into the slab
#' matC <- diff_1D_nm(Lx, Tt, nt, nx, D, C_ini, C_lim)

diff_1D_nm <- function(Lx,Tt,nt,nx,D,C_ini,C_lim){

  #1-Dimensional diffusion solver using implicit method, with border condition on its two extremities/sides

  ##Basic parameters
  if(length(C_lim) == 1){

    C_fin <- C_lim

  } else {

    if(length(C_lim == 2)){
      C_left <- C_lim[1]

      C_right <- C_lim[2]

    } else {

      stop('ERROR = Unfortunately, this function cannot take more than two condition limits coming from each side')

    }

  }

  dt <- Tt/nt #Size of one step in t

  dx <- Lx/(nx+1) #Size of one step in x

  ##Calculations
  s = D*dt/(dx^2)

  ##
  # Matrix A
  A = matrix(0, nrow = nx+1, ncol = nx+1)

  for (i in 2:nx) {

    A[i, i] <- (1+2*s) #Diagonal elements

    A[i, (i+1)] <- -s

    A[(i+1), i] <- -s

  }

  A[1,1] <- (1+2*s); A[nx+1,nx+1] <- A[1,1]

  A[2,1] <- -s; A[1,2] <- -s

  b <- matrix(0, nrow = nx+1, ncol = 1)

  #Conditioning
  if(length(C_lim) == 2){

    b[1] <- s*C_right #Final, right side

    b[nx+1] <- s*C_left #Final, left side

  } else {

    A[nx+1,nx] <- -2*s #For neumann boundary condition

    b[1] <- s*C_fin #Final, right side

    b[nx+1] <- 0 #Zero flux, neumann condition, left side

  }

  Cmat_i <- matrix(C_ini, nrow = nx+1, ncol = 1)

  Cmat_f <- matrix(0, nrow = nx+1, ncol = nt)

  t_n <- matrix(0, nrow = nt, ncol = 1)

  Cmat_f[,1] <- Cmat_i

  t_j <- 0 #Initial time

  t_n[1] <- t_j

  ######
  ##March solution in time

  for (j in 1:nt) {

    t_j <- t_j+dt

    if(j == 1){
      Cmat_j <- Cmat_i
    } else {
      Cmat_j <- Cmat_j_1
    }

    c <- Cmat_j + b

    Cmat_j_1 = solve(A,c)

    Cmat_f[,j] <- Cmat_j_1

    #Next time step
    t_n[j] <- t_j

  }

  return(Cmat_f)

}
