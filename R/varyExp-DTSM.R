#Set up grid and initial condition
initialize_xi <- function(xrange, chi, tau, age_max) {
  xi0 <-
    matrix(0, round((xrange[2] - xrange[1]) / chi + 1), round(age_max / tau + 1))
  xi0[round(0.5 * (xrange[2] - xrange[1]) / chi + 1), 1] <- 1
  xi0
}

#Evaluate Tail function and one-step survival probabilities at every grid point
survival_probabilities <-
  function(Psi, xrange, c, chi, tau, age_max, d) {
    Psi_with_d <- function(x,t) Psi(x - d(x)/c, t)
    TailF <-
      outer(seq(xrange[1], xrange[2], chi), seq(0, age_max, tau), Psi_with_d) / c
    TailF[TailF > 1] <- 1
    cbind(1, TailF[, 2:ncol(TailF)] / TailF[, 1:(ncol(TailF)-1)])
  }

jump_probabilities <- function(a, b, xrange, chi, t){
  left <- c(0, 
            (a(seq(xrange[1] + chi, xrange[2] - chi, chi), t) 
              - chi * b(seq(xrange[1] + chi, xrange[2] - chi, chi), t)) / 2, 
            0.5)
  center <- c(0.5, 
              1 - a(seq(xrange[1] + chi, xrange[2] - chi, chi), t), 
              0.5)
  right <- c(0.5, 
             (a(seq(xrange[1] + chi, xrange[2] - chi, chi), t) 
               + chi * b(seq(xrange[1] + chi, xrange[2] - chi, chi), t)) / 2, 
             0)
  list(left = left, center = center, right = right)
}
# evolve xi by one step
step_xi <- function(xi, i, Sprob) {
  #Evaluate jump probabilities at every spatial grid point + Set boundary condition
  Jprob <- jump_probabilities(a, b, xrange, chi, i * tau)
  
  #Evaluate survivals, escapes and jumps
  surviving <- xi * Sprob
  escaping <- rowSums(xi * (1 - Sprob))
  self_jumping  <- escaping * Jprob$center
  right_jumping <- escaping * Jprob$right
  left_jumping  <- escaping * Jprob$left
  
  #Update grid with survivals, escapes and jumps
  I <- dim(xi)[1]
  J <- dim(xi)[2]
  # save density of particles with maximum age
  oldest <- xi[, J]
  # increment age of surviving particles
  xi[ , 2:J] <- surviving[ , 1:J-1]
  xi[ , J] <- xi[ , J] + oldest
  # place escaping particles back on grid
  xi[ , 1] <- self_jumping + c(0, right_jumping[-I]) + c(left_jumping[-1], 0)
  xi
}

#The DTSM function implements a discrete time semi-Markov algorthim for calculating CTRW limit distributions
#Inputs are range of x (e.g. [0,1]), time T, Tail function Psi and coefficient functions a, b and d
#Outputs are 4 plots of density and age density (at times 0 and T) and a list containing the (x,y) values of the plot of the density at time = T


DTSM <- function (xrange, T, Psi, a, b, d) {
  xi <- initialize_xi(xrange, chi, tau, age_max)
  Sprob <- survival_probabilities(Psi, xrange, c, chi, tau, age_max, d)
  #Iterate master equation on grid
  k <- round(T/tau)
  chk1 <- sum(xi)
  for (i in 1:k) {
    xi <- step_xi(xi, i, Sprob)
    chk2 <- sum(xi)
  }
  list(x = seq(xrange[1],xrange[2],chi), Px = rowSums(xi)/chi)
}


