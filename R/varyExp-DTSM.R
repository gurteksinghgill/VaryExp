library(ggplot2)
library(tidyverse)

default_nuTail <- function(x, t) {
  if (t <= 0)
    Inf
  else
    t ^ (-0.7) / gamma(1 - 0.7)
}
default_a <- function(x,t) 
  1
default_b <- function(x,t)
  0
default_d <- function(x) 
  0.9

# arguments: 
#   xrange: a pair of numbers (xmin, xmax)
#   age_max: cutoff of age lattice
#   c: master scaling parameter
#   chi: spatial lattice spacing
#   tau: age lattice spacing
#   nuTail: space-dependent tail function of Levy measure
#   d: space-dependent temporal drift
# returns: a list of two items: 
#   * an (m,n) matrix of the CTRW states, with the first dimension
#     corresponding to space and the second to age
#   * an (m,n) matrix of survival probabilities
init_DTRM <- function(xrange,
                      age_max,
                      c,
                      chi,
                      tau,
                      nuTail = default_nuTail,
                      d = default_d) {
  # set up space-age-lattice
  m <- 2 * round((xrange[2] - xrange[1]) / (2 * chi)) + 1 # to make it an odd number
  n <- round(age_max / tau)
  xi0 <- matrix(0, m, n)
  # Put initial mass on center lattice point with age 0:
  midpoint_index <- (m + 1)/2
  xi0[midpoint_index, 1] <- 1 / chi
  
  # set up survival probability matrix
  nuTail_with_d <- function(x,t) {
    out <- min(1 - d(x), nuTail(x, t) / c)
    if (t <= tau)
      out <- out + d(x)
    out
  }
  x <- seq(from = xrange[1],
           to = xrange[2],
           length.out = m)
  age <- (1:(n+1)) * tau
  h <- outer(x, age, Vectorize(nuTail_with_d)) / c # see paper for definition of h
  h[h > 1] <- 1
  survival_probs <- h[ , -1] / h[ , -(n+1)]

  list(xi0 = xi0, survival_probs = survival_probs)
}

# calculates jump probabilities
# arguments: 
#   b, a: space- and time-dependent drift and diffusivity
#   x: a vector of locations
#   t: the current time (float)
# returns: 
#   a list with three items, each a vector of same length as x, 
#   for the probabilities to jump left, right and self-jumps
jump_probs <- function(x, t, a = default_a, b = default_b) {
  m <- length(x)
  chi <- diff(range(x)) / m
  a_vec <- sapply(x, function(x) a(x,t))
  b_vec <- sapply(x, function(x) b(x,t))
  left  <- (a_vec - chi * b_vec) / 2
  right <- (a_vec + chi * b_vec) / 2
  center <- 1 - a_vec
  # boundary conditions left end
  center[1] <- center[1] + left[1]
  left[1] <- 0
  #boundary conditions right end
  center[m] <- center[m] + right[m]
  right[m] <- 0
  list(left = left,
       center = center,
       right = right)
}

# evolve xi by one step
step_xi <- function(xi, Sprob, Jprob) {
  #Evaluate survivals, escapes and jumps
  surviving <- xi * Sprob
  escaping <- rowSums(xi - surviving)
  self_jumping  <- escaping * Jprob$center
  right_jumping <- escaping * Jprob$right
  left_jumping  <- escaping * Jprob$left
  
  #Update grid with survivals, escapes and jumps
  m <- dim(xi)[1]
  n <- dim(xi)[2]
  # increment age of surviving particles
  xi[ , -1] <- surviving[ , -n]
  # don't increment age of oldest particles
  xi[ , n] <- xi[ , n] + surviving[, n]
  # place escaping particles back on grid
  xi[ , 1] <- self_jumping + c(0, right_jumping[-m]) + c(left_jumping[-1], 0)
  xi
}


# Computes location-age densities for various snapshots in time.
# Arguments: 
#   xrange: 2-vector delineating the domain
#   snapshots: a vector of times at which location-age densities xi are computed
#   age_max: maximum age
#   c: master scalint parameter
#   chi: spatial grid parameter
#   tau: temporal grid parameter
#   a: diffusivity
#   b: drift
#   nuTail: tail function of Levy measure, truncated from above at 1
#   d: temporal drift
# Returns:
#   a list of the same length as snapthots, with containing the space-age
#   distributions xi   
DTSM <- function(xrange = c(-2, 2),
                 snapshots = c(0.5, 1, 2),
                 age_max = max(snapshots),
                 c = 1000,
                 chi = 1/sqrt(c),
                 tau = 1/c,
                 a = default_a,
                 b = default_b,
                 nuTail = default_nuTail,
                 d = default_d) {
  foo <-
    init_DTRM(
      xrange = xrange,
      age_max = age_max,
      c = c,
      chi = chi,
      tau = tau,
      nuTail = nuTail,
      d = d
    )
  xi   <- foo$xi0
  Sprob <- foo$survival_probs
  m <- dim(xi)[1]
  n <- dim(xi)[2]
  x <- seq(from = xrange[1],
           to = xrange[2],
           length.out = m)
  N <- length(snapshots)
  xi_list <- vector("list", N)
  t <- 0
  counter <- 0
  message("Need ", round(max(snapshots) / tau), " iterations.")
  for (i in 1:N) {
    while (t + tau <= snapshots[i]) {
      t <- t + tau
      counter <- counter + 1
      if (counter %% 1000 == 0)
        message("Finished ", counter, " iterations.")
      Jprob <- jump_probs(x = x, t = t)
      xi <- step_xi(xi = xi,
                    Sprob = Sprob,
                    Jprob = Jprob)
    }
    xi_list[[i]] <- xi
  }
  list(xi_list = xi_list, xrange = xrange, snapshots = snapshots)
}


plot_DTSM_output <- function(DTSM_output) {
  library(tibble)
  densities_df <- function(xi_list, xrange, snapshots) {
    N <- length(xi_list)
    xi <- xi_list[[1]]
    m <- dim(xi)[1]
    x <- seq(xrange[1], xrange[2], length.out = m)
    out <- tibble(x = x)
    for (i in 1:N) {
      rho <- rowSums(xi_list[[i]])
      out[paste0("t=", snapshots[i])] <- rho
    }
    out
  }
  xi_list <- DTSM_output[["xi_list"]]
  xrange <- DTSM_output[["xrange"]]
  snapshots <- DTSM_output[["snapshots"]]
  densities_df(xi_list, xrange, snapshots) %>%
    gather(snapshot, density,-x) %>%
    ggplot(aes(x = x, y = density, col = snapshot)) +
    geom_line()
}
