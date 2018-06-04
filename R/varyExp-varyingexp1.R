library(ggplot2)
library(tidyverse)

source("R/varyExp-DTSM.R")
#Parameters
alpha <- function(x)
  0.5 / (1 + exp(-x)) + 0.25
Psi <- function(x,t) {
  r <- (t^(-alpha(x))  )/gamma(1-alpha(x)) 
  r[r>1] <- 1
  r
}

#Set up grid size and spacing parameters
c <- 10
xrange <- c(-2,2)
chi <- 0.1 / sqrt(c)
tau <- 1/c
snapshots <- c(1,5,100)

out1 <- DTSM(xrange = xrange, snapshots = snapshots, c = c, chi = chi, tau = tau, Psi = Psi)

plot_DTSM_output(out1)


