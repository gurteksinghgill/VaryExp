source("R/varyExp-DTSM.R")
#Parameters
alpha <- function(x)
  0.5 / (1 + exp(-x)) + 0.25
Psi <- function(x,t) {
  if (t <= 0)
    Inf
  else 
    (t^(-alpha(x))  )/gamma(1-alpha(x)) 
}

out1 <- DTSM(Psi = Psi, xrange = c(-4,4), snapshots = c(1,2,4))
plot_DTSM_output(out1)
