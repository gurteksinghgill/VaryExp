source("R/varyExp-DTSM.R")
#Parameters
alpha <- function(x)
  0.55 * exp(-x^2) + 0.15 + 0.5 / (1 + exp(-2 * x))
Psi <- function(x,t) {
  if (t <= 0)
    Inf
  else
    (t^(-alpha(x))  )/gamma(1-alpha(x))
}

out <- DTSM(Psi = Psi, xrange = c(-1,1), snapshots = 2^(1:5), c = 25)
plot_DTSM_output(out)
