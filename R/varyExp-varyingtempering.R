#Parameters
alpha <- function(x)
  0.5 / (1 + exp(-x)) + 0.25
theta <- function(x)
  1 / (1 + exp(0.5 * x))
Psi <- function(x,t) {
  if (t <= 0)
    Inf
  else
    (t^(-alpha(x)) * exp(-t*theta(x)) )/gamma(1-alpha(x))
}

out <- DTSM(Psi = Psi, xrange = c(-1,1), snapshots = 2^(1:5), c = 25)
plot_DTSM_output(out)
