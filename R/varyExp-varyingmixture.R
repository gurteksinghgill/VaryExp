source("R/varyExp-DTSM.R")
#Parameters
alpha1 <- 0.25
alpha2 <- 0.75
p <- function(x)  {r <- 1/(1+exp(-10*x + 5)); return(r)}
Psi <- function(x,t) {
  if (t <= 0)
    Inf
  else
    p(x) * (t ^(-alpha1)) / gamma(1 - alpha1) + (1 - p(x)) * (t ^(-alpha2)) / gamma(1 - alpha2)
}

out <- DTSM(Psi = Psi, xrange = c(0,1), snapshots = 2^(1:5), c = 25)
plot_DTSM_output(out)
