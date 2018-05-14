c <- 40
chi <-
  sqrt(1 / c)                                              #x spacing
tau <-
  1 / c                                                    #Time spacing
vmax <-
  round(100 / tau + 1) * tau                                 #range of residence time v = (0, vmax)
xmax <-
  round(10 / chi) * chi                                     #range of x = (-xmax, xmax)
t <-
  round(100 / tau) * tau                                        #time at which evaluate density
#Initial grid (x,v) (v = residence time) with Initial condition = probability mass at one point (x = 0, v = 0) on grid
xi0 <- matrix(0, round(2 * xmax / chi + 1), round(vmax / tau + 1))
xi0[round(xmax / chi + 1), 1] <-
  1                              #Initial condition
xi <-
  xi0                                                     #iterate xi instead of xi0 in the loop below
alpha <-
  function(x)
    0.5 / (1 + exp(-x)) + 0.1  #varying exponent, alpha(-inf) = 0, alpha(+inf) = 1
tc <- 2
a <- function(x, s)
  (tc) ^ (-alpha(x))         #a(x,s) diffusivity
b <- function(x, s)
  0                        #b(x,s) drift function
Psi <-
  function(x, t) {
    r <-
      (t ^ (-alpha(x))) / gamma(1 - alpha(x))
    r[r > 1]
  }          #Tail function of psi measure
#Tail function probabilities at all grid points
TailF <- outer(seq(-xmax, xmax, chi), seq(0, vmax, tau), Psi) / c
TailF[TailF > 1] <- 1
#One-step Survival probabilities at all grid points
Sprob <-
  cbind(TailF[, 1], TailF[, seq(2, ncol(TailF), 1)] / TailF[, seq(1, ncol(TailF) -
                                                                    1, 1)])

k <-
  t / tau                                               #number of iterations
chk1 <-
  sum(xi0)                                         #check sum of probabilities in grid = 1
for (i in 1:k) {
  #jump probabilities for all x... Centre = 1 - a(x,t), Left = (a(x,t) - chi*b(x,t))/2, Right = (a(x,t) + chi*b(x,t))/2
  #First row contains Centre jump probabilities, second row contains Right jump probabilities, third row contains Left jump probabilities
  #Boundaries: At x = -xmax jump centre/right w.p. 0.5, At x = xmax jump centre/left w.p. 0.5
  Jprob <-
    rbind(c(0.5, 1 - a(seq(
      -xmax + chi, xmax - chi, chi
    ), i * tau)                                           , 0.5),
    c(0.5, (a(
      seq(-xmax + chi, xmax - chi, chi), i * tau
    ) + chi * b(
      seq(-xmax + chi, xmax - chi, chi), i * tau
    )) / 2 , 0),
    c(0, (a(
      seq(-xmax + chi, xmax - chi, chi), i * tau
    ) - chi * b(
      seq(-xmax + chi, xmax - chi, chi), i * tau
    )) / 2 , 0.5))
  
  S <-
    xi * Sprob                                            #amount that survive at each grid point
  E <-
    rowSums(xi * (1 - Sprob))                               #amount that escape at each x
  C_J <-
    E * Jprob[1, ]                                       #proportion of escapes at each x jumping centre
  R_J <-
    E * Jprob[2, ]                                       #proportion of escapes at each x jumping right
  L_J <-
    E * Jprob[3, ]                                       #proportion of escapes at each x jumping left
  
  #The amount that survives (S matrix) shifts right by one column (increase in residence time),
  #elements of the first column(v=0) calculated by summing amount that jumps centre and amount that jumps right/left from neighbouring rows,
  #elements of the last column (v=vmax) calculated by summing survivals S at both v = (vmax-tau, vmax)
  xi <-
    cbind(c(C_J[1] + L_J[2], R_J[seq(1, length(R_J) - 2, 1)] + C_J[seq(2, length(C_J) -
                                                                         1, 1)] + L_J[seq(3, length(L_J), 1)], R_J[length(R_J) - 1] + C_J[length(C_J)]), S[, seq(1, ncol(S) -
                                                                                                                                                                   2, 1)], S[, ncol(S) - 1] + S[, ncol(S)])
}
chk2 <-
  sum(xi)                                          #check sum of probabilities in grid = 1

#Plot density rho and v at time = 0 and t
par(mfrow = c(2, 2))
plot(
  seq(-xmax, xmax, chi),
  rowSums(xi0),
  type = "l",
  main = "Rho distribution at t = 0",
  xlab = "x",
  ylab = "Rho (x,0)"
)
grid()
plot(
  seq(0, vmax, tau),
  colSums(xi0),
  type = "l",
  main = "Residences times distribution at t = 0",
  xlab = "v",
  ylab = "Sum_{x} Xi (x,v,0)"
)
grid()
plot(
  seq(-xmax, xmax, chi),
  rowSums(xi) / chi,
  type = "l",
  main = paste("Rho distribution at t = ", t),
  xlab = "x",
  ylab = paste("Rho (x,", t, ")")
)
grid()
plot(
  seq(0, vmax, tau),
  colSums(xi),
  type = "l",
  main = paste("Residences times distribution at t = ", t),
  xlab = "v",
  ylab = paste("Sum_{x} Xi (x,v,", t, ")")
)
grid()
