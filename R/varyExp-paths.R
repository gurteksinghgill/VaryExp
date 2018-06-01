c <- 40;
chi <- sqrt(1/c);                                                                   #x spacing
tau <- 1/c;                                                                         #Time spacing
t <- round(10/tau)*tau;                                                             #time at which evaluate density
alpha <- function(x) {r <- 0.5/(1+exp(-x)) + 0.25; return(r)};                      #varying exponent, alpha(-inf) = 0, alpha(+inf) = 1
tc <- 2;                                                                            #time scale
a <- function(x,s) {r <- (tc)^(-alpha(x)); return(r)};                              #a(x,s) diffusivity
b <- function(x,s) {r <- 0; return(r)};                                             #b(x,s) drift function

#Simulates path until time t and stores final position x in M
M <- rep(0,1000);
for (i in 1:1000) {
  X <- 0;
  T <- 0;
  while (T < t) {
    u1 <- runif(1);
    tt <- (u1*c*gamma(1-alpha(X)))^(-1/alpha(X));
    C <- 1 - a(X,T); L <- (a(X,T) - chi*b(X,T))/2; R <- (a(X,T) + chi*b(X,T))/2;
    T <- T + tt;
    if (T < t) {
      u2 <- runif(1);
        if (u2 < L) {
          X <- X - rnorm(1,0,1/sqrt(c));
        }
        else if (u2 > L && u2 < (L+R)) {
          X <- X + rnorm(1,0,1/sqrt(c));
        }
    }
  }
  M[i] = X;
}

#plot
par(mfrow=c(1,1));
hist(M,freq=FALSE, col="lightgray");
lines(density(M),col="blue");
grid();
