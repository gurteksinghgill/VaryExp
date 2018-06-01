#Parameters
xrange <- c(-1,1);                                      
alpha <- function(x) {r <- 0.55*exp(-x^2) + 0.15 + 0.5/(1+exp(-2*x)); return(r)};   
Psi <- function(x,t) {r <- (t^(-alpha(x))  )/gamma(1-alpha(x)); r[r>1] <- 1; return(r)};           
a <- function(x,s) {r <- (1)^(-alpha(x)); return(r)};          
b <- function(x,s) {r <- 0; return(r)};                        
d <- function(x,s) {r <- 0; return(r)};                         

#Evaluate densities at multiple times
k0025 <- DTSM(xrange, T = 25, Psi, a, b, d);
k0050 <- DTSM(xrange, T = 50, Psi, a, b, d);
k0075 <- DTSM(xrange, T = 75, Psi, a, b, d);
k010 <- DTSM(xrange, T = 100, Psi, a, b, d);
k025 <- DTSM(xrange, T = 250, Psi, a, b, d);
k050 <- DTSM(xrange, T = 500, Psi, a, b, d);
k075 <- DTSM(xrange, T = 750, Psi, a, b, d);
k1 <- DTSM(xrange, T = 1000, Psi, a, b, d);
k2 <- DTSM(xrange, T = 2000, Psi, a, b, d);
k5 <- DTSM(xrange, T = 5000, Psi, a, b, d);
k10 <- DTSM(xrange, T = 10000, Psi, a, b, d);
k15 <- DTSM(xrange, T = 15000, Psi, a, b, d);
k25 <- DTSM(xrange, T = 25000, Psi, a, b, d);

#Plot densities
par(mfrow=c(1,1));
par(fig = c(0,1,0,1));
plot(k0050$x, k0050$Px, type="l", xlab = "x", ylab = paste("P(x,t)"), ylim=c(0,3.5));
lines(k010$x, k010$Px, type="l", col="blue");
lines(k050$x, k050$Px, type="l", col="green");
lines(k1$x, k1$Px, type="l", col="yellow");
lines(k10$x, k10$Px, type="l", col="red");
grid();
legend("topright",c(expression(50),expression(10^2),expression(5 %*% 10^2),expression(10^3),expression(10^4) ), title = "Densities at t =", lty = c(1,1,1,1,1), col=c("black","blue","green","yellow","red") )
par(fig = c(0.25,0.75, 0.55, 0.99), new = T) 
plot(seq(-5,5,0.01), alpha((seq(-5,5,0.01))), xlab="x", ylab=expression(alpha (x)), type="l")

