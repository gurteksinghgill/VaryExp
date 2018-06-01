#Parameters
xrange <- c(-1,1);                                      
alpha <- function(x) {r <- 0.5/(1+exp(-x)) + 0.25; return(r)};   
theta <- function(x) {r <- 1/(1+exp(0.5*x)); return(r)};   
Psi <- function(x,t) {r <- (t^(-alpha(x)) * exp(-t*theta(x)) )/gamma(1-alpha(x)); r[r>1] <- 1; return(r)};           
a <- function(x,s) {r <- (1)^(-alpha(x)); return(r)};          
b <- function(x,s) {r <- 0; return(r)};                         
d <- function(x,s) {r <- 0; return(r)};                         

#Evaluate densities at multiple times
m1 <- DTSM(xrange, T = 1000, Psi, a, b, d);
m2 <- DTSM(xrange, T = 2000, Psi, a, b, d);
m5 <- DTSM(xrange, T = 5000, Psi, a, b, d);
m10 <- DTSM(xrange, T = 10000, Psi, a, b, d);
m15 <- DTSM(xrange, T = 15000, Psi, a, b, d);
m25 <- DTSM(xrange, T = 25000, Psi, a, b, d);

#Plot densities
par(mfrow=c(1,1));
par(fig = c(0,1,0,1));
plot(m1$x, m1$Px, type="l", xlab = "x", ylab = paste("P(x,t)"));
lines(m2$x, m2$Px, type="l", col="blue");
lines(m5$x, m5$Px, type="l", col="pink");
lines(m10$x, m10$Px, type="l", col="yellow");
lines(m15$x, m15$Px, type="l", col="green");
lines(m25$x, m25$Px, type="l", col="red");
grid();
legend("topright",c(expression(10^3),expression(2 %*% 10^3),expression(5 %*% 10^3),expression(10 %*% 10^3),expression(15 %*% 10^3),expression(25 %*% 10^3) ), title = "Densities at t =", lty = c(1,1,1,1,1,1), col=c("black","blue","pink","yellow","green","red") )
par(fig = c(0.25,0.75, 0.05, 0.5), new = T) 
plot(seq(-5,5,0.01), theta((seq(-5,5,0.01))), xlab="x", ylab=expression(theta (x)), type="l")

