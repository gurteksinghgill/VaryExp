#Parameters
xrange <- c(-2,2);                                      
alpha <- function(x) {r <- 0.5/(1+exp(-x)) + 0.25; return(r)};   
Psi <- function(x,t) {r <- (t^(-alpha(x))  )/gamma(1-alpha(x)); r[r>1] <- 1; return(r)};           
a <- function(x,s) {r <- (1)^(-alpha(x)); return(r)};          
b <- function(x,s) {r <- 0; return(r)};                         
d <- function(x,s) {r <- 0; return(r)};                         

#Evaluate densities at multiple times 
p1 <- DTSM(xrange, T = 1000, Psi, a, b, d);
# p2 <- DTSM(xrange, T = 2000, Psi, a, b, d);
# p5 <- DTSM(xrange, T = 5000, Psi, a, b, d);
# p10 <- DTSM(xrange, T = 10000, Psi, a, b, d);
# p15 <- DTSM(xrange, T = 15000, Psi, a, b, d);
# p25 <- DTSM(xrange, T = 25000, Psi, a, b, d);

#plot densities
par(mfrow=c(1,1));
par(fig = c(0,1,0,1));
plot(p1$x, p1$Px, type="l", xlab = "x", ylab = paste("P(x,t)"));
# lines(p2$x, p2$Px, type="l", col="blue");
# lines(p5$x, p5$Px, type="l", col="pink");
# lines(p10$x, p10$Px, type="l", col="yellow");
# lines(p15$x, p15$Px, type="l", col="green");
# lines(p25$x, p25$Px, type="l", col="red");
grid();
legend("topleft",c(expression(10^3),expression(2 %*% 10^3),expression(5 %*% 10^3),expression(10 %*% 10^3),expression(15 %*% 10^3),expression(25 %*% 10^3) ), title = "Densities at t =", lty = c(1,1,1,1,1,1), col=c("black","blue","pink","yellow","green","red") )
par(fig = c(0.55,0.99, 0.55, 0.99), new = T) 
plot(seq(-5,5,0.01), alpha((seq(-5,5,0.01))), xlab="x", ylab=expression(alpha (x)), type="l")
