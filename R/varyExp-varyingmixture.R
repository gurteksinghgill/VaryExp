#Parameters
xrange <- c(0,1);                                      
alpha1 <- 0.25;   
alpha2 <- 0.75;   
p <- function(x)  {r <- 1/(1+exp(-10*x + 5)); return(r)};             
Psi <- function(x,t) {r <- p(x)*(t^(-alpha1)  )/gamma(1-alpha1) + (1-p(x))*(t^(-alpha2)  )/gamma(1-alpha2); r[r>1] <- 1; return(r)};           
a <- function(x,s) {r <- (1)*rep(1,length(x)); return(r)};          
b <- function(x,s) {r <- 0; return(r)};                         
d <- function(x,s) {r <- 0; return(r)};                         

#Evaluate densities at multiple times
n1 <- DTSM(xrange, T = 1000, Psi, a, b, d);
# n2 <- DTSM(xrange, T = 2000, Psi, a, b, d);
# n5 <- DTSM(xrange, T = 5000, Psi, a, b, d);
# n10 <- DTSM(xrange, T = 10000, Psi, a, b, d);
# n15 <- DTSM(xrange, T = 15000, Psi, a, b, d);
# n25 <- DTSM(xrange, T = 25000, Psi, a, b, d);
# n100 <- DTSM(xrange, T = 100000, Psi, a, b, d);

#Plot densities
par(mfrow=c(1,1));
par(fig = c(0,1,0,1));
plot(n1$x, n1$Px, type="l", xlab = "x", ylab = paste("P(x,t)"));
# lines(n2$x, n2$Px, type="l", col="blue");
# lines(n5$x, n5$Px, type="l", col="pink");
# lines(n10$x, n10$Px, type="l", col="yellow");
# lines(n15$x, n15$Px, type="l", col="green");
# lines(n25$x, n25$Px, type="l", col="red");
# lines(n100$x, n100$Px, type="l", col="brown");
grid();
legend("topleft",c(expression(10^3),expression(2 %*% 10^3),expression(5 %*% 10^3),expression(10 %*% 10^3),expression(15 %*% 10^3),expression(25 %*% 10^3), expression(10^5) ), title = "Densities at t =", lty = c(1,1,1,1,1,1,1), col=c("black","blue","pink","yellow","green","red","brown") )
par(fig = c(0.55,0.99, 0.55, 0.99), new = T) 
plot(seq(0,1,0.01), p((seq(0,1,0.01))), xlab="x", ylab=expression(p (x)), type="l")

