DTSM <- function (Q) {
c <- 100;
chi <- (2/c);                                               #x spacing
tau <- (10/c);                                                     #Time spacing
vmax <- round(100/tau + 1)*tau;                                  #range of residence time v = (0, vmax)
xmax <- round(1/chi)*chi;                                      #range of x = (-xmax, xmax)
t <- round(Q/tau)*tau;                                         #time at which evaluate density
#Initial grid (x,v) (v = residence time) with Initial condition = probability mass at one point (x = 0, v = 0) on grid
xi0 <- matrix(0,round(xmax/chi + 1), round(vmax/tau + 1)); 
xi0[round(0.5*xmax/chi + 1), 1] <- 1;                               #Initial condition
xi <- xi0;                                                      #iterate xi instead of xi0 in the loop below
alpha <- function(x) {r <- 0.9*rep(1,length(x)); return(r)};   #varying exponent, alpha(-inf) = 0, alpha(+inf) = 1
tc <- 1;
a <- function(x,s) {r <- (tc)^(-alpha(x)); return(r)};          #a(x,s) diffusivity
b <- function(x,s) {r <- 0; return(r)};                         #b(x,s) drift function
p <- function(x)  {r <- 1/(1+exp(-10*x + 5)); return(r)};             #{r <- 0.98*x + 0.01; r[r>0.99] <- 0.99; return(r)};
Psi <- function(x,t) {r <- p(x)*(t^(-0.25)  )/gamma(1-0.25) + (1-p(x))*(t^(-0.75)  )/gamma(1-0.75); r[r>1] <- 1; return(r)};           #Tail function of psi measure
#Tail function probabilities at all grid points
TailF <- outer(seq(0,xmax,chi),seq(0,vmax,tau),Psi)/c;
TailF[TailF > 1] <- 1;
#One-step Survival probabilities at all grid points
Sprob <- cbind(1, TailF[, seq(2,ncol(TailF),1)]/TailF[, seq(1,ncol(TailF)-1,1)] ); 

k <- t/tau;                                                #number of iterations
chk1 <- sum(xi0);                                          #check sum of probabilities in grid = 1
for (i in 1:k) {
  #jump probabilities for all x... Centre = 1 - a(x,t), Left = (a(x,t) - chi*b(x,t))/2, Right = (a(x,t) + chi*b(x,t))/2
  #First row contains Centre jump probabilities, second row contains Right jump probabilities, third row contains Left jump probabilities 
  #Boundaries: At x = -xmax jump centre/right w.p. 0.5, At x = xmax jump centre/left w.p. 0.5
  Jprob <- rbind(c(0.5, 1 - a(seq(chi,xmax-chi,chi),i*tau)                                           ,0.5),
                 c(0.5, (a(seq(chi,xmax-chi,chi),i*tau) + chi*b(seq(chi,xmax-chi,chi),i*tau))/2 , 0), 
                 c(0, (a(seq(chi,xmax-chi,chi),i*tau) - chi*b(seq(chi,xmax-chi,chi),i*tau))/2 , 0.5) );
  
  S <- xi*Sprob;                                             #amount that survive at each grid point
  E <- rowSums(xi*(1-Sprob));                                #amount that escape at each x
  C_J <- E*Jprob[1,];                                        #proportion of escapes at each x jumping centre
  R_J <- E*Jprob[2,];                                        #proportion of escapes at each x jumping right
  L_J <- E*Jprob[3,];                                        #proportion of escapes at each x jumping left
  
  #The amount that survives (S matrix) shifts right by one column (increase in residence time),
  #elements of the first column(v=0) calculated by summing amount that jumps centre and amount that jumps right/left from neighbouring rows, 
  #elements of the last column (v=vmax) calculated by summing survivals S at both v = (vmax-tau, vmax)
  xi <- cbind( c( C_J[1]+L_J[2], R_J[seq(1,length(R_J)-2,1)]+C_J[seq(2,length(C_J)-1,1)]+L_J[seq(3,length(L_J),1)], R_J[length(R_J)-1]+C_J[length(C_J)] ), S[, seq(1,ncol(S)-2,1)], S[, ncol(S)-1]+S[, ncol(S)] );
}
chk2 <- sum(xi);                                           #check sum of probabilities in grid = 1

#Plot density rho and v at time = 0 and t
# par(mfrow=c(2,2));
# plot(seq(0,xmax,chi),rowSums(xi0),type="l", main = "Rho distribution at t = 0", xlab = "x", ylab = "Rho (x,0)");
# grid();
# plot(seq(0,vmax,tau),colSums(xi0),type="l", main = "Residences times distribution at t = 0", xlab = "v", ylab = "Sum_{x} Xi (x,v,0)");
# grid();
# plot(seq(0,xmax,chi),rowSums(xi)/chi,type="l", main = paste("Rho distribution at t = ", t), xlab = "x", ylab = paste("Rho (x,", t, ")"));
# grid();
# plot(seq(0,vmax,tau),colSums(xi),type="l", main = paste("Residences times distribution at t = ", t), xlab = "v", ylab = paste("Sum_{x} Xi (x,v,", t, ")"));
# grid();

return(rowSums(xi)/chi);

}


