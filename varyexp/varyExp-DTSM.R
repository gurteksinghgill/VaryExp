#The DTSM function implements a discrete time semi-Markov algorthim for calculating CTRW limit distributions
#Inputs are range of x (e.g. [0,1]), time T, Tail function Psi and coefficient functions a, b and d
#Outputs are 4 plots of density and age density (at times 0 and T) and a list containing the (x,y) values of the plot of the density at time = T

DTSM <- function (xrange, T, Psi, a, b, d) {

  #Set up grid size and spacing parameters
  c <- 100;
  chi <- (3.125/c);                                               
  tau <- (10/c);                                                 
  age_max <- round(100/tau + 1)*tau;                             
  xrange <- round(xrange/chi)*chi;

  #Set up grid and initial condition
  xi0 <- matrix(0,round((xrange[2]-xrange[1])/chi + 1), round(age_max/tau + 1)); 
  xi0[round(0.5*(xrange[2]-xrange[1])/chi + 1), 1] <- 1;      
  xi <- xi0;                                                  
  
  #Evaluate Tail function and one-step survival probabilities at every grid point
  TailF <- outer(seq(xrange[1],xrange[2],chi),seq(0,age_max,tau) - tau*d(0,0),Psi)/c;
  TailF[TailF > 1] <- 1;
  Sprob <- cbind(1, TailF[, seq(2,ncol(TailF),1)]/TailF[, seq(1,ncol(TailF)-1,1)] ); 
  
  #Iterate master equation on grid
  k <- round(T/tau);
  chk1 <- sum(xi0);                                          
  for (i in 1:k) {
    #Evaluate jump probabilities at every spatial grid point + Set boundary condition 
    Jprob <- rbind(c(0.5, 1 - a(seq(xrange[1]+chi,xrange[2]-chi,chi),i*tau), 0.5),
                   c(0.5, (a(seq(xrange[1]+chi,xrange[2]-chi,chi),i*tau) + chi*b(seq(xrange[1]+chi,xrange[2]-chi,chi),i*tau))/2, 0), 
                   c(0, (a(seq(xrange[1]+chi,xrange[2]-chi,chi),i*tau) - chi*b(seq(xrange[1]+chi,xrange[2]-chi,chi),i*tau))/2, 0.5) );
    
    #Evaluate survivals, escapes and jumps
    S <- xi*Sprob;                                             
    E <- rowSums(xi*(1-Sprob));                                
    C_J <- E*Jprob[1,];                                        
    R_J <- E*Jprob[2,];                                        
    L_J <- E*Jprob[3,];
    
    #Update grid with survivals, escapes and jumps
    xi <- cbind( c( C_J[1]+L_J[2], R_J[seq(1,length(R_J)-2,1)]+C_J[seq(2,length(C_J)-1,1)]+L_J[seq(3,length(L_J),1)], R_J[length(R_J)-1]+C_J[length(C_J)] ), S[, seq(1,ncol(S)-2,1)], S[, ncol(S)-1]+S[, ncol(S)] );
  }
  chk2 <- sum(xi);
  
  #Plot density and age density at time = 0 and T
  par(mfrow=c(2,2));
  plot(seq(xrange[1],xrange[2],chi),rowSums(xi0),type="l", main = "Density at T = 0", xlab = "x", ylab = "P(x,0)");
  grid();
  plot(seq(0,age_max,tau),colSums(xi0),type="l", main = "Age Density at T = 0", xlab = "t", ylab = "P(V=t)");
  grid();
  plot(seq(xrange[1],xrange[2],chi),rowSums(xi)/chi,type="l", main = paste("Density at T = ", T), xlab = "x", ylab = paste("P(x,", T, ")"));
  grid();
  plot(seq(0,age_max,tau),colSums(xi),type="l", main = paste("Age Density at T = ", T), xlab = "t", ylab = "P(V=t)");
  grid();
  
  return(list(x = seq(xrange[1],xrange[2],chi), Px = rowSums(xi)/chi));
  
}


