n = 100;
chi = sqrt(1/n);                                          %x spacing
tau = 1/n;                                                %Time spacing
vmax = round(10/tau + 1)*tau;                             %range of residence time v = (0, vmax)
xmax = round(100/chi)*chi;                                  %range of x = (-xmax, xmax)
t = round(10/tau)*tau;                                    %time at which evaluate density
%Initial grid (x,v) (v = residence time) with Initial condition = probability mass at one point (x = 0, v = 0) on grid
xi0 = zeros(round(2*xmax/chi + 1), round(vmax/tau + 1)); 
xi0(round(xmax/chi + 1), 1) = 1;                          %Initial condition
xi = xi0;                                                 %iterate xi instead of xi0 in the loop below
alpha = @(x) atan(x)/pi + 1/2;                            %varying exponent, alpha(-inf) = 0, alpha(+inf) = 1
tc = 2;
a = @(x,s) (tc).^(-alpha(x));                             %a(x,s) diffusivity
b = @(x,s) 0;                                             %b(x,s) drift function
Psi = @(t,x) (t.^(-alpha(x)))./gamma(1-alpha(x));         %Tail function of psi measure
%Tail function probabilities at all grid points
TailF = min(1, Psi(repmat((tau:tau:vmax),size(xi,1),1), repmat((-xmax:chi:xmax)',1,size(xi,2)-1))/n );  
%One-step Survival probabilities at all grid points
Sprob = [TailF(:,1), TailF(:, 2:end)./TailF(:, 1:end-1), TailF(:,end)./TailF(:,end-1)]; 


k = t/tau;                                                %number of iterations
chk1 = sum(sum(xi0));                                     %check sum of probabilities in grid = 1
for i = 1:k
    %jump probabilities for all x... Centre = 1 - a(x,t), Left = (a(x,t) - chi*b(x,t))/2, Right = (a(x,t) + chi*b(x,t))/2
    %First row contains Centre jump probabilities, second row contains Right jump probabilities, third row contains Left jump probabilities 
    %Boundaries: At x = -xmax jump centre/right w.p. 0.5, At x = xmax jump centre/left w.p. 0.5
    Jprob = [0.5, 1 - a(-(xmax-chi):chi:(xmax-chi),i*tau)                                       , 0.5;
             0.5, (a(-(xmax-chi):chi:(xmax-chi),i*tau) + chi*b(-(xmax-chi):chi:(xmax-chi),i*tau))/2 , 0; 
             0, (a(-(xmax-chi):chi:(xmax-chi),i*tau) - chi*b(-(xmax-chi):chi:(xmax-chi),i*tau))/2 , 0.5];
    
    S = xi.*Sprob;                                        %amount that survive at each grid point
    E = sum(xi.*(1-Sprob), 2);                            %amount that escape at each x
    C_J = E.*Jprob(1,:)';                                 %proportion of escapes at each x jumping centre
    R_J = E.*Jprob(2,:)';                                 %proportion of escapes at each x jumping right
    L_J = E.*Jprob(3,:)';                                 %proportion of escapes at each x jumping left
    
    %The amount that survives in S matrix shifts right by one column (increase in residence time),
    %elements of the first column(v=0) calculated by summing amount that jumps centre and amount that jumps right/left from neighbouring rows, 
    %elements of the last column (v=vmax) calculated by summing survivals S at both v = (vmax-tau, vmax)
    xi = [ [ C_J(1)+L_J(2); R_J(1:end-2)+C_J(2:end-1)+L_J(3:end); R_J(end-1)+C_J(end) ], S(:, 1:end-2), S(:, end-1)+S(:, end)];
end
chk2 = sum(sum(xi));                                      %check sum of probabilities in grid = 1


%Plot density Xi at t = 0
figure(1);
surf(xi0);        
zlabel('Xi (x,v,0)');   
title('Initial Xi distribution at t = 0');
%Plot density Xi at t
figure(2);
surf(xi);
str1 = sprintf('Xi (x,v,%d)',t);
zlabel(str1);   
str2 = sprintf('Xi distribution at t = %d', t);
title(str2);

%Plot density rho and v at t = 0
figure(3);
subplot(2,1,1);
plot(-xmax:chi:xmax,sum(xi0, 2));
grid on;
xlabel('x');
ylabel('Rho (x,0)');
title('Initial Rho distribution across x at t = 0');
subplot(2,1,2);
plot(0:tau:vmax, sum(xi0));
grid on;
xlabel('v');
ylabel('Sum_{x} Xi (x,v,0)');
title('Initial distribution across residences times v');  

%Plot density rho and v at t
figure(4);
subplot(2,1,1);
plot((-xmax ):chi:(xmax ), sum(xi,2)/chi);
grid on;
xlabel('x');
str3 = sprintf('Rho (x,%d)',t);
ylabel(str3);
str4 = sprintf('Rho distribution across x at t = %d', t);
title(str4);
subplot(2,1,2);
plot(0:tau:vmax, sum(xi));
grid on;
xlabel('v');
str5 = sprintf('Sum_{x} Xi (x,v,%d)',t);
ylabel(str5);
str6 = sprintf('Distribution across residences times v at t = %d',t); 
title(str6);  




M = zeros(1000,1,1);
for i = 1:1000
X = 0;
T = 0;    
while (T < t)
u1 = rand;
tt = (u1*n*gamma(1-alpha(X)))^(-1/alpha(X));
C = 1 - a(X,T); L = (a(X,T) - chi*b(X,T))/2; R = (a(X,T) + chi*b(X,T))/2;
T = T + tt;
if (T < t)
u2 = rand;
if (u2 < L) 
    X = X - normrnd(0,1/sqrt(n)); 
elseif (u2 > L && u2 < (L+R)) 
    X = X + normrnd(0,1/sqrt(n)); 
end
end
end
M(i) = X;
end
figure(4);
subplot(2,1,1);
hold on;
histogram(M,'Normalization','pdf');
ksdensity(M);
hold off;
    