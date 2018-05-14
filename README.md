# VaryExp
varyexp.r implements a semi-Markov numeric scheme to approximate densities of CTRW (continuous time random walk) limits. The densities of these CTRW limits solve a variety of Fokker-Planck equations, Fractional Fokker-Plack equations and Variable order Fractional Fokker-Planck equations and so provide a model for various physical systems exhibiting diffusion or anomalous diffusion.

## Example
Below is a plot of the density P(x,t) of anomalous diffusion with a varying exponent. The tail function of the waiting time distribution is given by \Psi(x,t) = t^(-\alpha (x))/\Gamma(1-\alpha(x)) while the varying exponent is a logistic function \alpha (x) = 0.5/(1+e^(-x)) + 0.1

![CTRW Anomalous diffusion example](/README-example.png)
