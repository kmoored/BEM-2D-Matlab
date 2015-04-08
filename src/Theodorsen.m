% This function computes the lift coefficient, Cl = L/(1/2 rho*c*U^2), from
% the Theodorsen model undergoing a heave only motion.  
% The circulatory lift coefficient, Cl_c, and the added mass lift 
% coefficient, Cl_a, are also calculated.  The non-dimensional time, tau, 
% is given as well, which is the time normalized by the period of motion.  
% The lift coefficient is only a function of two non-dimensional 
% parameters: the reduced frequency, k = pi*f*c/U, and the Strouhal number,
% St = 2*f*h_0/U.  The number of time steps is Nstep.

function [Cl,Cl_c,Cl_a,tau] = Theodorsen(k,St,Nstep)

tau = linspace(0,4,Nstep)';

[F,G] = LiftDefFac(k);
C = F + 1i*G;
C_mag = abs(C);
C_phi = angle(C);

Cl_c = -2*pi^2*St*C_mag.*cos(2*pi*tau + C_phi);
Cl_a = pi^2*St*k*sin(2*pi*tau);
Cl = Cl_c + Cl_a;

