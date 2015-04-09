function [ U_nPlus, Udot_nPlus, UdotDot_nPlus ] = NEWMARK( beta, ...
    gamma, deltaT, M, K, U_n, Udot_n, UdotDot_n, Fext_nPlus )
%NEWMARK: Solves a dynamic system of equations using the NEWMARK Method.
%This is a transient, implicit method with numerical dissipation in the 
%high frequency domain. This method has first-order accuracy.
%   Solve for the displacement U_nPlus
%   Form the 'A' Matrix
    A = (M + beta*deltaT^2.*K) / (beta*deltaT^2);

%   Form the 'B' Matrix
    B = Fext_nPlus + 1/(beta*deltaT^2).*M*(U_n + ...
        deltaT*Udot_n + deltaT^2*(0.5-beta)*UdotDot_n);

%   Invert the matricies to get U_nPlus
    U_nPlus = A \ B;
  
%   Solve for the acceleration UdotDot_nPlus
    UdotDot_nPlus = (U_nPlus - (U_n + deltaT*Udot_n + ...
                  deltaT^2*(0.5-beta)*UdotDot_n)) / (beta*deltaT^2);

%   Solve for the velocity Udot_nPlus
    Udot_nPlus = (Udot_n + deltaT*(1-gamma)*UdotDot_n) + ...
               gamma*deltaT*UdotDot_nPlus;
end

