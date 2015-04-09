function [ U_nPlus, Udot_nPlus, UdotDot_nPlus ] = trapezoidalRuleFE( U_n, Udot_n, UdotDot_n, deltaT, M, K, Fext_nPlus )
%trapezoidalRuleFE: Solves for the system dynamics using the trapezoidal
%rule.
    U_nPlus = (K + (2/deltaT)^2.*M) \ (Fext_nPlus + M*((2/deltaT)^2.*U_n...
        + (4/deltaT).*Udot_n + UdotDot_n));
    
    Udot_nPlus = 2 * (U_nPlus - U_n) / deltaT - Udot_n;
    
    UdotDot_nPlus = 2 * (Udot_nPlus - Udot_n) / deltaT - UdotDot_n;
end

