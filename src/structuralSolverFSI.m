function [ U_nPlus, Udot_nPlus, UdotDot_nPlus ] = structuralSolverFSI( ...
    alpha, beta, gamma, Nelements, fixedNodes, l_0, rho, E, A, I, Fload,  ...
    U_n, Udot_n, endTime, fracDeltaT, method)
%structuralSolverFSI: Solves a dynamic system of equations using either the
% Trapdezoidal Rule, NEWMARK method, or Hilber-Hughes-Taylor (HHT) Method. 
% This is a transient, implicit method and returns the final displacement,
% velocity, and acceleration matrix.
%   Define the structural time-step size    
    deltaT    = fracDeltaT*endTime;

%   Matrix initialization
    M = zeros(3*(Nelements+1));
    K = zeros(3*(Nelements+1));

%   Assemble global mass and stiffness matricies
    for i = 1:Nelements
        k_e = zeros(6,6);
        m_e = zeros(6,6); 
        l_e = zeros(6,2*(Nelements+1));
        [ k_e ] = elementStiffnessMatrix( E, I(i), A(i), l_0 );
        [ m_e ] = elementMassMatrix( rho, A(i), l_0, 'consistent' );
        [ l_e ] = elementConnectivityMatrix( Nelements, i, -U_n(3*i,1));
        M = M + l_e' * m_e * l_e;
        K = K + l_e' * k_e * l_e;
    end
    
%   Set zero displacement constraint
    temp = 3*fixedNodes+1;
    M = M(temp:end,temp:end);
    K = K(temp:end,temp:end);
    Fload = Fload(temp:end,1);
    U_n = U_n(temp:end,1);
    Udot_n = Udot_n(temp:end,1);

%   Solve for the inital acceleration matrix
    Fext_n = Fload; 
    UdotDot_n = M \ (Fext_n - K*U_n);
    
%   March through time until the total simulated time has elapsed
    for i = deltaT:deltaT:endTime
%       Solve for the time-step using the specified method
        Fext_nPlus = Fext_n;
        if strcmp(method,'HHT')
            [ U_nPlus, Udot_nPlus, UdotDot_nPlus ] = HHT( alpha, beta, ...
                gamma, deltaT, M, K, U_n, Udot_n, UdotDot_n, Fext_n, Fext_nPlus );
        elseif strcmp(method,'NEWMARK')
            [ U_nPlus, Udot_nPlus, UdotDot_nPlus ] = NEWMARK( beta, gamma, ...
                deltaT, M, K, U_n, Udot_n, UdotDot_n, Fext_nPlus );
        elseif strcmp(method,'TRAPEZOIDAL')
            [ U_nPlus, Udot_nPlus, UdotDot_nPlus ] = trapezoidalRuleFE( U_n, ...
                Udot_n, UdotDot_n, deltaT, M, K, Fext_nPlus );
        else
            disp('ERROR! Invalid integrations scheme! Valid schemes are:')
            disp('    * HHT')
            disp('    * NEWMARK')
            disp('    * TRAPEZOIDAL')
        end

%       Update displacement, velocity, and acceleration if this is not the
%       last time-step for the inital conditions of the next time-step
        if i ~= endTime
            U_n = U_nPlus;
            Udot_n = Udot_nPlus;
            UdotDot_n = UdotDot_nPlus;
        end
    end
end

