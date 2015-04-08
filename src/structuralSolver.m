close all
clear all
clc

Nelements = 100;
Nnodes    = Nelements + 1;
Lbeam     = 1.0;
rho       = 1700;
E         = 25.5e6;
%rho       = 7750;
%E         = 2.1e11;
height    = 0.1;
width     = 1;
endTime   = 1/30;
deltaT    = 0.01*endTime;

% NEWMARK/HHT constants
alpha = 0.01;
beta = (1+alpha)^2/4;
gamma = 0.5+alpha;

% Calculated inputs
l = Lbeam / Nelements;
%A = height * width;
%I = width * height^3 / 12;
load('A.mat');
load('I.mat');
nodes_0(:,1) = 0:l:Lbeam;
nodes_0(:,2) = 0;
nodes = nodes_0;

% Matrix initialization
M = zeros(3*(Nelements+1));
K = zeros(3*(Nelements+1));
Fload = zeros(3*(Nelements+1),1);
Fext_n = zeros(3*(Nelements+1),1);
U_n = zeros(3*(Nelements+1),1);
Udot_n = zeros(3*(Nelements+1),1);
UdotDot_n = zeros(3*(Nelements+1),1);

% Define load condition
load('Fload.mat');

% Assemble global mass and stiffness matricies
for i = 1:Nelements
    k_e = zeros(4,4);
    m_e = zeros(4,4); 
    l_e = zeros(4,2*(Nelements+1));
    [ k_e ] = elementStiffnessMatrix( E, I(i), A(i), l );
    [ m_e ] = elementMassMatrix( rho, A(i), l, 'consistent' );
    [ l_e ] = elementConnectivityMatrix( Nelements, i, -U_n(3*i,1) );
    M = M + l_e' * m_e * l_e;
    K = K + l_e' * k_e * l_e;
end
% Set zero displacement constraint
fixedNodes = 31;
M = M(3*fixedNodes+1:end,3*fixedNodes+1:end);
K = K(3*fixedNodes+1:end,3*fixedNodes+1:end);
Fload = Fload(3*fixedNodes+1:end,1);
U_n = U_n(3*fixedNodes+1:end,1);
Udot_n = Udot_n(3*fixedNodes+1:end,1);
UdotDot_n = UdotDot_n(3*fixedNodes+1:end,1);

% Define inital displacement and velcotiy


% Solve for the inital acceleration matrix
%[ U_nPlus ] = trapezoidalRuleFE( U_n, Udot_n, UdotDot_n, deltaT, M, K, Fext_nPlus );
Fext_n = Fload; 
UdotDot_n = M \ (Fext_n - K*U_n);
for i = deltaT:deltaT:endTime
    % Solve for the timestep
    Fext_nPlus = Fext_n;
    [ U_nPlus, Udot_nPlus, UdotDot_nPlus ] = HHT( alpha, beta, ...
        gamma, deltaT, M, K, U_n, Udot_n, UdotDot_n, Fext_n, Fext_nPlus );
    %[ U_nPlus, Udot_nPlus, UdotDot_nPlus ] = NEWMARK( beta, ...
    %gamma, deltaT, M, K, U_n, Udot_n, UdotDot_n, Fext_nPlus );
    %[ U_nPlus, Udot_nPlus, UdotDot_nPlus ] = trapezoidalRuleFE( U_n, ...
    %    Udot_n, UdotDot_n, deltaT, M, K, Fext_nPlus );
    if i ~= endTime
        U_n = U_nPlus;
        Udot_n = Udot_nPlus;
        UdotDot_n = UdotDot_nPlus;
    end
end

% Update nodes
nodes(fixedNodes+1:end,2) = nodes(fixedNodes+1:end,2) + U_nPlus(2:3:end-1,1);
nodes(fixedNodes+1:end,1) = nodes(fixedNodes+1:end,1) + nodes(fixedNodes+1:end,2).*sin(-U_nPlus(3:3:end,1) + U_nPlus(1:3:end-2,1));

% Plot the initial and final results
plot(nodes_0(:,1),nodes_0(:,2),nodes(:,1),nodes(:,2))
legend('Original','Final')
