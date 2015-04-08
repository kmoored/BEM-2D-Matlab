% Geometry parameters
tmax      = 0.05*c;                     % Cylinder diameter, meters.
b         = 2*c;                        % Span length, meters.
d_c       = 0;                          % Distance from the ground, chords.
D         = d_c*c;                      % Distance from the ground, meters.

% Kinematic parameters
h_c       = 0.10;                       % Heave-to-chord ratio.
alpha_max = asin(1/2*A_c);              % Maximum pitching angle.
alpha     = 0*(pi/180);                 % Angle of attack.
A_TE      = 1/2*A_c*c;                  % Trailing edge amplitude of motion = 1/2 peak-to-peak amplitude.
phi       = pi;                         % Phase offset of actuation signal.

% Virtual body properties
BL        = 1;                          % Virtual body length, meters.
Beta      = 1;                          % Virtual mass coeficient.
M         = M_star*rho*(c*b)*(A_c*c);   % Virtual body mass, kg.
A_bod     = A_bp*(c*b);                 % Surface area of virtual body, m^2.
Cd_bod    = 0.045;                      % Virtual body drag coefficient.
Ct_prop   = 2.25;                       % Propulsor thrust coefficient.

% Panel parameters and cuttoff radii.
ep        = val*c;
epSC      = ep;
epBod     = ep;
epB       = val*c;
CptInval  = 0.15;
