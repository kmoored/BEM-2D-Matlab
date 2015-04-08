function [xt,zt,xb,zb,a,k,epsilon,theta1,theta2,theta,r1,r2] = NACA(c,npts,tmax)

tau = 12.4*(pi/180);        % Trailing edge angle.  This changes with the thickness.
epsilon = 0.075;            % Thickness parameter.

k = 2 - tau/pi;
a = c*(1 + epsilon)^(k - 1)*2^-k;

theta = linspace(0,pi,(npts+1)/2)';
theta1 = atan(a*sin(theta)./(a*cos(theta) - a)) + pi;
theta2 = atan2(a*sin(theta),(a*cos(theta) - epsilon*a));

r1 = sqrt((a*cos(theta) - a).^2 + a^2*sin(theta).^2);
r2 = sqrt((a*cos(theta) - epsilon*a).^2 + a^2*sin(theta).^2);

theta = linspace(pi,0,(npts+1)/2)';
theta1 = atan(a*sin(theta)./(a*cos(theta) - a)) + pi;
theta2 = atan2(a*sin(theta),(a*cos(theta) - epsilon*a));

r1 = sqrt((a*cos(theta) - a).^2 + a^2*sin(theta).^2);
r2 = sqrt((a*cos(theta) - epsilon*a).^2 + a^2*sin(theta).^2);

% tmax = 0.12;

% Stepping through each spanwise position to calculate the positions of the 
% fin neutral plane at the given time step.
xb = linspace(pi,0,(npts+1)/2)';
xt = linspace(0,pi,(npts+1)/2)';

% Creating Ny points along the chord for each spanwise position for a
% total of Nx*Ny points.  The points are spaced more closely at the 
% leading and trailing edges.
x_c = 1/2*(1 - cos(xb));
xb = x_c*c;

% NACA four series shape coefficients.
a0 = 0.2969;
a1 = -0.1260;
a2 = -0.3516;
a3 = 0.2843;
a4 = -0.1015;

% NACA four series shape equation.
zb = -tmax*c/0.2*(a0*sqrt(x_c) + a1*(x_c) + a2*(x_c).^2 + a3*(x_c).^3 + a4*(x_c).^4);
x_c = 1/2*(1 - cos(xt));
xt = x_c*c;

% NACA four series shape equation.
zt = tmax*c/0.2*(a0*sqrt(x_c) + a1*(x_c) + a2*(x_c).^2 + a3*(x_c).^3 + a4*(x_c).^4);


