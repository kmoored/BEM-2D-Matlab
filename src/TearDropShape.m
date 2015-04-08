function [xt,zt,xb,zb] = TearDropShape(c,npts,D)

% Stepping through each spanwise position to calculate the positions of the 
% fin neutral plane at the given time step.
xb = linspace(pi,0,(npts+1)/2)';
xt = linspace(0,pi,(npts+1)/2)';

% Slopes and intersects for the line segments
m = -D/2/(c - D/2);
b = D/2 + D^2/4/(c - D/2);

% Tear drop shape equation.
x_c = 1/2*(1 - cos(xb));
xb = x_c*c;
xb1 = xb(xb <= D/2);
xb2 = xb(xb > D/2);

zb2 = -m*xb2 - b;
zb1 = -sqrt((D/2)^2 - (xb1 - D/2).^2);
zb = [zb2; zb1];

% Tear drop shape equation.
x_c = 1/2*(1 - cos(xt));
xt = x_c*c;
xt1 = xt(xt <= D/2);
xt2 = xt(xt > D/2);

zt1 = sqrt((D/2)^2 - (xt1 - D/2).^2);
zt2 = m*xt2 + b;
zt = [zt1; zt2];

