% Written by Keith Moored, 12/14/10
%
%
% This code calculates the induced velocties (u,w) at a point (x,z) due to
% a source element of strength sigma located at another point (x0,z0).

function [dPhi_s,dL] = PhiCS(sigma,x,z,x1,z1,x2,z2,t,n)

R = [t; n];

pp = R*[x - x1; z - z1];
p2p = R*[x2 - x1; z2 - z1];

xp = pp(1);
zp = pp(2);

x2p = p2p(1);

r1 = sqrt(xp^2 + zp^2);
r2 = sqrt((xp - x2p)^2 + zp^2);
theta1 = atan2(zp,xp);
theta2 = atan2(zp,xp - x2p);

dPhi_s = sigma/4/pi*(xp*log(r1^2) - (xp - x2p)*log(r2^2) + 2*zp*(theta2 - theta1));
dL = x2p;


