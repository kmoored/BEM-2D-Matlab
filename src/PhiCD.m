% Written by Keith Moored, 12/14/10
%
%
% This code calculates the induced velocties (u,w) at a point (x,z) due to
% a source element of strength sigma located at another point (x0,z0).

function [dPhi_d] = PhiCD(mu,x,z,x1,z1,x2,z2,t,n)

R = [t; n];

pp = R*[x - x1; z - z1];
p2p = R*[x2 - x1; z2 - z1];

xp = pp(1);
zp = pp(2);

x2p = p2p(1);

theta1 = atan2(zp,xp);
theta2 = atan2(zp,xp-x2p);

dPhi_d = -mu/2/pi*(theta2 - theta1);
