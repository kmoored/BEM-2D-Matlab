% Written by Keith Moored, 12/14/10
%
%
% This code calculates the induced velocties (u,w) at a point (x,z) due to
% a source element of strength sigma located at another point (x0,z0).

function [dPhi_d] = PhiCD(mu,x,z,xp,zp,t,n)

N = length(x);
mu = kron(mu,ones(1,N))';

x1 = xp(1:end-1);
x2 = xp(2:end);
z1 = zp(1:end-1);
z2 = zp(2:end);

[pp] = TransformMatrix(x,z,t,n,x1,z1);
[p2] = TransformMatrix(x,z,t,n,x2,z2);

Xp = pp(1:2:end,:);
Zp = pp(2:2:end,:);

p2p = pp - p2;
X2p = p2p(1:2:end,:);

theta1 = atan2(Zp,Xp);
theta2 = atan2(Zp,Xp - X2p);

dPhi_d = -mu./2/pi.*(theta2 - theta1)';






