% Written by Keith Moored, 12/14/10
%
%
% This code calculates the induced velocties (u,w) at a point (x,z) due to
% a source element of strength sigma located at another point (x0,z0).

function [u_s,w_s,u_d,w_d] = DubV(mu,sigma,x,z,xp,zp,t,n)
N = length(x);

sig = kron(sigma,ones(1,N))';
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

r1 = sqrt(Xp.^2 + Zp.^2);
r2 = sqrt((Xp - X2p).^2 + Zp.^2);
theta1 = atan2(Zp,Xp);
theta2 = atan2(Zp,Xp - X2p);


up_s = sig/4/pi.*log(r1.^2./r2.^2)';
wp_s = sig/2/pi.*(theta2 - theta1)';

[ut_s,wt_s] = TransformMatrixBack(n,t,up_s',wp_s');
u_s = sum(ut_s',2);
w_s = sum(wt_s',2);


up_d = -mu/2/pi.*(Zp./r1.^2 -  Zp./r2.^2)';
wp_d = mu/2/pi.*(Xp./r1.^2 - (Xp - X2p)./r2.^2)';

[ut_d,wt_d] = TransformMatrixBack(n,t,up_d',wp_d');
u_d = sum(ut_d',2);
w_d = sum(wt_d',2);