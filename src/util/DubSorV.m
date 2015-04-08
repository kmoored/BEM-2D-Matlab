% Written by Keith Moored, 12/14/10
%
%
% This code calculates the induced velocties (u,w) at a point (x,z) due to
% a source element of strength sigma located at another point (x0,z0).

function [u_s,w_s,u_d,w_d] = DubSorV(mu,sigma,x,z,xp,zp,t,n,epSC,SC)

N = length(x);                  % N points
K = length(xp) - 1;             % K panels

sig = repmat(sigma',1,N);
mu = repmat(mu',1,N);

x1 = xp(1:end-1);
x2 = xp(2:end);
z1 = zp(1:end-1);
z2 = zp(2:end);

% Creating local cooridnate systems with the origin at p1 = (x1,z1) of each
% panel.

% Transforming the point P = (x,z) into the local coordinate system.
[P_p1] = TransformMatrix(x1,z1,t,n,x,z);

X = P_p1(1:2:end,:);
Z = P_p1(2:2:end,:);

% Transforming the point P2 = (x2,z2) into the local coordinate sytem.
[P2_p1] = TransformMatrix(x1,z1,t,n,x2,z2);

x2 = diag(P2_p1(1:2:end,:));
X2 = repmat(x2,1,N);

% Calculating vector magnitudes and angles from P1 -> P and P2 -> P. 
r1 = sqrt(X.^2 + Z.^2);
r2 = sqrt((X - X2).^2 + Z.^2);
theta1 = atan2(Z,X);
theta2 = atan2(Z,X - X2);

% Calculating the induced velocity from the source panels
if sum(sigma ~= 0) > 0
    up_s = sig/4/pi.*log(r1.^2./r2.^2);
    wp_s = sig/2/pi.*(theta2 - theta1);

    [u_s,w_s] = TransformMatrixBack(t,n,up_s,wp_s);
else
    u_s = zeros(1,N);
    w_s = zeros(1,N);
end



% Calculating the induced velocity from the doublet panels
V1u = Z./r1.^2;
V2u = Z./r2.^2;

V1w = X./r1.^2;
V2w = (X - X2)./r2.^2;

% Creating a zero velocity core around the vortex endpoints to avoid
% velocity singularities.
Mask1 = (r1 > epSC); 
if ~isempty(Mask1)
    V1u = Mask1.*V1u;
    V1w = Mask1.*V1w;
    V1u(isnan(V1u)) = 0;
    V1w(isnan(V1w)) = 0;
end

Mask2 = (r2 > epSC);
if ~isempty(Mask2)
    V2u = Mask2.*V2u;
    V2w = Mask2.*V2w;
    V2u(isnan(V2u)) = 0;
    V2w(isnan(V2w)) = 0;
end

up1_SC = -mu/2/pi/epSC^2.*r1.*sin(theta1);
wp1_SC = mu/2/pi/epSC^2.*r1.*cos(theta1);

up2_SC = mu/2/pi/epSC^2.*r2.*sin(theta2);
wp2_SC = -mu/2/pi/epSC^2.*r2.*cos(theta2);

up_d = -mu/2/pi.*(V1u - V2u);
wp_d = mu/2/pi.*(V1w - V2w);

up_d = up_d + SC*up1_SC.*(~Mask1) + SC*up2_SC.*(~Mask2);
wp_d = wp_d + SC*wp1_SC.*(~Mask1) + SC*wp2_SC.*(~Mask2);

% Transform back from local coordinates to global coordinates and summing
% contributiions from all panels at a single point of interest.
[u_d,w_d] = TransformMatrixBack(t,n,up_d,wp_d);

