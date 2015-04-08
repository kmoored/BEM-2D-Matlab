% Written by Keith Moored, 12/14/10
%
%
% This code calculates the induced velocties (u,w) at a point (x,z) due to
% a source element of strength sigma located at another point (x0,z0).

function [dPhi_s,dPhi_d,dL] = Phi(mu,sigma,x,z,xp,zp,t,n)

% Points of interest P = (x,z)
% Panel endpoints Pp = (xp,zp)
farfield = 5;
ep = 1e-10;

N = length(x);              % N points
K = length(xp) - 1;         % K panels

sig = repmat(sigma',1,N);
mu = repmat(mu',1,N);

x1 = xp(1:end-1);
x2 = xp(2:end);
z1 = zp(1:end-1);
z2 = zp(2:end);

% Creating local cooridnate systems with the origin at p1 = (x1,z1) of each panel.

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

r_ff = r1/2+r2/2;

r1 = reshape(r1,[1 N*K]);
r2 = reshape(r2,[1 N*K]);
r_ff = reshape(r_ff,[1 N*K]);
X = reshape(X,[1 N*K]);
Z = reshape(Z,[1 N*K]);
X2 = reshape(X2,[1 N*K]);
mu = reshape(mu,[1 N*K]);
sig = reshape(sig,[1 N*K]);

% Determining which points lie in the near-field or far-field regions in
% reference to a given panel.
ff = find(r_ff > farfield*abs(X2));
nf = find(r_ff <= farfield*abs(X2));
 
% Splitting panel properties into near and far-field components.
mu_ff = mu(ff);
mu_nf = mu(nf);
sig_ff = sig(ff);
sig_nf = sig(nf);
X_nf = X(nf);
Z_ff = Z(ff);
Z_nf = Z(nf);
X2_nf = X2(nf);

% Far-field Calculation
dPhi_s_ff = sig_ff.*abs(X2(ff))/2/pi.*log(r_ff(ff));
dPhi_d_ff = -mu_ff.*X2(ff)/2/pi.*Z_ff./(r_ff(ff).^2);

% Near-field Calculation
% theta1 = atan2(Z_nf,X_nf)
% theta2 = atan2(Z_nf,(X_nf - X2_nf))

theta1 = atan2(Z_nf,X_nf);
theta2 = atan2(Z_nf,(X_nf - X2_nf));

Mask = ~(abs(Z_nf) < ep);

if sum(sigma ~= 0) > 0
    dPhi_s_nf = sig_nf./4/pi.*Mask.*(2*X_nf.*log(r1(nf)) - 2*(X_nf - X2_nf).*log(r2(nf)) - 2*X2_nf + 2*Z_nf.*(theta2 - theta1));
else
    dPhi_s_nf = 0;
end

dPhi_d_nf = -mu_nf./2/pi.*Mask.*(theta2 - theta1);

% Recombining the near and far-field calculations.
dPhi_d = zeros(1,K*N);
dPhi_s = zeros(1,K*N);

dPhi_d(1,ff) = dPhi_d_ff;
dPhi_d(1,nf) = dPhi_d_nf;
dPhi_s(1,ff) = dPhi_s_ff;
dPhi_s(1,nf) = dPhi_s_nf;

dPhi_d = reshape(dPhi_d,[K N]);
dPhi_s = reshape(dPhi_s,[K N]);

dL = abs(x2);












