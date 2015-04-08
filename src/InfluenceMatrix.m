% Written by Keith Moored, 8/23/11.
%
%
% This code calculates the influence coefficient of K doublet and source 
% panels at N points.

function [Ac,Bc] = InfluenceMatrix(Xb,Yb,Zb,pc,cpts,vc,vn,vt,mu,sig,dl,dm,S)

farfield = 3;
ep = 1e-10;

[~,Kpan] = size(cpts);
[~,N] = size(Xb);

p = zeros(3,Kpan*N);

C1 = [-1  1  0  0;...
       0 -1  1  0;...
       0  0 -1  1;...
       1  0  0 -1];   
C = kron(C1,eye(3));

[Xp] = TransformMatrix(cpts,vc,vn,vt,Xb,Yb,Zb);

xp = Xp(1:3:3*Kpan,:);
yp = Xp(2:3:3*Kpan,:);
zp = Xp(3:3:3*Kpan,:);

for i = 1:Kpan
    vec = N*(i-1)+1:N*(i-1)+N;
    p(1,vec) = xp(i,:);
    p(2,vec) = yp(i,:);
    p(3,vec) = zp(i,:);
end

% Transforming the four corner points of the panel into the panel frame of
% reference.
Nbarpri = sparse(eye(3*Kpan));
Nbar1 = sparse([Nbarpri(1:3:3*Kpan,:); Nbarpri(2:3:3*Kpan,:); Nbarpri(3:3:3*Kpan,:)]');

clear Nbarpri

[Pt] = TransformMatrixPc(cpts,vc,vn,vt,pc(1,:),pc(2,:),pc(3,:));
xpt = Pt(1:3:3*Kpan,:);
ypt = Pt(2:3:3*Kpan,:);
zpt = Pt(3:3:3*Kpan,:);
pct1 = Nbar1*[xpt; ypt; zpt];

[Pt] = TransformMatrixPc(cpts,vc,vn,vt,pc(4,:),pc(5,:),pc(6,:));
xpt = Pt(1:3:3*Kpan,:);
ypt = Pt(2:3:3*Kpan,:);
zpt = Pt(3:3:3*Kpan,:);
pct2 = Nbar1*[xpt; ypt; zpt];

[Pt] = TransformMatrixPc(cpts,vc,vn,vt,pc(7,:),pc(8,:),pc(9,:));
xpt = Pt(1:3:3*Kpan,:);
ypt = Pt(2:3:3*Kpan,:);
zpt = Pt(3:3:3*Kpan,:);
pct3 = Nbar1*[xpt; ypt; zpt];

[Pt] = TransformMatrixPc(cpts,vc,vn,vt,pc(10,:),pc(11,:),pc(12,:));
xpt = Pt(1:3:3*Kpan,:);
ypt = Pt(2:3:3*Kpan,:);
zpt = Pt(3:3:3*Kpan,:);
pct4 = Nbar1*[xpt; ypt; zpt];

clear Nbar1

Nbarpri = sparse(eye(4*Kpan));
Nbar = sparse([Nbarpri(1:4:4*Kpan,:); Nbarpri(2:4:4*Kpan,:); Nbarpri(3:4:4*Kpan,:); Nbarpri(4:4:4*Kpan,:)]');

clear Nbarpri

pct = sparse(kron(Nbar,eye(3)))*[pct1; pct2; pct3; pct4];

clear Nbar

% Eliminating the z-coordinates of the panel corners in the panel
% coordinate system.  They are not needed for the influence calculations.
pct = sparse(kron(sparse(eye(4*Kpan)),[1 0 0; 0 1 0; 0 0 0]))*pct;

% Calculating the distance from a point of interest to the panel corner 
% points. 
Pctt = kron(ones(1,N),pct);

clear pct

Pct = zeros(12,Kpan*N);
for i = 1:Kpan
    Pct(:,N*(i-1)+1:N*(i-1)+N)  = Pctt(12*(i-1)+1:12*(i-1)+12,:); 
end

clear Pctt

% Calculating the distance from the collocation points to the points of
% interest.

Rx_ff = p(1,:);
Ry_ff = p(2,:);
Rz_ff = p(3,:);
R_ff = sqrt(Rx_ff.^2 + Ry_ff.^2+ Rz_ff.^2);

% Determining which points lie in the near-field or far-field regions in
% reference to a given panel.
ff = find(R_ff > farfield*(kron(dm,ones(1,N)) + kron(dl,ones(1,N))));
nf = find(R_ff <= farfield*(kron(dm,ones(1,N)) + kron(dl,ones(1,N))));

mut = kron(mu',ones(1,N));
sigt = kron(sig',ones(1,N));
St = kron(S,ones(1,N));

% Splitting panel properties into near and far-field components.
mu_ff = mut(:,ff);
mu_nf = mut(:,nf);
sig_ff = sigt(:,ff);
sig_nf = sigt(:,nf);
S_ff = St(:,ff);

% Calculating the velocity potential from the far-field panels.
Ac_ff = -1/4/pi*mu_ff.*S_ff.*p(3,ff)./R_ff(:,ff).^3;
Bc_ff = -1/4/pi*sig_ff.*S_ff./R_ff(:,ff);

% Near-field calculations.
P = kron(ones(4,1),p);
P_nf = P(:,nf);

clear P

Pct_nf = Pct(:,nf);

Rvec = P_nf - Pct_nf;
Dvec = C*Pct_nf;
Evec = kron(eye(4),[1 0 0; 0 0 0; 0 0 1])*Rvec;

% Block diagonalizing Rvec, Dvec and Evec.
Rx = Rvec(1:3:12,:);
Ry = Rvec(2:3:12,:);
Rz = Rvec(3:3:12,:);

R = sqrt(Rx.^2 + Ry.^2 + Rz.^2);

Dx = Dvec(1:3:12,:);
Dy = Dvec(2:3:12,:);
Dz = Dvec(3:3:12,:);

D = sqrt(Dx.^2 + Dy.^2 + Dz.^2);

Ex = Evec(1:3:12,:);
Ey = Evec(2:3:12,:);
Ez = Evec(3:3:12,:);

E = Ex.^2 + Ey.^2 + Ez.^2;

H = Rx.*Ry;

M1 = Dy;
M2 = Dx;

Nbarpri = eye(4);
Nbar = zeros(4);
Nbar(1,:) = Nbarpri(2,:);
Nbar(2,:) = Nbarpri(3,:);
Nbar(3,:) = Nbarpri(4,:);
Nbar(4,:) = Nbarpri(1,:);

At2 = M2.*(M1.*E - M2.*H).*kron(ones(4,1),p(3,nf));
Ab2 = R.*M2.*M2.*kron(ones(4,1),p(3,nf)).^2;

At1 = M2.*(M1.*(Nbar*E) - M2.*(Nbar*H)).*kron(ones(4,1),p(3,nf));
Ab1 = (Nbar*R).*M2.*M2.*kron(ones(4,1),p(3,nf)).^2;

theta2 = atan2(At2,Ab2);
theta1 = atan2(At1,Ab1);

Q = theta2 - theta1;

Beta = Rx.*Dy - Ry.*Dx;
Snew = log((R + (Nbar*R) + D)./(R + (Nbar*R) - D));
S = Beta.*Snew./D;

Mask = ~(D < ep);

% Calculating the influence coefficients.

Actt = -1/4/pi*sum(Mask.*Q,1);

% Mask2 = ~(abs(p(3,nf)) < ep);
% Actt = Mask2.*Actt;

Bctt = -1/4/pi*sum(Mask.*S,1) - p(3,nf).*Actt;

% Calculating Velocity Potential at N point due to K panels.

Ac_nf = mu_nf.*Actt;
Bc_nf = sig_nf.*Bctt;

% Recombining the near and far-field calculations.

Act = zeros(1,N*Kpan);
Bct = zeros(1,N*Kpan);

Act(:,ff) = Ac_ff;
Act(:,nf) = Ac_nf;
Bct(:,ff) = Bc_ff;
Bct(:,nf) = Bc_nf;

% Reorganizing matrices.
Ac = zeros(Kpan,N);
Bc = zeros(Kpan,N);

for i = 1:Kpan
    vec = N*(i-1)+1:N*(i-1)+N;
    Ac(i,:) = Act(1,vec);
    Bc(i,:) = Bct(1,vec);
end

