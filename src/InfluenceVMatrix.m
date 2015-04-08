function [ud,vd,wd,us,vs,ws] = InfluenceVMatrix(Xb,Yb,Zb,pc,cpts,vc,vn,vt,mu,sig,dl,dm,S,ep,epSC,SC)

farfield = 3;

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

% Calculating the induced velocity from the far-field panels.
us_ff = Rx_ff(:,ff).*sig_ff.*S_ff/4/pi./R_ff(:,ff).^3;
vs_ff = Ry_ff(:,ff).*sig_ff.*S_ff/4/pi./R_ff(:,ff).^3;
ws_ff = Rz_ff(:,ff).*sig_ff.*S_ff/4/pi./R_ff(:,ff).^3;

ud_ff = 3/4/pi*p(3,ff).*mu_ff.*S_ff.*Rx_ff(:,ff)./R_ff(:,ff).^5;
vd_ff = 3/4/pi*p(3,ff).*mu_ff.*S_ff.*Ry_ff(:,ff)./R_ff(:,ff).^5;
wd_ff = -(Rx_ff(:,ff).^2 + Ry_ff(:,ff).^2 - 2*p(3,ff).^2).*mu_ff.*S_ff/4/pi./R_ff(:,ff).^5;

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

% Calculating velocity induced at N points by K source panels
At2 = M2.*(M1.*E - M2.*H).*kron(ones(4,1),p(3,nf));
Ab2 = R.*M2.*M2.*kron(ones(4,1),p(3,nf)).^2;

At1 = M2.*(M1.*(Nbar*E) - M2.*(Nbar*H)).*kron(ones(4,1),p(3,nf));
Ab1 = (Nbar*R).*M2.*M2.*kron(ones(4,1),p(3,nf)).^2;

theta2 = atan2(At2,Ab2);
theta1 = atan2(At1,Ab1);

Q = theta2 - theta1;

Us = zeros(4,length(nf),3);

Us(:,:,1) = 1/4/pi*(Dy./D).*log((R + (Nbar*R) - D)./(R + (Nbar*R) + D));
Us(:,:,2) = -1/4/pi*(Dx./D).*log((R + (Nbar*R) - D)./(R + (Nbar*R) + D));
Us(:,:,3) = 1/4/pi*Q;

us_nf = sig_nf.*sum(Us(:,:,1),1);
vs_nf = sig_nf.*sum(Us(:,:,2),1);
ws_nf = sig_nf.*sum(Us(:,:,3),1);

% Enforcing the zero source contribution for TE and wake panels.  The
% problem is the Inf or NaN times zero is NaN.
[col] = find(sig_nf==0);
num = length(col);
for vec = 1:num
    us_nf(col(vec)) = 0;
    vs_nf(col(vec)) = 0;
    ws_nf(col(vec)) = 0;
end

% Calculating velocity induced at N points by K vortex ring panels (equivalent to doublet panels)
R0(:,:,1) = Dx;
R0(:,:,2) = Dy;
R0(:,:,3) = Dz;

R1(:,:,1) = Rx;
R1(:,:,2) = Ry;
R1(:,:,3) = Rz;

R2(:,:,1) = Nbar*Rx;
R2(:,:,2) = Nbar*Ry;
R2(:,:,3) = Nbar*Rz;

c12 = cross(R1,R2,3);
c10 = sqrt(sum((cross(R1,R0,3).^2),3))./sqrt(sum(R0.^2,3));
c20 = sqrt(sum((cross(R2,R0,3).^2),3))./sqrt(sum(R0.^2,3));

Mask0 = ~((sqrt(sum(R1.^2,3)) < ep) | (sqrt(sum(R2.^2,3)) < ep));

MaskSC = (c10 < epSC | c20 < epSC) & Mask0;

K_sc = -1/4/pi/ep^2*1./sum(R0.^2,3).*(sum(R0.*R1,3)./sqrt(sum(R1.^2,3)) - sum(R0.*R2,3)./sqrt(sum(R2.^2,3)));

[row,col] = find(~MaskSC);
num = length(row);
for vec = 1:num
    K_sc(row(vec),col(vec)) = 0;
end

K_SC(:,:,1) = K_sc;
K_SC(:,:,2) = K_sc;
K_SC(:,:,3) = K_sc;

Mask = ~MaskSC & Mask0;
K = -1/4/pi./sum(c12.^2,3).*(sum(R0.*R1,3)./sqrt(sum(R1.^2,3)) - sum(R0.*R2,3)./sqrt(sum(R2.^2,3)));

[row,col] = find(~Mask);
num = length(row);
for vec = 1:num
    K(row(vec),col(vec)) = 0;
end

K3(:,:,1) = K;
K3(:,:,2) = K;
K3(:,:,3) = K;

Ud = SC*K_SC.*c12 + K3.*c12;

ud_nf = mu_nf.*sum(Ud(:,:,1),1);
vd_nf = mu_nf.*sum(Ud(:,:,2),1);
wd_nf = mu_nf.*sum(Ud(:,:,3),1);

% Recombining the near and far-field calculations.

ust = zeros(1,N*Kpan);
vst = zeros(1,N*Kpan);
wst = zeros(1,N*Kpan);
udt = zeros(1,N*Kpan);
vdt = zeros(1,N*Kpan);
wdt = zeros(1,N*Kpan);

ust(:,ff) = us_ff;
ust(:,nf) = us_nf;
vst(:,ff) = vs_ff;
vst(:,nf) = vs_nf;
wst(:,ff) = ws_ff;
wst(:,nf) = ws_nf;

udt(:,ff) = ud_ff;
udt(:,nf) = ud_nf;
vdt(:,ff) = vd_ff;
vdt(:,nf) = vd_nf;
wdt(:,ff) = wd_ff;
wdt(:,nf) = wd_nf;

% Reorganizing matrices.
us = zeros(Kpan,N);
vs = zeros(Kpan,N);
ws = zeros(Kpan,N);
ud = zeros(Kpan,N);
vd = zeros(Kpan,N);
wd = zeros(Kpan,N);

for i = 1:Kpan
    vec = N*(i-1)+1:N*(i-1)+N;
    us(i,:) = ust(1,vec);
    vs(i,:) = vst(1,vec);
    ws(i,:) = wst(1,vec);
    ud(i,:) = udt(1,vec);
    vd(i,:) = vdt(1,vec);
    wd(i,:) = wdt(1,vec);
end

clear ust vst wst udt vdt wdt

% Transforming the induced velocities back into global coordinate systems.
[us,vs,ws] = TransformBack(vc,vn,vt,us,vs,ws);
[ud,vd,wd] = TransformBack(vc,vn,vt,ud,vd,wd);
