function [xp,zp,Vp] = Kinematics_HeavePitch2D(xp,zp,alpha_max,h_c,f,t,phi)

Npanels = length(xp)-1;
Nx = Npanels/2 + 1;

xpt = xp(Nx:end);
zpb = zp(1:Nx);
zpt = zp(Nx:end);

%c = xpt(end) - xpt(1);
c = 1;

v1 = zeros(2,Nx);
v2 = zeros(2,Nx);
Vp = zeros(2,2*Nx - 1);

% Calculating time before and after real-time for surface velocity 
% calcualtion.  The change in time is very small to accurately calculate
% the surface velocity.
Tstep = 1e-5;
Dstep = 1e-5;
t_plus = t + Tstep;
t_minus = t - Tstep;

% Calculating the normal vectors at each point on the neutral plane.


% Calculating position of the nuetral plane heaving and pitching about
% the leading edge.
x_n = (xpt - xpt(1))*cos(alpha_max*sin(2*pi*f*t + phi));
z_n = h_c*c*sin(2*pi*f*t) + (xpt - xpt(1))*sin(alpha_max*sin(2*pi*f*t + phi));

xplus = xpt + Dstep;
xminus = xpt - Dstep;

% Calculating positions of four points around position of neutral
% plane.
v1(1,:) = (xplus - xpt(1))*cos(alpha_max*sin(2*pi*f*t + phi));
v1(2,:) = h_c*c*sin(2*pi*f*t) + (xplus - xpt(1))*sin(alpha_max*sin(2*pi*f*t + phi));

v2(1,:) = (xminus - xpt(1))*cos(alpha_max*sin(2*pi*f*t + phi));
v2(2,:) = h_c*c*sin(2*pi*f*t) + (xminus - xpt(1))*sin(alpha_max*sin(2*pi*f*t + phi));

% Vector in y-direction is defined as v1-v2 while vector in
% x-direction is defined as v3-v4.
vt = v1-v2;
vn = [-vt(2,:); vt(1,:)];

for i = 1:Nx
    vn(:,i) = vn(:,i)/norm(vn(:,i));
end





% Calculating the position of a point at time t_plus.
x_np = (xpt - xpt(1))*cos(alpha_max*sin(2*pi*f*t_plus + phi));
z_np = h_c*c*sin(2*pi*f*t_plus) + (xpt - xpt(1))*sin(alpha_max*sin(2*pi*f*t_plus + phi));

% Calculating positions of four points around position of neutral
% plane for a time t_plus
v1(1,:) = (xplus - xpt(1))*cos(alpha_max*sin(2*pi*f*t_plus + phi));
v1(2,:) = h_c*c*sin(2*pi*f*t_plus) + (xplus - xpt(1))*sin(alpha_max*sin(2*pi*f*t_plus + phi));

v2(1,:) = (xminus - xpt(1))*cos(alpha_max*sin(2*pi*f*t_plus + phi));
v2(2,:) = h_c*c*sin(2*pi*f*t_plus) + (xminus - xpt(1))*sin(alpha_max*sin(2*pi*f*t_plus + phi));

% Vector in y-direction is defined as v1-v2 while vector in
% x-direction is defined as v3-v4.
vt = v1-v2;
vnp = [-vt(2,:); vt(1,:)];

for i = 1:Nx
    vnp(:,i) = vnp(:,i)/norm(vnp(:,i));
end




% Calculating the position of a point at time t_minus.
x_nm = (xpt - xpt(1))*cos(alpha_max*sin(2*pi*f*t_minus + phi));
z_nm = h_c*c*sin(2*pi*f*t_minus) + (xpt - xpt(1))*sin(alpha_max*sin(2*pi*f*t_minus + phi));


% Calculating positions of four points around position of neutral
% plane for a time t_minus
v1(1,:) = (xplus - xpt(1))*cos(alpha_max*sin(2*pi*f*t_minus + phi));
v1(2,:) = h_c*c*sin(2*pi*f*t_minus) + (xplus - xpt(1))*sin(alpha_max*sin(2*pi*f*t_minus + phi));

v2(1,:) = (xminus - xpt(1))*cos(alpha_max*sin(2*pi*f*t_minus + phi));
v2(2,:) = h_c*c*sin(2*pi*f*t_minus) + (xminus - xpt(1))*sin(alpha_max*sin(2*pi*f*t_minus + phi));

% Vector in y-direction is defined as v1-v2 while vector in
% x-direction is defined as v3-v4.
vt = v1-v2;
vnm = [-vt(2,:); vt(1,:)];

for i = 1:Nx
    vnm(:,i) = vnm(:,i)/norm(vnm(:,i));
end



% Creating a three dimensional array of the coordinates of the upper
% and lower surfaces.  The row dimension goes from the leading-edge to
% the trailing-edge in the chord direction.  The column dimension goes
% from the root to the tip is the spanwise direction. The slice (1, 2, 
% or 3) corresponds to the x, y or z coordinate.  Each array (u-upper,
% l-lower) is Ny by Nx by 3 in size.

xpt = x_n' + vn(1,:)*diag(zpt);
zpt = z_n' + vn(2,:)*diag(zpt);

xpb = x_n(end:-1:1)' + vn(1,end:-1:1)*diag(zpb);
zpb = z_n(end:-1:1)' + vn(2,end:-1:1)*diag(zpb);

xpt_p = x_np' + vnp(1,:)*diag(zpt);
zpt_p = z_np' + vnp(2,:)*diag(zpt);

xpb_p = x_np(end:-1:1)' + vnp(1,end:-1:1)*diag(zpb);
zpb_p = z_np(end:-1:1)' + vnp(2,end:-1:1)*diag(zpb);

xpt_m = x_nm' + vnm(1,:)*diag(zpt);
zpt_m = z_nm' + vnm(2,:)*diag(zpt);

xpb_m = x_nm(end:-1:1)' + vnm(1,end:-1:1)*diag(zpb);
zpb_m = z_nm(end:-1:1)' + vnm(2,end:-1:1)*diag(zpb);


% If the tip goes to a point the measured chord will go to zero.  When
% this happens the y_c coordinates become NaNs, generating NaNs for all
% of the y and z coordinates of the upper and lower surfaces.  Instead
% this condition replaces those NaNs with the appropriate values.  From
% a computational point of view, all of these points coincide and could
% be represented by a single separate value instead of Ny copies of the
% same coordinate.  However, it is simpler to keep all of the
% coordinates in a single array.


xp = [xpb xpt(2:end)];
zp = [zpb zpt(2:end)];

xp_p = [xpb_p xpt_p(2:end)];
zp_p = [zpb_p zpt_p(2:end)];

xp_m = [xpb_m xpt_m(2:end)];
zp_m = [zpb_m zpt_m(2:end)];


% Calculating velocity of the collocation point
Vp(1,:) = (xp_p - xp_m)/Tstep/2;
Vp(2,:) = (zp_p - zp_m)/Tstep/2;

xp = xp';
zp = zp';
Vp = Vp';
