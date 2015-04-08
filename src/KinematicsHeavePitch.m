function [pc,Vc,cpts,vc,vn,vt,S,P,C1,dl,dm,dm_hat,tp] = KinematicsHeavePitch2D(thick,c,Ny,Nx,yf,xf,t,f,phi,LE,alpha_max,h_c,c_r)

Npanels = (Ny - 1)*(Nx - 1); 

C1 = [-1  1  0  0;...
       0 -1  1  0;...
       0  0 -1  1;...
       1  0  0 -1;...
      -1  0  1  0;...
       0  1  0 -1];
Coll = kron([1 1 1 1],eye(3));

% Pre-allocating memory for increased computational speed.
P_n = zeros(Nx,Ny,3);
P_np = zeros(Nx,Ny,3);
P_nm = zeros(Nx,Ny,3);
vn_n = zeros(Nx,Ny,3);
vn_np = zeros(Nx,Ny,3);
vn_nm = zeros(Nx,Ny,3);
P_u = zeros(Nx,Ny,3);
P_l = zeros(Nx,Ny,3);
P_up = zeros(Nx,Ny,3);
P_lp = zeros(Nx,Ny,3);
P_um = zeros(Nx,Ny,3);
P_lm = zeros(Nx,Ny,3);
tp = zeros(Npanels,1);

% Calculating time before and after real-time for surface velocity 
% calcualtion.  The change in time is very small to accurately calculate
% the surface velocity.
Tstep = 1e-5;
Dstep = 1e-5;
t_plus = t + Tstep;
t_minus = t - Tstep;

% Calculating the normal vectors at each point on the neutral plane.
for i = 1:Nx
    for j = 1:Ny
        % Calculating position of point.
        xM = xf(i,j).*cos(alpha_max*sin(2*pi*f*t + phi));
        yM = yf(j);
        zM = h_c*c_r*sin(2*pi*f*t) + xf(i,j).*sin(alpha_max*sin(2*pi*f*t + phi)); 
        
        P_n(i,j,1) = xM;
        P_n(i,j,2) = yM;
        P_n(i,j,3) = zM;
        
        xplus = xf(i,j) + Dstep; 
        xminus = xf(i,j) - Dstep;
        yplus = yf(j) + Dstep;
        yminus = yf(j) - Dstep;
        
        % Calculating positions of four points around position of neutral
        % plane.
        v1(1,1) = xf(i,j).*cos(alpha_max*sin(2*pi*f*t + phi));
        v1(2,1) = yplus;
        v1(3,1) = h_c*c_r*sin(2*pi*f*t) + xf(i,j).*sin(alpha_max*sin(2*pi*f*t + phi)); 
        
        v2(1,1) = xf(i,j).*cos(alpha_max*sin(2*pi*f*t + phi));
        v2(2,1) = yminus;
        v2(3,1) = h_c*c_r*sin(2*pi*f*t) + xf(i,j).*sin(alpha_max*sin(2*pi*f*t + phi)); 
        
        v3(1,1) = xplus.*cos(alpha_max*sin(2*pi*f*t + phi));
        v3(2,1) = yf(j);
        v3(3,1) = h_c*c_r*sin(2*pi*f*t) + xplus.*sin(alpha_max*sin(2*pi*f*t + phi)); 
        
        v4(1,1) = xminus.*cos(alpha_max*sin(2*pi*f*t + phi));
        v4(2,1) = yf(j);
        v4(3,1) = h_c*c_r*sin(2*pi*f*t) + xminus.*sin(alpha_max*sin(2*pi*f*t + phi)); 
        
        % Vector in y-direction is defined as v1-v2 while vector in
        % x-direction is defined as v3-v4.
        vy = v1-v2;
        vx = v3-v4;
        vn = cross(vx,vy);
        
        vn_n(i,j,1) = vn(1)/norm(vn);
        vn_n(i,j,2) = vn(2)/norm(vn);
        vn_n(i,j,3) = vn(3)/norm(vn);
        
        
        
        % Calculating the position of a point at time t_plus.
        xM_p = xf(i,j).*cos(alpha_max*sin(2*pi*f*t_plus + phi));
        yM_p = yf(j);
        zM_p = h_c*c_r*sin(2*pi*f*t_plus) + xf(i,j).*sin(alpha_max*sin(2*pi*f*t_plus + phi)); 
        
        P_np(i,j,1) = xM_p;
        P_np(i,j,2) = yM_p;
        P_np(i,j,3) = zM_p;
        
        % Calculating positions of four points around position of neutral
        % plane for a time t_plus.
        v1(1,1) = xf(i,j).*cos(alpha_max*sin(2*pi*f*t_plus + phi));
        v1(2,1) = yplus;
        v1(3,1) = h_c*c_r*sin(2*pi*f*t_plus) + xf(i,j).*sin(alpha_max*sin(2*pi*f*t_plus + phi)); 
        
        v2(1,1) = xf(i,j).*cos(alpha_max*sin(2*pi*f*t_plus + phi));
        v2(2,1) = yminus;
        v2(3,1) = h_c*c_r*sin(2*pi*f*t_plus) + xf(i,j).*sin(alpha_max*sin(2*pi*f*t_plus + phi)); 
        
        v3(1,1) = xplus.*cos(alpha_max*sin(2*pi*f*t_plus + phi));
        v3(2,1) = yf(j);
        v3(3,1) = h_c*c_r*sin(2*pi*f*t_plus) + xplus.*sin(alpha_max*sin(2*pi*f*t_plus + phi)); 
        
        v4(1,1) = xminus.*cos(alpha_max*sin(2*pi*f*t_plus + phi));
        v4(2,1) = yf(j);
        v4(3,1) = h_c*c_r*sin(2*pi*f*t_plus) + xminus.*sin(alpha_max*sin(2*pi*f*t_plus + phi)); 
        
        % Vector in y-direction is defined as v1-v2 while vector in
        % x-direction is defined as v3-v4.
        vy = v1-v2;
        vx = v3-v4;
        vn = cross(vx,vy);
        
        vn_np(i,j,1) = vn(1)/norm(vn);
        vn_np(i,j,2) = vn(2)/norm(vn);
        vn_np(i,j,3) = vn(3)/norm(vn);
        
        
        
        % Calculating the position of a point at time t_minus.
        xM_m = xf(i,j).*cos(alpha_max*sin(2*pi*f*t_minus + phi));
        yM_m = yf(j);
        zM_m = h_c*c_r*sin(2*pi*f*t_minus) + xf(i,j).*sin(alpha_max*sin(2*pi*f*t_minus + phi)); 
        
        P_nm(i,j,1) = xM_m;
        P_nm(i,j,2) = yM_m;
        P_nm(i,j,3) = zM_m;
        
        % Calculating positions of four points around position of neutral
        % plane for a time t_plus.
        v1(1,1) = xf(i,j).*cos(alpha_max*sin(2*pi*f*t_minus + phi));
        v1(2,1) = yplus;
        v1(3,1) = h_c*c_r*sin(2*pi*f*t_minus) + xf(i,j).*sin(alpha_max*sin(2*pi*f*t_minus + phi)); 
        
        v2(1,1) = xf(i,j).*cos(alpha_max*sin(2*pi*f*t_minus + phi));
        v2(2,1) = yminus;
        v2(3,1) = h_c*c_r*sin(2*pi*f*t_minus) + xf(i,j).*sin(alpha_max*sin(2*pi*f*t_minus + phi)); 
        
        v3(1,1) = xplus.*cos(alpha_max*sin(2*pi*f*t_minus + phi));
        v3(2,1) = yf(j);
        v3(3,1) = h_c*c_r*sin(2*pi*f*t_minus) + xplus.*sin(alpha_max*sin(2*pi*f*t_minus + phi)); 
        
        v4(1,1) = xminus.*cos(alpha_max*sin(2*pi*f*t_minus + phi));
        v4(2,1) = yf(j);
        v4(3,1) = h_c*c_r*sin(2*pi*f*t_minus) + xminus.*sin(alpha_max*sin(2*pi*f*t_minus + phi)); 
        
        % Vector in y-direction is defined as v1-v2 while vector in
        % x-direction is defined as v3-v4.
        vy = v1-v2;
        vx = v3-v4;
        vn = cross(vx,vy);
        
        vn_nm(i,j,1) = vn(1)/norm(vn);
        vn_nm(i,j,2) = vn(2)/norm(vn);
        vn_nm(i,j,3) = vn(3)/norm(vn);
    end
end

% Stepping through each spanwise position.
for j = 1:Ny
    
    % Creating a three dimensional array of the coordinates of the upper
    % and lower surfaces.  The row dimension goes from the leading-edge to
    % the trailing-edge in the chord direction.  The column dimension goes
    % from the root to the tip is the spanwise direction. The slice (1, 2, 
    % or 3) corresponds to the x, y or z coordinate.  Each array (u-upper,
    % l-lower) is Ny by Nx by 3 in size.
    P_u(:,j,1) = P_n(:,j,1) + diag(thick(:,j))*vn_n(:,j,1);
    P_u(:,j,2) = P_n(:,j,2) + diag(thick(:,j))*vn_n(:,j,2);
    P_u(:,j,3) = P_n(:,j,3) + diag(thick(:,j))*vn_n(:,j,3);
    
    P_l(:,j,1) = P_n(end:-1:1,j,1) - diag(thick(end:-1:1,j))*vn_n(end:-1:1,j,1);
    P_l(:,j,2) = P_n(end:-1:1,j,2) - diag(thick(end:-1:1,j))*vn_n(end:-1:1,j,2);
    P_l(:,j,3) = P_n(end:-1:1,j,3) - diag(thick(end:-1:1,j))*vn_n(end:-1:1,j,3);
    
    P_up(:,j,1) = P_np(:,j,1) + diag(thick(:,j))*vn_np(:,j,1);
    P_up(:,j,2) = P_np(:,j,2) + diag(thick(:,j))*vn_np(:,j,2);
    P_up(:,j,3) = P_np(:,j,3) + diag(thick(:,j))*vn_np(:,j,3);
    
    P_lp(:,j,1) = P_np(end:-1:1,j,1) - diag(thick(end:-1:1,j))*vn_np(end:-1:1,j,1);
    P_lp(:,j,2) = P_np(end:-1:1,j,2) - diag(thick(end:-1:1,j))*vn_np(end:-1:1,j,2);
    P_lp(:,j,3) = P_np(end:-1:1,j,3) - diag(thick(end:-1:1,j))*vn_np(end:-1:1,j,3);
    
    P_um(:,j,1) = P_nm(:,j,1) + diag(thick(:,j))*vn_nm(:,j,1);
    P_um(:,j,2) = P_nm(:,j,2) + diag(thick(:,j))*vn_nm(:,j,2);
    P_um(:,j,3) = P_nm(:,j,3) + diag(thick(:,j))*vn_nm(:,j,3);
    
    P_lm(:,j,1) = P_nm(end:-1:1,j,1) - diag(thick(end:-1:1,j))*vn_nm(end:-1:1,j,1);
    P_lm(:,j,2) = P_nm(end:-1:1,j,2) - diag(thick(end:-1:1,j))*vn_nm(end:-1:1,j,2);
    P_lm(:,j,3) = P_nm(end:-1:1,j,3) - diag(thick(end:-1:1,j))*vn_nm(end:-1:1,j,3);
    
    % If the tip goes to a point the measured chord will go to zero.  When
    % this happens the y_c coordinates become NaNs, generating NaNs for all
    % of the y and z coordinates of the upper and lower surfaces.  Instead
    % this condition replaces those NaNs with the appropriate values.  From
    % a computational point of view, all of these points coincide and could
    % be represented by a single separate value instead of Ny copies of the
    % same coordinate.  However, it is simpler to keep all of the
    % coordinates in a single array.
    if c(j) == 0
        P_u(:,j,1) = LE(j)*ones(Nx,1);
        P_u(:,j,3) = zeros(Nx,1);
        
        P_l(:,j,1) = LE(j)*ones(Nx,1);
        P_l(:,j,3) = zeros(Nx,1);
        
        P_up(:,j,1) = LE(j)*ones(Nx,1);
        P_up(:,j,3) = zeros(Nx,1);
        
        P_lp(:,j,1) = LE(j)*ones(Nx,1);
        P_lp(:,j,3) = zeros(Nx,1);
        
        P_um(:,j,1) = LE(j)*ones(Nx,1);
        P_um(:,j,3) = zeros(Nx,1);
        
        P_lm(:,j,1) = LE(j)*ones(Nx,1);
        P_lm(:,j,3) = zeros(Nx,1);
    end
end

P = [P_l;P_u(2:end,:,:)];
Pp = [P_lp;P_up(2:end,:,:)];
Pm = [P_lm;P_um(2:end,:,:)];
thick = [thick(end:-1:1,:);thick(2:end,:)];
% Calculating panel corner points, collocation points and local panel axes.
k = 1;

for j = 1:Ny - 1
    for i = 1:2*(Nx - 1)
        pc(:,k) = [P(i,j,1);P(i,j,2);P(i,j,3);  P(i,j+1,1);P(i,j+1,2);P(i,j+1,3);  P(i+1,j+1,1);P(i+1,j+1,2);P(i+1,j+1,3);  P(i+1,j,1);P(i+1,j,2);P(i+1,j,3)];
        pcp(:,k) = [Pp(i,j,1);Pp(i,j,2);Pp(i,j,3);  Pp(i,j+1,1);Pp(i,j+1,2);Pp(i,j+1,3);  Pp(i+1,j+1,1);Pp(i+1,j+1,2);Pp(i+1,j+1,3);  Pp(i+1,j,1);Pp(i+1,j,2);Pp(i+1,j,3)];
        pcm(:,k) = [Pm(i,j,1);Pm(i,j,2);Pm(i,j,3);  Pm(i,j+1,1);Pm(i,j+1,2);Pm(i,j+1,3);  Pm(i+1,j+1,1);Pm(i+1,j+1,2);Pm(i+1,j+1,3);  Pm(i+1,j,1);Pm(i+1,j,2);Pm(i+1,j,3)];

        % Determining the vertical displacement from the neutral surface of
        % the kth panel.  This is used to move the collocation points
        % inside the body without having the them cross each other.
        tp(k) = 1/4*(thick(i,j) + thick(i,j+1) + thick(i+1,j+1) + thick(i+1,j));
        
        k = k + 1;
    end
end

% Calculating panel properties.
[cptsp,~,~,~,~,~,~,Sp] = Panel(pcp);
[cptsm,~,~,~,~,~,~,Sm] = Panel(pcm);
[cpts,vc,vn,vt,~,~,~,S] = Panel(pc);

% Calculating velocity of the collocation point
Vc = (cptsp - cptsm)/Tstep/2;

% Calculating the panel length in two directions for on-body velocity
% calculation.
cf = kron([1 1 0 0],eye(3))*pc/2;
cb = kron([0 0 1 1],eye(3))*pc/2;
tl = kron([1 0 0 1],eye(3))*pc/2;
tr = kron([0 1 1 0],eye(3))*pc/2;

dl = (sqrt(sum((cpts - cf).^2,1)) + sqrt(sum((cb - cpts).^2,1)))/2';
dm_vec = ((cpts - tl) + (tr - cpts))/2;
dm = (sqrt(sum((cpts - tl).^2,1)) + sqrt(sum((tr - cpts).^2,1)))/2';
dm_hat = dm_vec*diag(dm)^-1;

% Moving the collocation points inside the body by a distance thick/2 along the
% normal direction.
% for k = 1:2*Npanels
%     del(:,k) = 0.005*tp(k)*vn(:,k);
%     cpts(:,k) = cpts(:,k) - del(:,k);
% end