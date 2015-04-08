% Written by Keith Moored, 12/15/10
%
%
% This program is a 2D panel code using source elements, following
% example 1 in section 11.2.2 of Katz and Plotkin.

clear
clc

% function [xp,Cp,Cd,xc,Qt,dL] = PanelCode2D_Cp(alpha,Npanels,Re)

% Parameters: number of panels, chord length, freestream velocity 
% magnitude, angle of attack, density of fluid.

% Fluid parameters
rho = 998;                                     % kg/m^3
nu = 1.004e-6;                                  % m^2/s
Re = 2*10^6;

Npanels = 200;
c = 8*(2.54/100);                      % m
b = 2/3*c;
alpha = 0*(pi/180);         % Enter in degrees, converts to radians.
tmax = 0.12;
Qinf = Re*nu/c;             % m/s


% Specifying the surface shape function (Van de Vooren airfoil) and panel 
% corner points.
% [xpt,zpt,xpb,zpb,a_vdv,k,epsilon,theta1,theta2,theta,r1,r2] = VanDeVooren(c,Npanels+1,tmax/2);
[xpt,zpt,xpb,zpb,~,~,~,~,~,~,~,r2] = NACA(c,Npanels+1,tmax);
xp = [xpb;xpt(2:end)];
% zp = [zpb;zpt(2:end)]
zp = [zpb(1:end-1)-zpb(1);0;zpt(2:end)-zpt(end)];
xw = [xp(end);101*c];
zw = [0;0];

% Locations of collocation points (mid-panel on the surface).
xc = (xp(2:end) + xp(1:end-1))/2;
zc = (zp(2:end) + zp(1:end-1))/2;
xwc = (xw(2)+xw(1))/2;
zwc = 0;

% Normal vectors for the panels at collocation points.  Each row is a new
% collocation point.
t = [(xp(2:end) - xp(1:end-1)) (zp(2:end) - zp(1:end-1))];
t = diag(1./sqrt(t(:,1).^2 + t(:,2).^2))*t;
n = [-t(:,2) t(:,1)];
tw = [1 0];
nw = [0 1];

% Moving collocation points inward 5% of their zc location.
% dLapprox = 2*c/Npanels;
% num1 = dLapprox*ones(Npanels,1) <= abs(zc);
% num2 = not(num1);
% xc = xc - abs(0.05*dLapprox*ones(Npanels,1)).*num1.*n(:,1) - abs(0.25*zc).*num2.*n(:,1);
% zc = zc - abs(0.05*dLapprox*ones(Npanels,1)).*num1.*n(:,2) - abs(0.25*zc).*num2.*n(:,2);
val = 0.25;
xc = xc - abs(val*zc).*n(:,1);
zc = zc - abs(val*zc).*n(:,2);

% Calculating induced velocities and influence coefficients (a = q*n, q = 
% [u, w]') at each collocation point.  The row denotes the collocation
% point while the column denotes the vortex element.  Also calculating RHS,
% RHS = -Qinf*n.
Uinf = Qinf*cos(alpha);
Winf = Qinf*sin(alpha);


%     % Influence of the body panels
%     [dPhi_s,dPhi_d,dLtemp] = Phi(ones(Npanels,1),ones(Npanels,1),Xc,Zc,xp(:,i_t),zp(:,i_t),vt(:,:,i_t),vn(:,:,i_t));
%     dL(:,i_t) = dLtemp(:,1);
%     B = dPhi_s;
%     Cb = dPhi_d;
%     
%      
%     % Influence of the trailing edge panel
%     Cte = zeros(Npanels,Npanels);
%    
%     [~,dPhi_d,~] = Phi(1,0,Xc,Zc,xTE(:,i_t),zTE(:,i_t),vtTE(1,:,i_t),vnTE(1,:,i_t));
%     Cte(:,Npanels) = dPhi_d;
%     Cte(:,1) = -dPhi_d;
%     
%     
%     % Influence of the wake panels
%     if t > 0
%         xwtemp = [xw(2,2);xw(1,2:end)'];
%         zwtemp = [zw(2,2);zw(1,2:end)'];
%         vtwtemp = zeros(Ncyc*Nstep,2);
%         vnwtemp = zeros(Ncyc*Nstep,2);
%         for i = 2:Ncyc*Nstep+1
%             vtwtemp(i-1,:) = vtw(1,:,i);
%             vnwtemp(i-1,:) = vnw(1,:,i);
%         end
%         [~,dPhi_d,~] = Phi(ones(length(muW)-1,1),0*muW',Xc,Zc,xwtemp,zwtemp,vtwtemp,vnwtemp);
%         Cw(:,2:end) = dPhi_d;
%     end
% 
%     
%     % Calculating source strengths.  
%     Vt = kron(ones(Npanels,1),[Uinf Winf]);
%     sigma(:,i_t) = sum(vn.*Vt,2);



for i = 1:Npanels
    for j = 1:Npanels
%         if i == j
% %             [dPhi_s,dL_temp] = PhiCS(1,xc(i),zc(i),xp(j),zp(j),xp(j+1),zp(j+1),t(j,:),n(j,:));
% %             B(i,j) = dPhi_s;
% %             C(i,j) = 1/2;
%             dL(j) = dL_temp;
%         else

            [dPhi_s,dL_temp] = PhiCS(1,xc(i),zc(i),xp(j),zp(j),xp(j+1),zp(j+1),t(j,:),n(j,:));
            [dPhi_d] = PhiCD(1,xc(i),zc(i),xp(j),zp(j),xp(j+1),zp(j+1),t(j,:),n(j,:));
            B(i,j) = dPhi_s;
            C(i,j) = dPhi_d;
            if i == j
                dL(j,1) = dL_temp;
            end
                
%         end
    end
    sigma(i,1) = -n(i,:)*[Uinf Winf]';
    [dPhi_d] = PhiCD(1,xc(i),zc(i),xw(1),zw(1),xw(2),zw(2),tw,nw);
    C(i,Npanels+1) = dPhi_d;
end
sigma;

% Applying the Kutta Condition
A = [C(:,1) - C(:,Npanels + 1), C(:,2:Npanels - 1), C(:,Npanels) + C(:,Npanels + 1)];
RHS = -B*sigma;

% Solving for doublet strengths.
mu = A\RHS;

% Calculating on-body velocities and pressures. Pressures are calculated
% from Bernoulli's equation.  The normal velocities on the body should be
% zero.  
for i = 2:Npanels - 1
    Qt(i,1) = (mu(i-1) - mu(i + 1))/(dL(i+1)/2 + dL(i) + dL(i-1)/2) + Uinf*t(i,1)+ Winf*t(i,2);
%     Qt(i,1) = 2*(mu(i) - mu(i + 1))/(dL(i+1) + dL(i)) + Uinf*(t(i,1) + t(i+1,1))/2 + Winf*(t(i,2) + t(i+1,2))/2;
end


% First and last panel velocity calculation
Qt(1,1) = (mu(1) - mu(2))/(dL(2)/2 + 3/2*dL(1)) + Uinf*t(1,1)/2 + Winf*t(1,2)/2;
Qt(Npanels,1) = (mu(Npanels - 1) - mu(Npanels))/(dL(Npanels - 1)/2 + 3/2*dL(Npanels)) + Uinf*t(Npanels,1)/2 + Winf*t(Npanels,2)/2;

% Calculating the pressure coefficient.
Cp = 1 - Qt.^2/Qinf^2;

% Calculating viscous drag correction.               
[dFi_pct,dL_pct,factor,xbc_sep,xtc_sep] = SkinFrictionSolver(xc,Qt,nu,rho,c,tmax,Cp,Qinf);

dFi = dFi_pct.*kron(dL*b,ones(factor,1));
dFshear = kron(eye(length(xc)),ones(1,factor))*dFi;
dtau = dFshear./(dL*b);
   
% Calculating lift on the body and pitching moment about the leading edge.
for j = 1:Npanels
%     alphaj(j) = acos((n(j,2) + n(j+1,2))/2);
%     delCl(j,1) = -Cp(j)*(dL(j+1) + dL(j))/2*cos(alpha + alphaj(j))/c;
    delF(:,j) = -Cp(j)*1/2*rho*Qinf^2*dL(j)*b*n(j,:)' + dFshear(j)*t(j,:)';
%     delF(:,j) = -Cp(j)*1/2*rho*Qinf^2*dL(j)*n(j,:)';
    %     delL(j,1) = delCl(j)*1/2*rho*Qinf^2*c;
end

% D_visc = sum(dDf);

F = sum(delF,2);
CF = norm(F)/(1/2*rho*Qinf^2*c);

L = F(2)*cos(alpha) - F(1)*sin(alpha);
D = F(2)*sin(alpha) + F(1)*cos(alpha)


% M0 = -sum(delL*cos(alpha).*xp(2:end-1));

% Calculating non-dimensional coefficients.
% Cl = sum(delCl)
% Cm = M0/(1/2*rho*Qinf^2*c^2)
Gamma = mu(1) - mu(Npanels);
Gamma_ideal = pi*c*Qinf*alpha;
Cl = L/(1/2*rho*Qinf^2*c)
Cl_ideal = 2*pi*sin(alpha);
Cd = D/(1/2*rho*Qinf^2*c)
Cl_KJ = rho*Qinf*Gamma/(1/2*rho*Qinf^2*c);

thickness = (max(zp) - min(zp))/c;



% Analytical solution.

[xpt,zpt,xpb,zpb,a_vdv,k,epsilon,theta1,theta2,theta,r1,r2] = VanDeVooren(c,200,tmax/2);

xp_a = [xpb;xpt(2:end)];
% zp = [zpb;zpt(2:end)];
zp_a = [zpb(1:end-1)-zpb(1);0;zpt(2:end)-zpt(end)];

Aa = cos((k-1)*theta1).*cos(k*theta2) + sin((k-1)*theta1).*sin(k*theta2);
Ba = sin((k-1)*theta1).*cos(k*theta2) - cos((k-1)*theta1).*sin(k*theta2);
D0 = a_vdv*(1 - k + k*epsilon);
D1 = Aa.*(a_vdv*cos(theta) - D0) - Ba.*a_vdv.*sin(theta);
D2 = Aa.*a_vdv.*sin(theta) + Ba.*(a_vdv.*cos(theta) - D0);

u_anltc = 2*Qinf*r2.^k./r1.^(k-1).*(sin(alpha) - sin(alpha - theta))./(D1.^2 + D2.^2).*(D1.*sin(theta) + D2.*cos(theta));
w_anltc = -2*Qinf*r2.^k./r1.^(k-1).*(sin(alpha) - sin(alpha - theta))./(D1.^2 + D2.^2).*(D1.*cos(theta) - D2.*sin(theta));

Cp_anltc_t = 1 - (u_anltc.^2 + w_anltc.^2)/Qinf^2;

alpha_b = -alpha;
[xpt,zpt,xpb,zpb,a_vdv,k,epsilon,theta1,theta2,theta,r1,r2] = VanDeVooren(c,200,tmax/2);

xp_a = [xpb;xpt(2:end)];
% zp = [zpb;zpt(2:end)];
zp_a = [zpb(1:end-1)-zpb(1);0;zpt(2:end)-zpt(end)];

Aa = cos((k-1)*theta1).*cos(k*theta2) + sin((k-1)*theta1).*sin(k*theta2);
Ba = sin((k-1)*theta1).*cos(k*theta2) - cos((k-1)*theta1).*sin(k*theta2);
D0 = a_vdv*(1 - k + k*epsilon);
D1 = Aa.*(a_vdv*cos(theta) - D0) - Ba.*a_vdv.*sin(theta);
D2 = Aa.*a_vdv.*sin(theta) + Ba.*(a_vdv.*cos(theta) - D0);

u_anltc = 2*Qinf*r2.^k./r1.^(k-1).*(sin(alpha_b) - sin(alpha_b - theta))./(D1.^2 + D2.^2).*(D1.*sin(theta) + D2.*cos(theta));
w_anltc = -2*Qinf*r2.^k./r1.^(k-1).*(sin(alpha_b) - sin(alpha_b - theta))./(D1.^2 + D2.^2).*(D1.*cos(theta) - D2.*sin(theta));

Cp_anltc_b = 1 - (u_anltc.^2 + w_anltc.^2)/Qinf^2;


% % Plotting discretized airfoil with normals and lift vectors.
% figure
% hold on
% plot(xp,zp,'+-k')
% plot(xc,zc,'xr')
% quiver(xp(1:end-1),zp(1:end-1),n(:,1),n(:,2),'k')
% % quiver(xp(1:end-1),zp(1:end-1),t(:,1),t(:,2),'k')
% axis([-c/2 3*c/2 -c c])
% axis equal
% xlabel('$$x$$','interpreter','latex','fontsize',12,'fontname','TimesNewRoman')
% ylabel('$$z$$','interpreter','latex','fontsize',12,'fontname','TimesNewRoman')

% Plotting pressure coefficient both numerical and analytical solutions.
figure
hold on
plot(xp_a(round(end/2):end),-Cp_anltc_t,'-k')
plot(xp_a(round(end/2):end),-Cp_anltc_b,'-k')
plot(xc(1:Npanels/2),-Cp(1:Npanels/2),'sk')
plot(xc(Npanels/2 + 1:end),-Cp(Npanels/2 + 1:end),'^k')

% plot(xp(2:round((Npanels + 1)/2)),-Cp(1:round((Npanels+1)/2)-1),'sk')
% plot(xp(round((Npanels + 1)/2):end-1),-Cp(round((Npanels+1)/2)-1:end),'^k')
xlabel('$$\frac{x}{c}$$','interpreter','latex','fontname','TimesNewRoman','fontsize',12)
ylabel('$$C_p$$','interpreter','latex','fontname','TimesNewRoman','fontsize',12)
% axis([0 c -1 2])
% axis([-1.6 -1.3 -1 2])
legend('Analytical','Computational','Location','South')
% 



% % Calculating the flowfield around the wing.
% nx = 40;
% ny = 40;
% 
% U = Uinf*ones(ny,nx);
% W = Winf*ones(ny,nx);
% 
% xf = linspace(-c,3*c,nx)';
% yf = linspace(-2*c,2*c,ny)';
% 
% [Xf,Yf] = meshgrid(xf,yf);
% 
% for i = 1:ny
%     for j = 1:nx
%         for k = 1:Npanels
%             [u_newS,w_newS] = Sor2DC(sigma(k),Xf(i,j),Yf(i,j),xp(k),zp(k),xp(k+1),zp(k+1),t(k,:),n(k,:));
%             [u_newD,w_newD] = Dub2DC(mu(k),Xf(i,j),Yf(i,j),xp(k),zp(k),xp(k+1),zp(k+1),t(k,:),n(k,:));
%             u_temp(k,1) = u_newD;
%             u_temp(k,2) = u_newS;
%             w_temp(k,1) = w_newD;
%             w_temp(k,2) = w_newS;
%         end
%         [u_new,w_new] = Dub2DC(mu(Npanels)-mu(1),Xf(i,j),Yf(i,j),xw(1),zw(1),xw(2),zw(2),tw,nw);
%         u_temp(Npanels+1,1) = u_new;
%         w_temp(Npanels+1,1) = w_new;
%         u_temp(Npanels+1,2) = 0;
%         w_temp(Npanels+1,2) = 0;
%         u_d(i,j) = sum(u_temp(:,1));
%         w_d(i,j) = sum(w_temp(:,1));
%         u_s(i,j) = sum(u_temp(:,2));
%         w_s(i,j) = sum(w_temp(:,2));
%         
%         U(i,j) = U(i,j) + u_d(i,j) + u_s(i,j);
%         W(i,j) = W(i,j) + w_d(i,j) + w_s(i,j);
%         u_temp = [];
%         w_temp = [];
%     end
% end
% 
% % Plotting discretized airfoil with normals and lift vectors.
% figure
% hold on
% plot(xp,zp,'+-k')
% plot(xc,zc,'xb')
% streamline(Xf,Yf,U,W,Xf(:,1),Yf(:,1))
% % quiver(Xf,Yf,U,W)
% axis equal
% axis([-c 3*c -c c])
% title('2D Solution','interpreter','latex','fontsize',14,'fontname','TimesNewRoman')
% xlabel('$$x$$','interpreter','latex','fontsize',12,'fontname','TimesNewRoman')
% ylabel('$$z$$','interpreter','latex','fontsize',12,'fontname','TimesNewRoman')
% 
% % Plotting discretized airfoil with normals and lift vectors.
% figure
% hold on
% plot(xp,zp,'+-k')
% plot(xc,zc,'xb')
% quiver(Xf,Yf,U,W,'b')
% axis equal
% axis([-c 4*c -c c])
% title('2D Solution','interpreter','latex','fontsize',14,'fontname','TimesNewRoman')
% xlabel('$$x$$','interpreter','latex','fontsize',12,'fontname','TimesNewRoman')
% ylabel('$$z$$','interpreter','latex','fontsize',12,'fontname','TimesNewRoman')
% 
% % Plotting discretized airfoil with normals and lift vectors.
% figure
% hold on
% plot(xp,zp,'+-k')
% plot(xc,zc,'xb')
% % quiver(Xf,Yf,u_d,w_d,'b')
% axis equal
% axis([-c 4*c -c c])
% title('2D Solution','interpreter','latex','fontsize',14,'fontname','TimesNewRoman')
% xlabel('$$x$$','interpreter','latex','fontsize',12,'fontname','TimesNewRoman')
% ylabel('$$z$$','interpreter','latex','fontsize',12,'fontname','TimesNewRoman')
% 
% % Plotting discretized airfoil with normals and lift vectors.
% figure
% hold on
% plot(xp,zp,'+-k')
% plot(xc,zc,'xb')
% quiver(Xf,Yf,u_s,w_s,'b')
% axis equal
% axis([-c 4*c -c c])
% title('2D Solution','interpreter','latex','fontsize',14,'fontname','TimesNewRoman')
% xlabel('$$x$$','interpreter','latex','fontsize',12,'fontname','TimesNewRoman')
% ylabel('$$z$$','interpreter','latex','fontsize',12,'fontname','TimesNewRoman')

