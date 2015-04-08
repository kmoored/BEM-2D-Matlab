% clear
% clc

function [dFi,dL_pct,factor,xbc_sep,xtc_sep] = SkinFrictionSolver(xc,Qt,nu,rho,c,tmax,Cp,Qinf)

% Re = 6e6;
% trans = 5e5;
% alpha = 15*(pi/180);
% N = 150;
% nu = 1.004e-6;              % m^2/s
% 
% [xp,~,~,~,xc,Qt,dL] = PanelCode2D_Cp(alpha,N,Re,trans);

%% Interpolating points for a more refined analysis

% factor = 5;
% Npan = length(xc);

%% Interpolating points for a more refined analysis

Npan = length(xc);

if Npan >= 200
    factor = 1;
else
    factor = round(200/Npan);
end

%% Calculating xp, zp

% NACA four series shape coefficients.
a0 = 0.2969;
a1 = -0.1260;
a2 = -0.3516;
a3 = 0.2843;
a4 = -0.1015;

% NACA four series shape equation.
xb = linspace(pi,0,(Npan/2) + 1)';
x_c = 1/2*(1 - cos(xb));
xb = x_c*c;
zb = -tmax*c/0.2*(a0*sqrt(x_c) + a1*(x_c) + a2*(x_c).^2 + a3*(x_c).^3 + a4*(x_c).^4);

% NACA four series shape equation.
xt = linspace(0,pi,(Npan/2) + 1)';
x_c = 1/2*(1 - cos(xt));
xt = x_c*c;
zt = tmax*c/0.2*(a0*sqrt(x_c) + a1*(x_c) + a2*(x_c).^2 + a3*(x_c).^3 + a4*(x_c).^4);

xp = [xb;xt(2:end)];
zp = [zb;zt(2:end)];

%% Calculating xpi, zpi

% NACA four series shape coefficients.
a0 = 0.2969;
a1 = -0.1260;
a2 = -0.3516;
a3 = 0.2843;
a4 = -0.1015;

% NACA four series shape equation.
xb = linspace(pi,0,factor*(Npan/2) + 1)';
x_c = 1/2*(1 - cos(xb));
xb = x_c*c;
zb = -tmax*c/0.2*(a0*sqrt(x_c) + a1*(x_c) + a2*(x_c).^2 + a3*(x_c).^3 + a4*(x_c).^4);

% NACA four series shape equation.
xt = linspace(0,pi,factor*(Npan/2) + 1)';
x_c = 1/2*(1 - cos(xt));
xt = x_c*c;
zt = tmax*c/0.2*(a0*sqrt(x_c) + a1*(x_c) + a2*(x_c).^2 + a3*(x_c).^3 + a4*(x_c).^4);

xpi = [xb;xt(2:end)];
zpi = [zb;zt(2:end)];

%% Calculating interpolated velocity field
dLi = sqrt((xpi(2:end) - xpi(1:end-1)).^2 + (zpi(2:end) - zpi(1:end-1)).^2);

xci = (xpi(2:end) + xpi(1:end-1))/2;
Qtp = (Qt(2:end) + Qt(1:end-1))/2;
Qtp = [0;Qtp;0];

XC = xc;
QT = Qt;

Xp = xp.*sign(zp);
Xi = xpi.*sign(zpi);

Qti = interp1([-c;Xp(2:end-1);c],Qtp,[-c;Xi(2:end-1);c],'cubic');

Qt = (Qti(2:end) + Qti(1:end-1))/2;
xc = xci;
dL = dLi;

% % Plotting_________________________________________
% figure
% set(gcf,'DefaultAxesfontsize',16,'DefaultAxesfontname','TimesNewRoman','DefaultAxesGridLineStyle','-.')
% xlabel('$$x/c$$','interpreter','latex','fontname','TimesNewRoman','fontsize',30)
% ylabel('$$Q$$','interpreter','latex','fontname','TimesNewRoman','fontsize',30,'rotation',0)
% 
% hold on
% 
% plot(XC,QT,'ok','linewidth',2,'markersize',12)
% plot(xc,Qt,'-k','linewidth',2,'markersize',12)
% 
% legend('Top Surface','Bottom Surface','Location','NorthEast')
% % axis([-0.1 1.1 0 max(Qt)*1.2])


%% Finding the stagnation point location.
stagpts = find((Qt(1:end-1).*Qt(2:end)) <= 0);
if isempty(stagpts)
    [~,stagpts] = min((Qt(1:end-1).*Qt(2:end)));
end
LE = Npan*factor/2;
[val,ind] = min(abs(stagpts - LE*ones(length(stagpts),1)));
stagpt = stagpts(ind);

% Splitting the surface into the top surface above the stagantion point and
% the bottom surface below the stagnation point.
Qt_b = -Qt(stagpt:-1:1);
Qt_t = Qt(stagpt + 1:end);


% Calculating the length steps along the surface for both the top and
% bottom surface.
dels_b = dL(stagpt:-1:1);
dels_t = dL(stagpt + 1:end);

% Calculating the top surface coordinate and the bottom surface coordinate.
for i = 1:length(dels_b)
    if i == 1
        s_b(i,1) = dels_b(i)/2;
    else 
        s_b(i,1) = dels_b(i)/2 + sum(dels_b(1:i-1));
    end
end
for i = 1:length(dels_t)
    if i == 1
        s_t(i,1) = dels_t(i)/2;
    else 
        s_t(i,1) = dels_t(i)/2 + sum(dels_t(1:i-1));
    end
end

% Calculating the Reynolds number along the top and bottom surface
Re_xb = Qt_b.*s_b/nu;
Re_xt = Qt_t.*s_t/nu;


%% Laminar Boundary Layer Calculation (Method of Thwaites)

xc_b = xc(stagpt:-1:1);
xc_t = xc(stagpt + 1:end);

% % Plotting_________________________________________
% figure
% set(gcf,'DefaultAxesfontsize',16,'DefaultAxesfontname','TimesNewRoman','DefaultAxesGridLineStyle','-.')
% xlabel('$$x/c$$','interpreter','latex','fontname','TimesNewRoman','fontsize',30)
% ylabel('$$Q$$','interpreter','latex','fontname','TimesNewRoman','fontsize',30,'rotation',0)
% 
% hold on
% 
% plot(s_t,Qt_t,'-k','linewidth',2,'markersize',12)
% plot(s_b,Qt_b,':k','linewidth',2,'markersize',12)
% 
% legend('Top Surface','Bottom Surface','Location','NorthEast')
% axis([-0.1 1.1 0 max(Qt)*1.2])

% Calculations______________________________________

Ibot_pri = Qt_b.^5.*dels_b;
Itop_pri = Qt_t.^5.*dels_t;


for i = 1:length(s_b)
    Ibot(i,1) = sum(Ibot_pri(1:i));
end
for i = 1:length(s_t)
    Itop(i,1) = sum(Itop_pri(1:i));
end

theta_b = sqrt(0.45*nu*Ibot./Qt_b.^6);
theta_t = sqrt(0.45*nu*Itop./Qt_t.^6);

Re_thetab = Qt_b.*theta_b/nu;
Re_thetat = Qt_t.*theta_t/nu;


% % Plotting_________________________________________
% figure
% set(gcf,'DefaultAxesfontsize',16,'DefaultAxesfontname','TimesNewRoman','DefaultAxesGridLineStyle','-.')
% xlabel('$$s$$','interpreter','latex','fontname','TimesNewRoman','fontsize',30)
% ylabel('$$Re_x$$','interpreter','latex','fontname','TimesNewRoman','fontsize',30,'rotation',0)
% 
% hold on
% 
% plot(s_t(1:end-5),2.9*Re_xt(1:end-5).^(0.4),'-k','linewidth',2,'markersize',12)
% plot(s_b(1:end-5),2.9*Re_xb(1:end-5).^(0.4),':k','linewidth',2,'markersize',12)
% plot(s_t(1:end-5),Re_thetat(1:end-5),'-b','linewidth',2,'markersize',12)
% plot(s_b(1:end-5),Re_thetab(1:end-5),':b','linewidth',2,'markersize',12)
% 
% legend('Top Surface','Bottom Surface','Location','NorthWest')
% 

% % Plotting_________________________________________
% figure
% set(gcf,'DefaultAxesfontsize',16,'DefaultAxesfontname','TimesNewRoman','DefaultAxesGridLineStyle','-.')
% xlabel('$$x/c$$','interpreter','latex','fontname','TimesNewRoman','fontsize',30)
% ylabel('$$\theta$$','interpreter','latex','fontname','TimesNewRoman','fontsize',30,'rotation',0)
% 
% hold on
% 
% plot(s_t(1:end-1),theta_t(1:end-1),'-k','linewidth',2,'markersize',12)
% plot(s_b(1:end-1),theta_b(1:end-1),':k','linewidth',2,'markersize',12)
% 
% legend('Top Surface','Bottom Surface','Location','NorthWest')
% axis([-0.1 1.1 0 0.04])

% Calculating the velocity gradients______________________________________
dUedx_b = zeros(length(s_b),1);
dUedx_t = zeros(length(s_t),1);

if numel(Qt_b) > 2
    dUedx_b(1,1) = (Qt_b(2) - Qt_b(1))./(dels_b(2)/2 + dels_b(1)/2); 
    dUedx_b(2:end-1) = (Qt_b(3:end) - Qt_b(1:end-2))./(dels_b(3:end)/2 + dels_b(2:end-1) + dels_b(1:end-2)/2); 
    dUedx_b(end,1) = (Qt_b(end) - Qt_b(end - 1))./(dels_b(end)/2 + dels_b(end - 1)/2);
end

if numel(Qt_t) > 2
    dUedx_t(1,1) = (Qt_t(2) - Qt_t(1))./(dels_t(2)/2 + dels_t(1)/2); 
    dUedx_t(2:end-1) = (Qt_t(3:end) - Qt_t(1:end-2))./(dels_t(3:end)/2 + dels_t(2:end-1) + dels_t(1:end-2)/2); 
    dUedx_t(end,1) = (Qt_t(end) - Qt_t(end - 1))./(dels_t(end)/2 + dels_t(end - 1)/2);
end


lambda_t = theta_t.^2.*dUedx_t/nu;
lambda_b = theta_b.^2.*dUedx_b/nu;
lambda_t = lambda_t.*(lambda_t >= -0.09) + -0.09*(lambda_t <= -0.09);
lambda_b = lambda_b.*(lambda_b >= -0.09) + -0.09*(lambda_b <= -0.09);


% Calculating the laminar seperation point_________
SepIndex_b = find(lambda_b == -0.09,1,'first');
SepIndex_t = find(lambda_t == -0.09,1,'first');

xbc_sep = xc_b(SepIndex_b)/c;
xtc_sep = xc_t(SepIndex_t)/c;


% Calculating the shear function____________________
S_t = (lambda_t + 0.09).^0.62;
S_b = (lambda_b + 0.09).^0.62;


% % Plotting__________________________________________
% figure
% set(gcf,'DefaultAxesfontsize',16,'DefaultAxesfontname','TimesNewRoman','DefaultAxesGridLineStyle','-.')
% xlabel('$$x/c$$','interpreter','latex','fontname','TimesNewRoman','fontsize',30)
% ylabel('$$S$$','interpreter','latex','fontname','TimesNewRoman','fontsize',30,'rotation',0)
% 
% hold on
% 
% plot(s_t(1:end-1),S_t(1:end-1),'-k','linewidth',2,'markersize',12)
% plot(s_b(1:end-1),S_b(1:end-1),':k','linewidth',2,'markersize',12)
% 
% legend('Top Surface','Bottom Surface','Location','NorthEast')

% Calculating friction coefficient___________________
Cf_t = 2*nu.*S_t./(Qt_t.*theta_t);
Cf_b = 2*nu.*S_b./(Qt_b.*theta_b);


% % Plotting friction coefficient______________________
% figure
% set(gcf,'DefaultAxesfontsize',16,'DefaultAxesfontname','TimesNewRoman','DefaultAxesGridLineStyle','-.')
% xlabel('$$x/c$$','interpreter','latex','fontname','TimesNewRoman','fontsize',30)
% ylabel('$$C_f$$','interpreter','latex','fontname','TimesNewRoman','fontsize',30,'rotation',0)
% 
% hold on
% 
% plot(s_t(1:end-1),Cf_t(1:end-1),'-k','linewidth',2,'markersize',12)
% plot(s_b(1:end-1),Cf_b(1:end-1),':k','linewidth',2,'markersize',12)
% 
% legend('Top Surface','Bottom Surface','Location','NorthEast')
% axis([0 1.1 0 0.2])







%% Transition Calculations

bot_trans = find(2.9*Re_xb.^(0.4) - Re_thetab <= 0,1,'first');
top_trans = find(2.9*Re_xt.^(0.4) - Re_thetat <= 0,1,'first');

sb_trans = s_b(bot_trans);
st_trans = s_t(top_trans);

sb_sep = s_b(SepIndex_b);
st_sep = s_t(SepIndex_t);

% If the laminar separation point is upstream of the turbulent transition
% point then turbulent transition is assumed at the separation point.
if sb_sep < sb_trans
    bot_trans = SepIndex_b;
end
if st_sep < st_trans
    top_trans = SepIndex_t;
end



%% Turbulent Boundary Layer Calculations (Approximate Polynomial Solution)

eta = linspace(0,1,1000)';
u = eta.^(1/7);
H = 1.286;

if (isempty(bot_trans)) 
    delta_b_tur = [];
    delta_b_star = [];
    theta_b_tur = [];
    Re_thetab_tur = [];
    Cf_b_tur = [];
elseif (bot_trans == length(s_b)) 
    delta_b_tur = [];
    delta_b_star = [];
    theta_b_tur = [];
    Re_thetab_tur = [];
    Cf_b_tur = [];
else
    % Calculating turbulent BL parameters for the bottom surface.
    Xspan = [s_b(bot_trans) s_b(end)];
    IC = 72/7*theta_b(bot_trans);
%     Qt_b = Qt_b.*(Qt_b > 0)
%     Qt
    
    [X Delta] = ode45(@(x,delta) BL_ode(x,delta,Qt_b,dUedx_b,s_b,H,nu),Xspan,IC);

    delta_b_tur = interp1(X,Delta,s_b(bot_trans:end));
    delta_b_star = 1/8*delta_b_tur;
    theta_b_tur = 7/72*delta_b_tur;
    Re_thetab_tur = Qt_b(bot_trans:end).*theta_b_tur/nu;
    Cf_b_tur = (0.3*exp(-1.33*H))./(log10(Re_thetab_tur)).^(1.74 + 0.31*H);
end

if isempty(top_trans) 
    delta_t_tur = [];
    delta_t_star = [];
    theta_t_tur = [];
    Re_thetat_tur = [];
    Cf_t_tur = [];
elseif top_trans == length(s_t)
    delta_t_tur = [];
    delta_t_star = [];
    theta_t_tur = [];
    Re_thetat_tur = [];
    Cf_t_tur = [];
else
    % Calculating turbulent BL parameters for the top surface.
    Xspan = [s_t(top_trans) s_t(end)];
    IC = 72/7*theta_t(top_trans); 
    [X Delta] = ode45(@(x,delta) BL_ode(x,delta,Qt_t,dUedx_t,s_t,H,nu),Xspan,IC);

    delta_t_tur = interp1(X,Delta,s_t(top_trans:end));
    delta_t_star = 1/8*delta_t_tur;
    theta_t_tur = 7/72*delta_t_tur;
    Re_thetat_tur = Qt_t(top_trans:end).*theta_t_tur/nu;
    Cf_t_tur = (0.3*exp(-1.33*H))./(log10(Re_thetat_tur)).^(1.74 + 0.31*H);
end
    
% Combine laminar and turbulent boundary layer calculations.
if isempty(bot_trans)
elseif bot_trans == length(s_b)
else
    theta_b(bot_trans:end) = theta_b_tur;
    Cf_b(bot_trans:end) = Cf_b_tur;
end

if isempty(top_trans)
elseif top_trans == length(s_t)
else
    theta_t(top_trans:end) = theta_t_tur;
    Cf_t(top_trans:end) = Cf_t_tur;
end

% % Plotting momentum thickness_________________________________
% figure
% set(gcf,'DefaultAxesfontsize',16,'DefaultAxesfontname','TimesNewRoman','DefaultAxesGridLineStyle','-.')
% xlabel('$$x/c$$','interpreter','latex','fontname','TimesNewRoman','fontsize',30)
% ylabel('$$\theta$$','interpreter','latex','fontname','TimesNewRoman','fontsize',30,'rotation',0)
% 
% hold on
% 
% plot(s_t(1:end-1),theta_t(1:end-1),'-k','linewidth',2,'markersize',12)
% plot(s_b(1:end-1),theta_b(1:end-1),':k','linewidth',2,'markersize',12)
% 
% legend('Top Surface','Bottom Surface','Location','NorthWest')
% axis([-0.1 1.1 0 0.04])
% 
% % Plotting friction coefficient______________________
% figure
% set(gcf,'DefaultAxesfontsize',16,'DefaultAxesfontname','TimesNewRoman','DefaultAxesGridLineStyle','-.')
% xlabel('$$x/c$$','interpreter','latex','fontname','TimesNewRoman','fontsize',30)
% ylabel('$$C_f$$','interpreter','latex','fontname','TimesNewRoman','fontsize',30,'rotation',0)
% 
% hold on
% 
% plot(s_t(1:end-1),Cf_t(1:end-1),'-k','linewidth',2,'markersize',12)
% plot(s_b(1:end-1),Cf_b(1:end-1),':k','linewidth',2,'markersize',12)
% 
% legend('Top Surface','Bottom Surface','Location','NorthEast')
% axis([0 1.1 0 0.02])




%% Calculating turbulent separation
% Cpp = (Cp(2:end) + Cp(1:end-1))/2;
% Cpp = [0;Cpp;0];
% 
% Cpi = interp1([-c;Xp(2:end-1);c],Cpp,[-c;Xi(2:end-1);c],'cubic');
% Cp = (Cpi(2:end) + Cpi(1:end-1))/2;
% 
% % Splitting the surface into the top surface above the stagnation point and
% % the bottom surface below the stagnation point.
% Cp_b = Cp(stagpt:-1:1);
% Cp_t = Cp(stagpt + 1:end);
% 
% [~,ind_b] = min(Cp_b);
% [~,ind_t] = min(Cp_t);
% 
% % Calculating the points of minimum pressure and the velocity at those
% % points.
% xb_m = s_b(ind_b);
% xt_m = s_t(ind_t);
% 
% ub_m = Qt_b(ind_b);
% ut_m = Qt_t(ind_t);
% 
% % Calculating the points of turbulent transition and the velocity at those
% % points.
% xb_t = s_b(bot_trans);
% xt_t = s_t(top_trans);
% 
% ub_t = Qt_b(bot_trans);
% ut_t = Qt_t(top_trans);
% 
% % Calculations______________________________________
% 
% Ibot1_pri = (Qt_b/ub_m).^5.*dels_b;
% Itop1_pri = (Qt_t/ut_m).^5.*dels_t;
% 
% Ibot2_pri = (Qt_b/ub_m).^4.*dels_b;
% Itop2_pri = (Qt_t/ut_m).^4.*dels_t;
% 
% Ibot1 = sum(Ibot1_pri(1:bot_trans));
% Itop1 = sum(Itop1_pri(1:top_trans));
% 
% if xb_t < xb_m    
%     Ibot2 = -sum(Ibot2_pri(bot_trans:ind_b));
% elseif xb_t > xb_m
%     Ibot2 = sum(Ibot2_pri(ind_b:bot_trans));
% end
% 
% if xt_t < xt_m    
%     Itop2 = -sum(Itop2_pri(top_trans:ind_t));
% elseif xt_t > xt_m
%     Itop2 = sum(Itop2_pri(ind_t:top_trans));
% end
% 
% xpri_b = xb_m - 58*(nu/ub_m)*(ub_t/nu*Ibot1).^(3/5) + Ibot2; 
% xpri_t = xt_m - 58*(nu/ut_m)*(ut_t/nu*Itop1).^(3/5) + Itop2;
% 
% xbar_b = s_b - xpri_b;
% xbar_t = s_t - xpri_t;
% 
% Cpbar_b = 1 + (Cp_b - 1)*Qinf^2/ub_m^2;
% Cpbar_t = 1 + (Cp_t - 1)*Qinf^2/ut_m^2;
% 
% 
% 
% % Calculating the pressure gradients______________________________________
% dCpdx_b = zeros(length(s_b),1);
% dCpdx_t = zeros(length(s_t),1);
% 
% dCpdx_b(1,1) = (Cp_b(2) - Cp_b(1))./(dels_b(2)/2 + dels_b(1)/2); 
% dCpdx_b(2:end-1) = (Cp_b(3:end) - Cp_b(1:end-2))./(dels_b(3:end)/2 + dels_b(2:end-1) + dels_b(1:end-2)/2); 
% dCpdx_b(end,1) = (Cp_b(end) - Cp_b(end - 1))./(dels_b(end)/2 + dels_b(end - 1)/2);
% 
% dCpdx_t(1,1) = (Cp_t(2) - Cp_t(1))./(dels_t(2)/2 + dels_t(1)/2); 
% dCpdx_t(2:end-1) = (Cp_t(3:end) - Cp_t(1:end-2))./(dels_t(3:end)/2 + dels_t(2:end-1) + dels_t(1:end-2)/2); 
% dCpdx_t(end,1) = (Cp_t(end) - Cp_t(end - 1))./(dels_t(end)/2 + dels_t(end - 1)/2);
% 
% Re_xbar_b = ub_m*xbar_b/nu;
% Re_xbar_t = ut_m*xbar_t/nu;
% 
% F_b = Cpbar_b.*sqrt(xbar_b.*dCpdx_b)./(10^(-6).*Re_xbar_b).^(0.1);
% F_t = Cpbar_t.*sqrt(xbar_t.*dCpdx_t)./(10^(-6).*Re_xbar_t).^(0.1);
% 
% F_b = F_b.*(xbar_b > 0)
% F_t = F_t.*(xbar_t > 0)
% 
% max(Cpbar_b.*(xbar_b > 0))
% max(Cpbar_t.*(xbar_t > 0))
% 
% % if F_b >= 0.5


%% Gathering data
Cf = [-Cf_b(end:-1:1); Cf_t];
Cf = (imag(Cf) == 0).*Cf;
dL_act = kron(eye(Npan),ones(1,factor))*dL;
dL_pct = kron(1./dL_act,ones(factor,1)).*dL;

dFi = 1/2*rho*Qt.^2.*Cf.*dL_pct;





