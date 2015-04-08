% clear
% clc

function [dFi,dL_pct,factor,xbc_sep,xtc_sep] = SkinFrictionSolver(xc,Qt,nu,rho,c,tmax,Cp,Qinf,trip)


%% Interpolating points for a more refined analysis

Npan = length(xc);
factor = ceil(150/Npan);


%% Calculating xp, zp

% [xt,zt,xb,zb] = TearDropShape(c,Npan + 1,tmax);
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

% [xti,zti,xbi,zbi] = TearDropShape(c,factor*(Npan) + 1,tmax);
% NACA four series shape coefficients.
a0 = 0.2969;
a1 = -0.1260;
a2 = -0.3516;
a3 = 0.2843;
a4 = -0.1015;

% NACA four series shape equation.
xb = linspace(pi,0,factor*(Npan/2) + 1)';
x_c = 1/2*(1 - cos(xb));
xbi = x_c*c;
zbi = -tmax*c/0.2*(a0*sqrt(x_c) + a1*(x_c) + a2*(x_c).^2 + a3*(x_c).^3 + a4*(x_c).^4);

% NACA four series shape equation.
xt = linspace(0,pi,factor*(Npan/2) + 1)';
x_c = 1/2*(1 - cos(xt));
xti = x_c*c;
zti = tmax*c/0.2*(a0*sqrt(x_c) + a1*(x_c) + a2*(x_c).^2 + a3*(x_c).^3 + a4*(x_c).^4);

xpi = [xbi;xti(2:end)];
zpi = [zbi;zti(2:end)];

%% Calculating interpolated velocity field
dLi = sqrt((xpi(2:end) - xpi(1:end-1)).^2 + (zpi(2:end) - zpi(1:end-1)).^2);

xci = (xpi(2:end) + xpi(1:end-1))/2;
Qtp = (Qt(2:end) + Qt(1:end-1))/2;
Qtp = [0;Qtp;0];

XC = xc;
QT = Qt;
XP = xp;
ZP = zp;

Xp = xp.*sign(zp);
Xi = xpi.*sign(zpi);

Qti = interp1([-c;Xp(2:end-1);c],Qtp,[-c;Xi(2:end-1);c],'cubic');

Qt = (Qti(2:end) + Qti(1:end-1))/2;
xc = xci;
dL = dLi;

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


% Calculating friction coefficient___________________
Cf_t = 2*nu.*S_t./(Qt_t.*theta_t);
Cf_b = 2*nu.*S_b./(Qt_b.*theta_b);



%% Transition Calculations
if trip == 1
    SepIndex_b = 1;
    SepIndex_t = 1;
end

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
    [X Delta] = ode23(@(x,delta) BL_ode(x,delta,Qt_b,dUedx_b,s_b,H,nu),Xspan,IC);

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
    [X Delta] = ode23(@(x,delta) BL_ode(x,delta,Qt_t,dUedx_t,s_t,H,nu),Xspan,IC);

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


% %% Calculating turbulent separation
% Cpp = (Cp(2:end) + Cp(1:end-1))/2;
% Cpp = [Cpp(1);Cpp;Cpp(end)];
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
% Cpbar_b = 1 + (Cp_b - 1)*Qinf^2/ub_m^2;
% Cpbar_t = 1 + (Cp_t - 1)*Qinf^2/ut_m^2;
% 
% 
% % Calculating the points of turbulent transition and the velocity at those
% % points.
% xb_t = s_b(bot_trans);
% xt_t = s_t(top_trans);
% 
% ub_t = Qt_b(bot_trans);
% ut_t = Qt_t(top_trans);
% 
% 
% 
% % Calculating the pressure gradients______________________________________
% dCpbardx_b = zeros(length(s_b),1);
% dCpbardx_t = zeros(length(s_t),1);
% 
% dCpbardx_b(1,1) = (Cpbar_b(2) - Cpbar_b(1))./(dels_b(2)/2 + dels_b(1)/2); 
% dCpbardx_b(2:end-1) = (Cpbar_b(3:end) - Cpbar_b(1:end-2))./(dels_b(3:end)/2 + dels_b(2:end-1) + dels_b(1:end-2)/2); 
% dCpbardx_b(end,1) = (Cpbar_b(end) - Cpbar_b(end - 1))./(dels_b(end)/2 + dels_b(end - 1)/2);
% 
% dCpbardx_t(1,1) = (Cpbar_t(2) - Cpbar_t(1))./(dels_t(2)/2 + dels_t(1)/2); 
% dCpbardx_t(2:end-1) = (Cpbar_t(3:end) - Cpbar_t(1:end-2))./(dels_t(3:end)/2 + dels_t(2:end-1) + dels_t(1:end-2)/2); 
% dCpbardx_t(end,1) = (Cpbar_t(end) - Cpbar_t(end - 1))./(dels_t(end)/2 + dels_t(end - 1)/2);
% 
% 
% 
% bot = real(sqrt(Cp_b).*s_b.*dCpbardx_b);
% top = real(-sqrt(Cp_t).*s_t.*dCpbardx_t);
% turbsep_b = find(bot - 0.102 >= 0,1,'first');
% turbsep_t = find(top - 0.102 >= 0,1,'first');
% 
% 
% Stagpts = find((QT(1:end-1).*QT(2:end)) <= 0);
% LE = Npan/2;
% [~,ind] = min(abs(Stagpts - LE*ones(length(Stagpts),1)));
% Stagpt = Stagpts(ind);
% 
% Xpb = XP(Stagpt:-1:1);
% Xpt = XP(Stagpt:end);
% 
% if isempty(turbsep_b)
%     Indsep_b = 1;
% else
%     Indsep_b = find(Xpb >= xc_b(turbsep_b),1,'first');
%     Indsep_b = Stagpt - Indsep_b;
% end
% 
% if isempty(turbsep_t)
%     Indsep_t = Npan + 1;
% else
%     Indsep_t = find(Xpt >= xc_t(turbsep_t),1,'first');
%     Indsep_t = Stagpt + Indsep_t;
% end


%% Gathering data
Cf = [-Cf_b(end:-1:1); Cf_t];
Cf = (imag(Cf) == 0).*Cf;
dL_act = kron(eye(Npan),ones(1,factor))*dL;
dL_pct = kron(1./dL_act,ones(factor,1)).*dL;

dFi = 1/2*rho*Qt.^2.*Cf.*dL_pct;


