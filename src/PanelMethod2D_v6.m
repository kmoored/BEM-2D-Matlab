% Written by Keith Moored, 2/22/12


clear
close all
clc

% profile -memory on
% setpref('profiler','showJitLines',1);

% function [Ct_avg,Tnet_avg,D_avg,Cl_avg,Cpow_avg,Pow_avg,np,Ec,f,A_TE,rho,c,Ucruise,V_t] = PanelMethod2D_v6(A_bp,M_star,f,A_c,c,Qinf,Npanels,Nstep,Ncyc,Nlump,TEfac,val)
% Npanels
% Nstep
% Ncyc
% Nlump
% Qinf
% f
% A_c
% c
% A_bp
% M_star


%%%% Fixed V %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% St = 0.4;
% A_c = 0.5;
% Qinf = 5;
% M_star = 1/5; 
% A_bp = 1;
% c = 1;
% Npanels = 150;
% Nstep = 150;
% Ncyc = 8;
% Nlump = 4;
% val = 2.5e-4;
% TEfac = 1;




%%%% Free-Swimming %%%%%%%%%%%%%%%%%%%%%%%%%
f = 0.25;
A_c = 0.125;
Qinf = 0.01;
c = 1;
A_bp = 1;
M_star = 1/5;
Npanels = 150;
Nstep = 150;
Ncyc = 15;
Nlump = 4;
val = 2.5e-4;
TEfac = 1;


%% Parameters: 

% Fluid parameters
rho = 998;                 % kg/m^3
nu = 1.004e-6;             % m^2/s
Re = Qinf*c/nu



% Geometry parameters
% tmax = 0.15;
% tmax = 5/16*(0.0254);       % Cylinder diameter
% c = 0.12;
tmax = 0.1*c;       % Cylinder diameter
BL = 1;
% c_D = 10;
% c = c_D*tmax;               % meters
b = 2*c;                    % meters
d_c = 1/2;                    % Distance from the ground in chords
D = d_c*c;                  % Distance from the ground in m

% Paneling parameters
LEfac = 1;
% TEfac = 1;
kappa = 8;
CircRed = 1;

% 2P Kin: f = 1-2, 2S Kin: f = 2.5, h_c = 1 c = 1/10, phi = -pi/2,
% alpha_max = 30.

% Kinematic parameters
% A_c = 0.1667;
% A_c = 0.25;
% h_c = 1/2*A_c;
h_c = 0;
A_TE = 1/2*A_c*c;
% alpha_max = 0;
alpha_max = asin(1/2*A_c)
alpha = 0*(pi/180);        % Enter in degrees, converts to radians.
% f = St*Qinf/(A_c*c)

phi = 0; %pi-(8/10)*pi/2;      % Phase offset of actuation signal.
% M = 1000;

% Virtual body properties
% M_star = 2;
Beta = 1;
C_1 = 1;
M = M_star*rho*(c*b)*(A_c*c)
% A_bp = 50;
A_bod = A_bp*(c*b);
Cd_bod = 0.015;
Ct_prop = 2.25;

% Free-swimming speed and number of cycles to reach steady-state.
U_swim = f*(A_c*c)*sqrt(2*Beta/A_bp*(Ct_prop/Cd_bod))
% % U_swim = (f*(A_c*c))^(4/3)*((Ct_prop/Cd_bod)*(1/A_bp))^(2/3)*(BL/nu)^(1/3)
N_a = C_1*M_star/Ct_prop*sqrt(2/A_bp*(Ct_prop/Cd_bod))
% 
Qinf = 0.01*U_swim;             % m/s 

% Settings
LES = 0;                                        % 0 - off, 1 - on. 
SC = 1;                                         % 0 - off, 1 - on. 
Rollup = 1;                                     % 0 - off, 1 - on.
Flowfield = 0;                                  % 0 - off, 1 - on.
wakeinf = 1;                                    % 0 - off, 1 - on.
bod = 0;                                        % 0 - off, 1 - on.
grd = 1;                                        % 0 - off, 1 - on.
free = 1;                                       % 0 - off, 1 - on.
SurfaceColor = 0;                               % 0 - shear stress, 1 - pressure
ViscDrag = 1;                                   % 0 - off, 1 - on.
PlotTimeStepFig = 1;                            % 0 - off, 1 - on.
PlotVelFig = 1;                                 % 0 - off, 1 - on.
SaveData = 0;                                   % 0 - off, 1 - on.
BodDrag = 0;                                    % 0 - off, 1 - on.

% Cuttoff radii
% val = 1e-3;
ep = val*c;
epSC = ep;
epBod = ep;
epB = val*c;
CptInval = 0.25;



%% Calculations
% h_c = 2*pi*St/f_red                     % Heave to chord ratio
% St = 
% A_TE = h_c*c;
% Qinf = Re*nu/c                       % m/s
Re = c*Qinf/nu
% f = 100;
% f = St*Qinf/(A_c*c)
% Ncyc = max([6 round(6*c/(Qinf/f))])
% f = 1;
% f = Qinf*f_red/2/pi/c                      % Hz
AoA_max = abs(atan2(pi*f*A_c*c,Qinf) + alpha_max)*180/pi
St = f*2*A_TE/Qinf
PressureNum = 2*pi*f^2*c*A_TE/Qinf^2
if f == 0
    delT = 1/0.001/Nstep           % s
else
    delT = 1/f/Nstep           % s
end

NonDimTimeStep = delT*Qinf/c

% Free-stream velocity
Uinf = Qinf*cos(alpha);
Winf = Qinf*sin(alpha);

%% Opening txt file for data storage and saving parameters
if SaveData == 1
    
    savefilename = ['_Heave_k',num2str(k),...
            '_A_c',num2str(A_c),...
            '_Qinf',num2str(Qinf)];
%     savefilename = ['_Pitch_Mstar',num2str(M_star),...
%         '_A_bp',num2str(A_bp),...
%         '_A_c',num2str(A_c),...
%         '_f',num2str(f)];
        
    save(['FixedV/Parameters',savefilename,'.mat'],'-v7.3')

    fid_Data = fopen(['FixedV/Data',savefilename,'.txt'],'w');
    fid_PanelProp = fopen(['FixedV/PanelProp',savefilename,'.txt'],'w');
    fid_WakeProp = fopen(['FixedV/WakeProp',savefilename,'.txt'],'w');
end

%% Initializing Matrices
xw = zeros(1,Nlump*Nstep+1);
zw = zeros(1,Nlump*Nstep+1);
xl = zeros(1,2);
zl = zeros(1,2);

vtTE = zeros(1,2);
vnTE = zeros(1,2);
vtw = zeros(Nlump*Nstep,2);
vnw = zeros(Nlump*Nstep,2);
vtlw = zeros(1,2);
vnlw = zeros(1,2);

sigma = zeros(Npanels,1);
mu = zeros(Npanels,3);
muTE = zeros(1,2);
muW = zeros(1,Nlump*Nstep);
muLump = zeros(1,2);
Gamma = zeros(1,1);
GammaW = zeros(1,Nlump*Nstep+1);
GammaTot = 0;

Qp = zeros(Npanels,1);

Cw = zeros(Npanels,Nlump*Nstep);
Clw = zeros(Npanels,1);
dPhi_wd = zeros(Nlump*Nstep,Npanels);
dPhi_wd_2 = zeros(Nlump*Nstep,Npanels);

dFshear = zeros(Npanels,1);

fighand = zeros(Ncyc*Nstep+1,1);
wakeInd = 0;
Tnet_avg = 0;
D_avg = 0;
Cl_avg = 0;
Ct_avg = 0;
Pow_avg = 0;
Cpow_avg = 0;
np = 0;
Ec = 0;

% Initializing matrices for the LES calcs.
if LES == 1
    xpbod = zeros(Npanels + 1,Ncyc*Nstep + 1);
    zpbod = zeros(Npanels + 1,Ncyc*Nstep + 1);

    xpt_LES = zeros(2,Ncyc*Nstep);
    zpt_LES = zeros(2,Ncyc*Nstep);
    xpb_LES = zeros(2,Ncyc*Nstep);
    zpb_LES = zeros(2,Ncyc*Nstep);
    
    vtLEt = zeros(Ncyc*Nstep,2);
    vnLEt = zeros(Ncyc*Nstep,2);
    vtLEb = zeros(Ncyc*Nstep,2);
    vnLEb = zeros(Ncyc*Nstep,2);
    
    Vshear_t = zeros(2,Ncyc*Nstep+1);
    Vshear_b = zeros(2,Ncyc*Nstep+1);
    
    muLEt = zeros(1,Ncyc*Nstep);
    muLEb = zeros(1,Ncyc*Nstep);
    
    CLEwt = zeros(Npanels,Ncyc*Nstep-1);
    CLEwb = zeros(Npanels,Ncyc*Nstep-1);
    
    Sep_t_Store = ones(1,Ncyc*Nstep+1);
    Sep_b_Store = ones(1,Ncyc*Nstep+1);
    Sep_t = 0;
    Sep_b = 0;
else
    xpbod = [];
    zpbod = [];
    
    xpt_LES = [];
    zpt_LES = [];
    xpb_LES = [];
    zpb_LES = [];
    
    vtLEt = [];
    vnLEt = [];
    vtLEb = [];
    vnLEb = [];
    
    Vshear_t = [];
    Vshear_b = [];
    
    muLEt = [];
    muLEb = [];
    
    CLEwt = [];
    CLEwb = [];
    
    Sep_t_Store = [];
    Sep_b_Store = [];
    Sep_t = 0;
    Sep_b = 0;
end


%% Geometry
% Specifying the surface shape function (Van de Vooren airfoil) and panel 
% corner points.
% [xpt,zpt,xpb,zpb,a,k,epsilon,theta1,theta2,theta,r1,r2] = VanDeVooren(c,Npanels+1,tmax/2);
[xpt,zpt,xpb,zpb] = TearDropShape(c,Npanels+1,tmax);
xp = [xpb;xpt(2:end)];
xp_0 = xp - min(xp);
zp_0 = [zpb(1:end-1)-zpb(1);0;zpt(2:end)-zpt(end)];

xp = xp_0;
zp = zp_0;


[~,Indsep_tmax] = min(abs(xp(Npanels/2 + 1:end)  - 0.05*c)); 
Indsep_tmax = Indsep_tmax + Npanels/2;
[~,Indsep_bmax] = min(abs(xp(1:Npanels/2 + 1) - 0.05*c));

% Locations of collocation points (mid-panel on the surface).
xc_0 = (xp(2:end) + xp(1:end-1))/2;
zc_0 = (zp(2:end) + zp(1:end-1))/2;
zcval = zc_0(Npanels,1);

xc = xc_0;
zc = zc_0;

% Normal vectors for the panels at collocation points.  Each row is a new
% collocation point.
vttemp = [(xp(2:end,1) - xp(1:end-1,1)) (zp(2:end,1) - zp(1:end-1,1))];
vt(:,:,1) = diag(1./sqrt(vttemp(:,1).^2 + vttemp(:,2).^2))*vttemp;
vn(:,:,1) = [-vt(:,2,1) vt(:,1,1)];


%% Flight trajectory: Translation

if free == 1
    Q0 = kron(-[Uinf 0]',ones(1,Ncyc*Nstep+1));
else
    Q0 = kron(-[Uinf Winf]',ones(1,Ncyc*Nstep+1));
end
x_b = zeros(1,3);
z_b = D*ones(1,3);
a_b = 0;

%% Time-stepping solution
scrsz = get(0,'ScreenSize');
% wbar = waitbar(0,'Calculating...');
delTime = zeros(Ncyc*Nstep+1,1);


for i_t = 1:Ncyc*Nstep+1
    % Start timer.
    tic
    
    t = (i_t - 1)*delT;
    
    if i_t > 1
        muW(1,2:end) = muW(1,1:end-1);
        muW(1,1) = muTE(1,1);
    end 
    
    x_b(2:end) = x_b(1:end-1);
    z_b(2:end) = z_b(1:end-1);
    
    % Calculating shift in panel positions with the swimming velocity. 
    delxp = x_b(2)*ones(Npanels+1,1);
    delzp = z_b(2)*ones(Npanels+1,1);
    delxc = x_b(2)*ones(Npanels,1);
    delzc = z_b(2)*ones(Npanels,1);
    
    if LES == 1 && t > 0
        
        
        % Translating all body panels by Q0 without the kinematics update.
        xp_nokin = xpbod(:,i_t-1) + delxp;
        zp_nokin = zpbod(:,i_t-1) + delzp;
        
%         if i_t > 2
%         close(fighand)
%     end
%     % Plotting airfoil with LE vortex sheets and the TE vortex sheet.
%                 fighand = figure;
%                 hold on
%                 axis equal
%                 plot(xp_nokin,zp_nokin,'-k')
%                 plot(xTE(:,i_t-1),zTE(:,i_t-1),'.-b','linewidth',2)
% 
%                 if t > 0
%                     plot(xpt_LES(1,Ncyc*Nstep + 2 - i_t),zpt_LES(1,Ncyc*Nstep + 2 - i_t),'og','linewidth',2)
%                     plot(xpb_LES(1,Ncyc*Nstep + 2 - i_t),zpb_LES(1,Ncyc*Nstep + 2 - i_t),'or','linewidth',2)
% 
%                     for j = i_t-1:-1:2
%             %             if muLEt(j) == 0
%             %             else
%                             plot(xpt_LES(:,Ncyc*Nstep + 1 - j),zpt_LES(:,Ncyc*Nstep + 1 - j),'.-g','linewidth',1)
%             %             end
%                     end
% 
%                     for j = i_t-1:-1:2
%             %             if muLEb(j) == 0   
%             %             else
%                             plot(xpb_LES(:,Ncyc*Nstep + 1 - j),zpb_LES(:,Ncyc*Nstep + 1 - j),'.-r','linewidth',1)
%             %             end
%                     end
% 
%                     plot(xp(Stagpt,i_t-1),zp(Stagpt,i_t-1),'xk','linewidth',2)
%                     plot(xw(Ncyc*Nstep + 3 - i_t:end),zw(Ncyc*Nstep + 3 - i_t:end),'.-b')
%                 end
% 
%                 xlabel('$$x$$','interpreter','latex','fontsize',12,'fontname','TimesNewRoman')
%                 ylabel('$$z$$','interpreter','latex','fontsize',12,'fontname','TimesNewRoman')

%         
%         % Translating the trailing edge panel.
%         xp_avg = xp_nokin(end-1)/2 + xp_nokin(2)/2;
%         zp_avg = zp_nokin(end-1)/2 + zp_nokin(2)/2;
%         TEvec = [xp_nokin(end) - xp_avg; zp_nokin(end) - zp_avg]/norm([xp_nokin(end) - xp_avg; zp_nokin(end) - zp_avg]);
%         xTE(:,i_t) = [xp_nokin(end); xp_nokin(end) + TEvec(1)*Qinf*delT*TEfac];
%         zTE(:,i_t) = [zp_nokin(end); zp_nokin(end) + TEvec(2)*Qinf*delT*TEfac];
% % 
%         
%         % Shedding a wake panel after first time step.
%         % Calculating the downstream position of the newly shed wake panel
%         xw(1,Ncyc*Nstep + 3 - i_t) = xTE(2,i_t-1);
%         zw(1,Ncyc*Nstep + 3 - i_t) = zTE(2,i_t-1);
%         
%         % Calculating the upstream position of the newly shed wake panel
%         xw(1,Ncyc*Nstep + 2 - i_t) = xTE(2,i_t);
%         zw(1,Ncyc*Nstep + 2 - i_t) = zTE(2,i_t);
%         
%         vttemp = [(xw(1,Ncyc*Nstep + 3 - i_t:end) - xw(1,Ncyc*Nstep + 2 - i_t:end-1))' (zw(1,Ncyc*Nstep + 3 - i_t:end) - zw(1,Ncyc*Nstep + 2 - i_t:end-1))'];
%         vtw(Ncyc*Nstep + 2 - i_t:end,:) = diag(1./sqrt(vttemp(:,1).^2 + vttemp(:,2).^2))*vttemp;
%         vnw(Ncyc*Nstep + 2 - i_t:end,:) = [-vtw(Ncyc*Nstep + 2 - i_t:end,2) vtw(Ncyc*Nstep + 2 - i_t:end,1)]; 


        % Normal vector at separation point
        if Indsep_t == Npanels + 1
            vnLESt = vn(1,:,i_t-1)/2 + vn(end,:,i_t-1)/2;
            vtLESt = vt(1,:,i_t-1)/2 + vt(end,:,i_t-1)/2;
            Q_LESt = Qt(1,i_t-1)/2 + Qt(end,i_t-1)/2;
        else
            vnLESt = vn(Indsep_t,:,i_t-1)/2 + vn(Indsep_t-1,:,i_t-1)/2;
            vtLESt = vt(Indsep_t,:,i_t-1)/2 + vt(Indsep_t-1,:,i_t-1)/2;
            Q_LESt = Qt(Indsep_t,i_t-1)/2 + Qt(Indsep_t-1,i_t-1)/2;
        end
        
        if Indsep_b == 1
            vnLESb = vn(1,:,i_t-1)/2 + vn(end,:,i_t-1)/2;
            vtLESb = vt(1,:,i_t-1)/2 + vt(end,:,i_t-1)/2;
            Q_LESb = Qt(1,i_t-1)/2 + Qt(end,i_t-1)/2;
        else
            vnLESb = vn(Indsep_b,:,i_t-1)/2 + vn(Indsep_b-1,:,i_t-1)/2;
            vtLESb = vt(Indsep_b,:,i_t-1)/2 + vt(Indsep_b-1,:,i_t-1)/2;
            Q_LESb = Qt(Indsep_b,i_t-1)/2 + Qt(Indsep_b-1,i_t-1)/2;
        end
        
        %% Defining LE panel.
        LESxt_out = epB*vnLESt(1,1);
        LESzt_out = epB*vnLESt(1,2);
        
        LESxb_out = epB*vnLESb(1,1);
        LESzb_out = epB*vnLESb(1,2);
        
%         LESxt_out = 20*epBod*vnLESt(1,1) + Q_LESt*vtLESt(1,1);
%         LESzt_out = 20*epBod*vnLESt(1,2) + Q_LESt*vtLESt(1,2);
%         
%         LESxb_out = 20*epBod*vnLESb(1,1) + Q_LESb*vtLESb(1,1);
%         LESzb_out = 20*epBod*vnLESb(1,2) + Q_LESb*vtLESb(1,2);

        LEt_vec = [LESxt_out LESzt_out]; 
        LEb_vec = [LESxb_out LESzb_out]; 
        
        xpt_LES(:,Ncyc*Nstep + 2 - i_t) = [xp_nokin(Indsep_t); xp_nokin(Indsep_t) + LEfac*LEt_vec(1)*delT];
        zpt_LES(:,Ncyc*Nstep + 2 - i_t) = [zp_nokin(Indsep_t); zp_nokin(Indsep_t) + LEfac*LEt_vec(2)*delT];

        xpb_LES(:,Ncyc*Nstep + 2 - i_t) = [xp_nokin(Indsep_b); xp_nokin(Indsep_b) + LEfac*LEb_vec(1)*delT];
        zpb_LES(:,Ncyc*Nstep + 2 - i_t) = [zp_nokin(Indsep_b); zp_nokin(Indsep_b) + LEfac*LEb_vec(2)*delT];
        
%         xpt_LES(:,Ncyc*Nstep + 2 - i_t) = [xp_nokin(Indsep_t); xp_nokin(Indsep_t) + LEfac*Vshear_t(1,i_t)*delT];
%         zpt_LES(:,Ncyc*Nstep + 2 - i_t) = [zp_nokin(Indsep_t); zp_nokin(Indsep_t) + LEfac*Vshear_t(2,i_t)*delT];
% 
%         xpb_LES(:,Ncyc*Nstep + 2 - i_t) = [xp_nokin(Indsep_b); xp_nokin(Indsep_b) + LEfac*Vshear_b(1,i_t)*delT];
%         zpb_LES(:,Ncyc*Nstep + 2 - i_t) = [zp_nokin(Indsep_b); zp_nokin(Indsep_b) + LEfac*Vshear_b(2,i_t)*delT];
        
        
        if Sep_t == 0 && i_t > 2
            xpt_LES(1,Ncyc*Nstep + 3 - i_t) = xpt_LES(2,Ncyc*Nstep + 2 - i_t);
            zpt_LES(1,Ncyc*Nstep + 3 - i_t) = zpt_LES(2,Ncyc*Nstep + 2 - i_t);
        else
        end
        
        if Sep_b == 0 && i_t > 2
            xpb_LES(1,Ncyc*Nstep + 3 - i_t) = xpb_LES(2,Ncyc*Nstep + 2 - i_t);
            zpb_LES(1,Ncyc*Nstep + 3 - i_t) = zpb_LES(2,Ncyc*Nstep + 2 - i_t);
        else
        end
            
%         xpt_LES(1,Ncyc*Nstep + 3 - i_t:Ncyc*Nstep + 4 - i_t) = [xp_nokin(Indsep_t), xp_nokin(Indsep_t) + LEfac*LEt_vec(1)*delT];
%         zpt_LES(1,Ncyc*Nstep + 3 - i_t:Ncyc*Nstep + 4 - i_t) = [zp_nokin(Indsep_t), zp_nokin(Indsep_t) + LEfac*LEt_vec(2)*delT];
% 
%         xpb_LES(1,Ncyc*Nstep + 3 - i_t:Ncyc*Nstep + 4 - i_t) = [xp_nokin(Indsep_b), xp_nokin(Indsep_b) + LEfac*LEb_vec(1)*delT];
%         zpb_LES(1,Ncyc*Nstep + 3 - i_t:Ncyc*Nstep + 4 - i_t) = [zp_nokin(Indsep_b), zp_nokin(Indsep_b) + LEfac*LEb_vec(2)*delT];
%         
        
 


        
        %% Translating all free wake points by body trajectory
        if i_t > 2
            if Sep_t == 1
                % Translating free wake points.
                xpt_LES(:,Ncyc*Nstep + 3 - i_t:end) = xpt_LES(:,Ncyc*Nstep + 3 - i_t:end) + Q0(1,i_t-1)*delT;
                zpt_LES(:,Ncyc*Nstep + 3 - i_t:end) = zpt_LES(:,Ncyc*Nstep + 3 - i_t:end) + Q0(2,i_t-1)*delT;

                % Applying fencing scheme to free wake points of top LES sheet. 
                [xstar_t1,zstar_t1,fence1] = fencing(c,0*Vp(:,:,i_t-1),xpt_LES(1,Ncyc*Nstep + 3 - i_t:end),zpt_LES(1,Ncyc*Nstep + 3 - i_t:end),-Q0(1,i_t)*ones(1,i_t-2),-Q0(2,i_t)*ones(1,i_t-2),xp_nokin,zp_nokin,vn(:,:,i_t-1),x_b(1,i_t),z_b(1,i_t),delT,epB);           
                [xstar_t2,zstar_t2,fence2] = fencing(c,0*Vp(:,:,i_t-1),xpt_LES(2,Ncyc*Nstep + 3 - i_t:end),zpt_LES(2,Ncyc*Nstep + 3 - i_t:end),-Q0(1,i_t)*ones(1,i_t-2),-Q0(2,i_t)*ones(1,i_t-2),xp_nokin,zp_nokin,vn(:,:,i_t-1),x_b(1,i_t),z_b(1,i_t),delT,epB);           
                xstar_t = [xstar_t1; xstar_t2];
                zstar_t = [zstar_t1; zstar_t2];
                fence = [fence1; fence2];

                xstar_t = xstar_t.*fence + (xstar_t - Q0(1,i_t-1)*delT).*(~fence);
                zstar_t = zstar_t.*fence + (zstar_t - Q0(2,i_t-1)*delT).*(~fence);

            else
                % Translating free wake points (p1).
                xpt_LES(1,Ncyc*Nstep + 4 - i_t:end) = xpt_LES(1,Ncyc*Nstep + 4 - i_t:end) + Q0(1,i_t-1)*delT;
                zpt_LES(1,Ncyc*Nstep + 4 - i_t:end) = zpt_LES(1,Ncyc*Nstep + 4 - i_t:end) + Q0(2,i_t-1)*delT;

                % Translating free wake points (p2).
                xpt_LES(2,Ncyc*Nstep + 3 - i_t:end) = xpt_LES(2,Ncyc*Nstep + 3 - i_t:end) + Q0(1,i_t-1)*delT;
                zpt_LES(2,Ncyc*Nstep + 3 - i_t:end) = zpt_LES(2,Ncyc*Nstep + 3 - i_t:end) + Q0(2,i_t-1)*delT;

                % Applying fencing scheme to free wake points of top LES
                % sheet.     
                [xstar_t1,zstar_t1,fence1] = fencing(c,0*Vp(:,:,i_t-1),xpt_LES(1,Ncyc*Nstep + 4 - i_t:end),zpt_LES(1,Ncyc*Nstep + 4 - i_t:end),-Q0(1,i_t)*ones(1,i_t-3),-Q0(2,i_t)*ones(1,i_t-3),xp_nokin,zp_nokin,vn(:,:,i_t-1),x_b(1,i_t),z_b(1,i_t),delT,epB);           
                [xstar_t2,zstar_t2,fence2] = fencing(c,0*Vp(:,:,i_t-1),xpt_LES(2,Ncyc*Nstep + 3 - i_t:end),zpt_LES(2,Ncyc*Nstep + 3 - i_t:end),-Q0(1,i_t)*ones(1,i_t-2),-Q0(2,i_t)*ones(1,i_t-2),xp_nokin,zp_nokin,vn(:,:,i_t-1),x_b(1,i_t),z_b(1,i_t),delT,epB);           
                xstar_t = [xpt_LES(1,Ncyc*Nstep + 3 - i_t) xstar_t1; xstar_t2];
                zstar_t = [zpt_LES(1,Ncyc*Nstep + 3 - i_t) zstar_t1; zstar_t2];
                fence = [1 fence1; fence2];

                xstar_t = xstar_t.*fence + (xstar_t - Q0(1,i_t-1)*delT).*(~fence);
                zstar_t = zstar_t.*fence + (zstar_t - Q0(2,i_t-1)*delT).*(~fence);

            end

            if Sep_b == 1
                % Translating free wake points.
                xpb_LES(:,Ncyc*Nstep + 3 - i_t:end) = xpb_LES(:,Ncyc*Nstep + 3 - i_t:end) + Q0(1,i_t-1)*delT;
                zpb_LES(:,Ncyc*Nstep + 3 - i_t:end) = zpb_LES(:,Ncyc*Nstep + 3 - i_t:end) + Q0(2,i_t-1)*delT;

                 % Closing previous figure 
%                 if i_t > 2
%                     close(fighand)
%                 end

%                 % Plotting airfoil with LE vortex sheets and the TE vortex sheet.
%                 fighand = figure;
%                 hold on
%                 axis equal
%                 plot(xp_nokin,zp_nokin,'-k')
%                 plot(xTE(:,i_t-1),zTE(:,i_t-1),'.-b','linewidth',2)
% 
%                 if t > 0
%                     plot(xpt_LES(1,Ncyc*Nstep + 2 - i_t),zpt_LES(1,Ncyc*Nstep + 2 - i_t),'og','linewidth',2)
%                     plot(xpb_LES(1,Ncyc*Nstep + 2 - i_t),zpb_LES(1,Ncyc*Nstep + 2 - i_t),'or','linewidth',2)
% 
%                     for j = i_t-1:-1:2
%             %             if muLEt(j) == 0
%             %             else
%                             plot(xpt_LES(:,Ncyc*Nstep + 1 - j),zpt_LES(:,Ncyc*Nstep + 1 - j),'.-g','linewidth',1)
%             %             end
%                     end
% 
%                     for j = i_t-1:-1:2
%             %             if muLEb(j) == 0   
%             %             else
%                             plot(xpb_LES(:,Ncyc*Nstep + 1 - j),zpb_LES(:,Ncyc*Nstep + 1 - j),'.-r','linewidth',1)
%             %             end
%                     end
% 
%                     plot(xp(Stagpt,i_t-1),zp(Stagpt,i_t-1),'xk','linewidth',2)
%                     plot(xw(Ncyc*Nstep + 3 - i_t:end),zw(Ncyc*Nstep + 3 - i_t:end),'.-b')
%                 end
% 
%                 xlabel('$$x$$','interpreter','latex','fontsize',12,'fontname','TimesNewRoman')
%                 ylabel('$$z$$','interpreter','latex','fontsize',12,'fontname','TimesNewRoman')

                
                % Applying fencing scheme to free wake points of bot LES sheet. 
                [xstar_b1,zstar_b1,fence1] = fencing(c,0*Vp(:,:,i_t-1),xpb_LES(1,Ncyc*Nstep + 3 - i_t:end),zpb_LES(1,Ncyc*Nstep + 3 - i_t:end),-Q0(1,i_t)*ones(1,i_t-2),-Q0(2,i_t)*ones(1,i_t-2),xp_nokin,zp_nokin,vn(:,:,i_t-1),x_b(1,i_t),z_b(1,i_t),delT,epB);           
                [xstar_b2,zstar_b2,fence2] = fencing(c,0*Vp(:,:,i_t-1),xpb_LES(2,Ncyc*Nstep + 3 - i_t:end),zpb_LES(2,Ncyc*Nstep + 3 - i_t:end),-Q0(1,i_t)*ones(1,i_t-2),-Q0(2,i_t)*ones(1,i_t-2),xp_nokin,zp_nokin,vn(:,:,i_t-1),x_b(1,i_t),z_b(1,i_t),delT,epB);           
                xstar_b = [xstar_b1; xstar_b2];
                zstar_b = [zstar_b1; zstar_b2];
                fence = [fence1; fence2];

                xstar_b = xstar_b.*fence + (xstar_b - Q0(1,i_t-1)*delT).*(~fence);
                zstar_b = zstar_b.*fence + (zstar_b - Q0(2,i_t-1)*delT).*(~fence);

            else
                % Translating free wake points (p1).
                xpb_LES(1,Ncyc*Nstep + 4 - i_t:end) = xpb_LES(1,Ncyc*Nstep + 4 - i_t:end) + Q0(1,i_t-1)*delT;
                zpb_LES(1,Ncyc*Nstep + 4 - i_t:end) = zpb_LES(1,Ncyc*Nstep + 4 - i_t:end) + Q0(2,i_t-1)*delT;

                % Translating free wake points (p2).
                xpb_LES(2,Ncyc*Nstep + 3 - i_t:end) = xpb_LES(2,Ncyc*Nstep + 3 - i_t:end) + Q0(1,i_t-1)*delT;
                zpb_LES(2,Ncyc*Nstep + 3 - i_t:end) = zpb_LES(2,Ncyc*Nstep + 3 - i_t:end) + Q0(2,i_t-1)*delT;

                % Applying fencing scheme to free wake points of bot LES sheet. 
                [xstar_b1,zstar_b1,fence1] = fencing(c,0*Vp(:,:,i_t-1),xpb_LES(1,Ncyc*Nstep + 4 - i_t:end),zpb_LES(1,Ncyc*Nstep + 4 - i_t:end),-Q0(1,i_t-1)*ones(1,i_t-3),-Q0(2,i_t-1)*ones(1,i_t-3),xp_nokin,zp_nokin,vn(:,:,i_t-1),x_b(1,i_t),z_b(1,i_t),delT,epB);           
                [xstar_b2,zstar_b2,fence2] = fencing(c,0*Vp(:,:,i_t-1),xpb_LES(2,Ncyc*Nstep + 3 - i_t:end),zpb_LES(2,Ncyc*Nstep + 3 - i_t:end),-Q0(1,i_t-1)*ones(1,i_t-2),-Q0(2,i_t-1)*ones(1,i_t-2),xp_nokin,zp_nokin,vn(:,:,i_t-1),x_b(1,i_t),z_b(1,i_t),delT,epB);           
                xstar_b = [xpb_LES(1,Ncyc*Nstep + 3 - i_t) xstar_b1; xstar_b2];
                zstar_b = [zpb_LES(1,Ncyc*Nstep + 3 - i_t) zstar_b1; zstar_b2];
                fence = [1 fence1; fence2];

                xstar_b = xstar_b.*fence + (xstar_b - Q0(1,i_t-1)*delT).*(~fence);
                zstar_b = zstar_b.*fence + (zstar_b - Q0(2,i_t-1)*delT).*(~fence);
            end
        

%             xpt_LES(1,Ncyc*Nstep + 5 - i_t:end) = xpt_LES(1,Ncyc*Nstep + 5 - i_t:end) + Q0(1,i_t-1)*delT;
%             zpt_LES(1,Ncyc*Nstep + 5 - i_t:end) = zpt_LES(1,Ncyc*Nstep + 5 - i_t:end) + Q0(2,i_t-1)*delT;
%             
%             xpb_LES(1,Ncyc*Nstep + 5 - i_t:end) = xpb_LES(1,Ncyc*Nstep + 5 - i_t:end) + Q0(1,i_t-1)*delT;
%             zpb_LES(1,Ncyc*Nstep + 5 - i_t:end) = zpb_LES(1,Ncyc*Nstep + 5 - i_t:end) + Q0(2,i_t-1)*delT;

        
%         % Applying fencing scheme to top LES sheet     
%     
%         xstar_t = xstar_t.*fence + (xstar_t - Q0(1,i_t-1)*delT).*(~fence);
%         zstar_t = zstar_t.*fence + (zstar_t - Q0(2,i_t-1)*delT).*(~fence);

%         % Applying fencing scheme to bot LES sheet
%         [xstar_b,zstar_b,fence] = fencing(c,Vp(:,:,i_t-1),xpb_LES(Ncyc*Nstep + 5 - i_t:end),zpb_LES(Ncyc*Nstep + 5 - i_t:end),-Q0(1,i_t-1)*ones(1,i_t-2),-Q0(2,i_t-1)*ones(1,i_t-2),xp_nokin,zp_nokin,vn(:,:,i_t-1),x_b(1,i_t),z_b(1,i_t),delT,epB);           
% 
%         xstar_b = xstar_b.*fence + (xstar_b - Q0(1,i_t-1)*delT).*(~fence);
%         zstar_b = zstar_b.*fence + (zstar_b - Q0(2,i_t-1)*delT).*(~fence);

            % Applying fencing to LES wake points.
            xpt_LES(:,Ncyc*Nstep + 3 - i_t:end) = xstar_t;
            zpt_LES(:,Ncyc*Nstep + 3 - i_t:end) = zstar_t;

            xpb_LES(:,Ncyc*Nstep + 3 - i_t:end) = xstar_b;
            zpb_LES(:,Ncyc*Nstep + 3 - i_t:end) = zstar_b;
        end
        
        % Define LES panel tangential and normal vectors.
        vttemp = [(xpt_LES(2,Ncyc*Nstep + 2 - i_t:end) - xpt_LES(1,Ncyc*Nstep + 2 - i_t:end))' (zpt_LES(2,Ncyc*Nstep + 2 - i_t:end) - zpt_LES(1,Ncyc*Nstep + 2 - i_t:end))'];
        vtLEt(Ncyc*Nstep + 2 - i_t:end,:) = diag(1./sqrt(vttemp(:,1).^2 + vttemp(:,2).^2))*vttemp;
        vnLEt(Ncyc*Nstep + 2 - i_t:end,:) = [-vtLEt(Ncyc*Nstep + 2 - i_t:end,2) vtLEt(Ncyc*Nstep + 2 - i_t:end,1)];

        vttemp = [(xpb_LES(2,Ncyc*Nstep + 2 - i_t:end) - xpb_LES(1,Ncyc*Nstep + 2 - i_t:end))' (zpb_LES(2,Ncyc*Nstep + 2 - i_t:end) - zpb_LES(1,Ncyc*Nstep + 2 - i_t:end))'];
        vtLEb(Ncyc*Nstep + 2 - i_t:end,:) = diag(1./sqrt(vttemp(:,1).^2 + vttemp(:,2).^2))*vttemp;
        vnLEb(Ncyc*Nstep + 2 - i_t:end,:) = [-vtLEb(Ncyc*Nstep + 2 - i_t:end,2) vtLEb(Ncyc*Nstep + 2 - i_t:end,1)]; 
    end  
    
    
   
    
    %% Updating the kinematics.    
    
    % Ramping up the amplitude of motion.
%     if i_t < Nstep+1
%         alpha_max_ru = alpha_max*t*f;
%     else
        alpha_max_ru = alpha_max;
%     end
        
%     [xp,zp,Vp] = Kinematics_Pitch2D(xp_0,zp_0,alpha_max_ru,f,t,phi,tmax);
    [xp,zp,Vp] = Kinematics_HeavePitch2D(xp_0,zp_0,alpha_max_ru,h_c,f,t,phi);
    xpbod = xp;
    zpbod = zp;
    Vc = 1/2*Vp(2:end,:) + 1/2*Vp(1:end-1,:);
    
%     if LES == 1 && i_t > 1
%         Vbod(:,1,i_t) = (xpbod(:,i_t) - xpbod(:,i_t-1))/delT;
%         Vbod(:,2,i_t) = (zpbod(:,i_t) - zpbod(:,i_t-1))/delT;
%     else
%         Vbod = Vp;
%     end
    
    % Calculating locations of collocation points (mid-panel on the surface).
    xc = (xp(2:end) + xp(1:end-1))/2;
    zc = (zp(2:end) + zp(1:end-1))/2;
    
    % Superposing the kinematics and swimming translations.
    xp = xp + delxp;
    zp = zp + delzp;
    xc = xc + delxc;
    zc = zc + delzc;
   
    % Normal vectors for the panels at collocation points.  Each row is a new
    % collocation point.
    vttemp = [(xp(2:end) - xp(1:end-1)) (zp(2:end) - zp(1:end-1))];
    vt = diag(1./sqrt(vttemp(:,1).^2 + vttemp(:,2).^2))*vttemp;
    vn = [-vt(:,2) vt(:,1)];
    
    %% Applying fencing after updating the kinematics
    if LES == 1 && i_t > 1

        % Normal vector at separation point
        if Indsep_t == Npanels + 1
            vnLESt = vn(1,:,i_t)/2 + vn(end,:,i_t)/2;
            vtLESt = vt(1,:,i_t)/2 + vt(end,:,i_t)/2;
            Q_LESt = Qt(1,i_t-1)/2 + Qt(end,i_t-1)/2;
        else
            vnLESt = vn(Indsep_t,:,i_t)/2 + vn(Indsep_t-1,:,i_t)/2;
            vtLESt = vt(Indsep_t,:,i_t)/2 + vt(Indsep_t-1,:,i_t)/2;
            Q_LESt = Qt(Indsep_t,i_t-1)/2 + Qt(Indsep_t-1,i_t-1)/2;
        end
        
        if Indsep_b == 1
            vnLESb = vn(1,:,i_t)/2 + vn(end,:,i_t)/2;
            vtLESb = vt(1,:,i_t)/2 + vt(end,:,i_t)/2;
            Q_LESb = Qt(1,i_t-1)/2 + Qt(end,i_t-1)/2;
        else
            vnLESb = vn(Indsep_b,:,i_t)/2 + vn(Indsep_b-1,:,i_t)/2;
            vtLESb = vt(Indsep_b,:,i_t)/2 + vt(Indsep_b-1,:,i_t)/2;
            Q_LESb = Qt(Indsep_b,i_t-1)/2 + Qt(Indsep_b-1,i_t-1)/2;
        end
        
        %% Defining LE panel after kinematics update.
        LESxt_out = epB*vnLESt(1,1);
        LESzt_out = epB*vnLESt(1,2);
        
        LESxb_out = epB*vnLESb(1,1);
        LESzb_out = epB*vnLESb(1,2);
        
%         LESxt_out = 20*epBod*vnLESt(1,1) + Q_LESt*vtLESt(1,1);
%         LESzt_out = 20*epBod*vnLESt(1,2) + Q_LESt*vtLESt(1,2);
%         
%         LESxb_out = 20*epBod*vnLESb(1,1) + Q_LESb*vtLESb(1,1);
%         LESzb_out = 20*epBod*vnLESb(1,2) + Q_LESb*vtLESb(1,2);

        LEt_vec = [LESxt_out LESzt_out];
        LEb_vec = [LESxb_out LESzb_out]; 
        
%         xpt_LES(:,Ncyc*Nstep + 2 - i_t) = [xp(Indsep_t,i_t); xp(Indsep_t,i_t) + LEfac*Vshear_t(1,i_t)*delT];
%         zpt_LES(:,Ncyc*Nstep + 2 - i_t) = [zp(Indsep_t,i_t); zp(Indsep_t,i_t) + LEfac*Vshear_t(2,i_t)*delT];
% 
%         xpb_LES(:,Ncyc*Nstep + 2 - i_t) = [xp(Indsep_b,i_t); xp(Indsep_b,i_t) + LEfac*Vshear_b(1,i_t)*delT];
%         zpb_LES(:,Ncyc*Nstep + 2 - i_t) = [zp(Indsep_b,i_t); zp(Indsep_b,i_t) + LEfac*Vshear_b(2,i_t)*delT];
        
        xpt_LES(:,Ncyc*Nstep + 2 - i_t) = [xp(Indsep_t,i_t); xp(Indsep_t,i_t) + LEfac*LEt_vec(1)*delT];
        zpt_LES(:,Ncyc*Nstep + 2 - i_t) = [zp(Indsep_t,i_t); zp(Indsep_t,i_t) + LEfac*LEt_vec(2)*delT];

        xpb_LES(:,Ncyc*Nstep + 2 - i_t) = [xp(Indsep_b,i_t); xp(Indsep_b,i_t) + LEfac*LEb_vec(1)*delT];
        zpb_LES(:,Ncyc*Nstep + 2 - i_t) = [zp(Indsep_b,i_t); zp(Indsep_b,i_t) + LEfac*LEb_vec(2)*delT];
        
        if Sep_t == 0 && i_t > 2
            xpt_LES(1,Ncyc*Nstep + 3 - i_t) = xpt_LES(2,Ncyc*Nstep + 2 - i_t);
            zpt_LES(1,Ncyc*Nstep + 3 - i_t) = zpt_LES(2,Ncyc*Nstep + 2 - i_t);
        else
        end
        
        if Sep_b == 0 && i_t > 2
            xpb_LES(1,Ncyc*Nstep + 3 - i_t) = xpb_LES(2,Ncyc*Nstep + 2 - i_t);
            zpb_LES(1,Ncyc*Nstep + 3 - i_t) = zpb_LES(2,Ncyc*Nstep + 2 - i_t);
        else
        end
        
%         xpt_LES(1,Ncyc*Nstep + 3 - i_t:Ncyc*Nstep + 4 - i_t) = [xp_nokin(Indsep_t), xp_nokin(Indsep_t) + LEfac*LEt_vec(1)*delT];
%         zpt_LES(1,Ncyc*Nstep + 3 - i_t:Ncyc*Nstep + 4 - i_t) = [zp_nokin(Indsep_t), zp_nokin(Indsep_t) + LEfac*LEt_vec(2)*delT];
% 
%         xpb_LES(1,Ncyc*Nstep + 3 - i_t:Ncyc*Nstep + 4 - i_t) = [xp_nokin(Indsep_b), xp_nokin(Indsep_b) + LEfac*LEb_vec(1)*delT];
%         zpb_LES(1,Ncyc*Nstep + 3 - i_t:Ncyc*Nstep + 4 - i_t) = [zp_nokin(Indsep_b), zp_nokin(Indsep_b) + LEfac*LEb_vec(2)*delT];
%         
   
        %% Applying fencing after updating kinematics
        if i_t > 2
            if Sep_t == 1

                % Applying fencing scheme to free wake points of top LES sheet. 
                [xstar_t1,zstar_t1,~] = fencing(c,Vbod(:,:,i_t),xpt_LES(1,Ncyc*Nstep + 3 - i_t:end),zpt_LES(1,Ncyc*Nstep + 3 - i_t:end),0*ones(1,i_t-2),0*ones(1,i_t-2),xp(:,i_t),zp(:,i_t),vn(:,:,i_t),x_b(1,i_t),z_b(1,i_t),delT,epB);           
                [xstar_t2,zstar_t2,~] = fencing(c,Vbod(:,:,i_t),xpt_LES(2,Ncyc*Nstep + 3 - i_t:end),zpt_LES(2,Ncyc*Nstep + 3 - i_t:end),0*ones(1,i_t-2),0*ones(1,i_t-2),xp(:,i_t),zp(:,i_t),vn(:,:,i_t),x_b(1,i_t),z_b(1,i_t),delT,epB);           
                xstar_t = [xstar_t1; xstar_t2];
                zstar_t = [zstar_t1; zstar_t2];

            else

                % Applying fencing scheme to free wake points of top LES sheet. 
                [xstar_t1,zstar_t1,~] = fencing(c,Vbod(:,:,i_t),xpt_LES(1,Ncyc*Nstep + 4 - i_t:end),zpt_LES(1,Ncyc*Nstep + 4 - i_t:end),0*ones(1,i_t-3),0*ones(1,i_t-3),xp(:,i_t),zp(:,i_t),vn(:,:,i_t),x_b(1,i_t),z_b(1,i_t),delT,epB);           
                [xstar_t2,zstar_t2,~] = fencing(c,Vbod(:,:,i_t),xpt_LES(2,Ncyc*Nstep + 3 - i_t:end),zpt_LES(2,Ncyc*Nstep + 3 - i_t:end),0*ones(1,i_t-2),0*ones(1,i_t-2),xp(:,i_t),zp(:,i_t),vn(:,:,i_t),x_b(1,i_t),z_b(1,i_t),delT,epB);           
                xstar_t = [xpt_LES(1,Ncyc*Nstep + 3 - i_t) xstar_t1; xstar_t2];
                zstar_t = [zpt_LES(1,Ncyc*Nstep + 3 - i_t) zstar_t1; zstar_t2];

            end

            if Sep_b == 1

                % Applying fencing scheme to free wake points of bot LES sheet. 
                [xstar_b1,zstar_b1,~] = fencing(c,Vbod(:,:,i_t),xpb_LES(1,Ncyc*Nstep + 3 - i_t:end),zpb_LES(1,Ncyc*Nstep + 3 - i_t:end),0*ones(1,i_t-2),0*ones(1,i_t-2),xp(:,i_t),zp(:,i_t),vn(:,:,i_t),x_b(1,i_t),z_b(1,i_t),delT,epB);           
                [xstar_b2,zstar_b2,~] = fencing(c,Vbod(:,:,i_t),xpb_LES(2,Ncyc*Nstep + 3 - i_t:end),zpb_LES(2,Ncyc*Nstep + 3 - i_t:end),0*ones(1,i_t-2),0*ones(1,i_t-2),xp(:,i_t),zp(:,i_t),vn(:,:,i_t),x_b(1,i_t),z_b(1,i_t),delT,epB);           
                xstar_b = [xstar_b1; xstar_b2];
                zstar_b = [zstar_b1; zstar_b2];

            else

                % Applying fencing scheme to free wake points of bot LES sheet. 
                [xstar_b1,zstar_b1,~] = fencing(c,Vbod(:,:,i_t),xpb_LES(1,Ncyc*Nstep + 4 - i_t:end),zpb_LES(1,Ncyc*Nstep + 4 - i_t:end),0*ones(1,i_t-3),0*ones(1,i_t-3),xp(:,i_t),zp(:,i_t),vn(:,:,i_t),x_b(1,i_t),z_b(1,i_t),delT,epB);           
                [xstar_b2,zstar_b2,~] = fencing(c,Vbod(:,:,i_t),xpb_LES(2,Ncyc*Nstep + 3 - i_t:end),zpb_LES(2,Ncyc*Nstep + 3 - i_t:end),0*ones(1,i_t-2),0*ones(1,i_t-2),xp(:,i_t),zp(:,i_t),vn(:,:,i_t),x_b(1,i_t),z_b(1,i_t),delT,epB);           
                xstar_b = [xpb_LES(1,Ncyc*Nstep + 3 - i_t) xstar_b1; xstar_b2];
                zstar_b = [zpb_LES(1,Ncyc*Nstep + 3 - i_t) zstar_b1; zstar_b2];

            end

            % Applying fencing to LES wake points.
            xpt_LES(:,Ncyc*Nstep + 3 - i_t:end) = xstar_t;
            zpt_LES(:,Ncyc*Nstep + 3 - i_t:end) = zstar_t;

            xpb_LES(:,Ncyc*Nstep + 3 - i_t:end) = xstar_b;
            zpb_LES(:,Ncyc*Nstep + 3 - i_t:end) = zstar_b;
        end
        
        % Define LES panel tangential and normal vectors.
        vttemp = [(xpt_LES(2,Ncyc*Nstep + 2 - i_t:end) - xpt_LES(1,Ncyc*Nstep + 2 - i_t:end))' (zpt_LES(2,Ncyc*Nstep + 2 - i_t:end) - zpt_LES(1,Ncyc*Nstep + 2 - i_t:end))'];
        vtLEt(Ncyc*Nstep + 2 - i_t:end,:) = diag(1./sqrt(vttemp(:,1).^2 + vttemp(:,2).^2))*vttemp;
        vnLEt(Ncyc*Nstep + 2 - i_t:end,:) = [-vtLEt(Ncyc*Nstep + 2 - i_t:end,2) vtLEt(Ncyc*Nstep + 2 - i_t:end,1)]; 

        vttemp = [(xpb_LES(2,Ncyc*Nstep + 2 - i_t:end) - xpb_LES(1,Ncyc*Nstep + 2 - i_t:end))' (zpb_LES(2,Ncyc*Nstep + 2 - i_t:end) - zpb_LES(1,Ncyc*Nstep + 2 - i_t:end))'];
        vtLEb(Ncyc*Nstep + 2 - i_t:end,:) = diag(1./sqrt(vttemp(:,1).^2 + vttemp(:,2).^2))*vttemp;
        vnLEb(Ncyc*Nstep + 2 - i_t:end,:) = [-vtLEb(Ncyc*Nstep + 2 - i_t:end,2) vtLEb(Ncyc*Nstep + 2 - i_t:end,1)]; 
    end
   

    
    
    %% Calculating the new trailing edge panels.
    if i_t > 1
        xTE_old = xTE;
        zTE_old = zTE;
    end
    xp_avg = xp(end-1)/2 + xp(2)/2;
    zp_avg = zp(end-1)/2 + zp(2)/2;
    TEvec = [xp(end) - xp_avg; zp(end) - zp_avg]/norm([xp(end) - xp_avg; zp(end) - zp_avg]);
    xTE = [xp(end); xp(end) + TEvec(1)*sqrt(Q0(:,i_t)'*Q0(:,i_t))*delT*TEfac];
    zTE = [zp(end); zp(end) + TEvec(2)*sqrt(Q0(:,i_t)'*Q0(:,i_t))*delT*TEfac];

    dLTE = sqrt((xTE(2)-xTE(1))^2 + (zTE(2)-zTE(1))^2);
    
    vttemp = [(xTE(2) - xTE(1)) (zTE(2) - zTE(1))];
    vtTE = diag(1./sqrt(vttemp(1).^2 + vttemp(2).^2))*vttemp;
    vnTE = [-vtTE(2) vtTE(1)];
    

    %% Shedding a wake panel after first time step.
    if i_t > 1
       
        muLump(2) = muLump(1);
        xl(2) = xl(1);
        zl(2) = zl(1);
        
        % Lumped wake doublet elements 
        if i_t > Nlump*Nstep + 1
            muLump(1) = muLump(2) - GammaW(Nlump*Nstep + 1);
            GammaTot_new = GammaTot + abs(-GammaW(Nlump*Nstep + 1));

            xl(1) = abs(GammaTot/GammaTot_new)*xl(2) + abs(-GammaW(Nlump*Nstep + 1)/GammaTot_new)*xw(Nlump*Nstep + 1);
            zl(1) = abs(GammaTot/GammaTot_new)*zl(2) + abs(-GammaW(Nlump*Nstep + 1)/GammaTot_new)*zw(Nlump*Nstep + 1);

            GammaTot = GammaTot_new;
            
            vttemp = [(xl(1) - xw(Nlump*Nstep))' (zl(1) - zw(Nlump*Nstep))'];
            vtlw = diag(1./sqrt(vttemp(1).^2 + vttemp(2).^2))*vttemp;
            vnlw = [-vtlw(2) vtlw(1)]; 
        end
        
        
        % Shifting wake panel coordinates
        xw(1,2:end) = xw(1,1:end-1);
        zw(1,2:end) = zw(1,1:end-1);
        
        % Calculating the downstream position of the newly shed wake panel
        xw(1,2) = xTE_old(2);
        zw(1,2) = zTE_old(2);
        
        % Calculating the upstream position of the newly shed wake panel
        xw(1,1) = xTE(2);
        zw(1,1) = zTE(2);
        
        vttemp = [(xw(1,2:end) - xw(1,1:end-1))' (zw(1,2:end) - zw(1,1:end-1))'];
        vtw = diag(1./sqrt(vttemp(:,1).^2 + vttemp(:,2).^2))*vttemp;
        vnw = [-vtw(:,2) vtw(:,1)]; 
          
    end
       
    %% Calculating source strengths. 
    Vt = kron(ones(Npanels,1),Q0(:,i_t)') + Vc;
    sigma = sum(vn.*Vt,2);
    
 
    
    %% Moving collocation points inward.
    
    InVal = abs(CptInval*zcval);
    Xc = xc - InVal.*vn(:,1);
    Zc = zc - InVal.*vn(:,2);


    %% Calculating influence coefficients 
    % (a = q*n, q = [u, w]') at each collocation point.  The row denotes the collocation
    % point while the column denotes the doublet element.  Also calculating RHS,
    % RHS = -Qinf*n.
    
    % Influence of the body panels
    Cb = zeros(Npanels,Npanels);
    dPhi_s_2 = zeros(Npanels,Npanels);
    dPhi_d_2 = zeros(Npanels,Npanels);
    
    [dPhi_s,dPhi_d,dL] = Phi(ones(1,Npanels),ones(1,Npanels),Xc',Zc',xp',zp',vt',vn');
    
    if grd == 1
        [dPhi_s_2,dPhi_d_2,dLtemp2] = Phi(ones(1,Npanels),ones(1,Npanels),Xc',Zc',xp',-zp',[vt(:,1) -vt(:,2)]',[vn(:,1) -vn(:,2)]');
    end
    
    B = dPhi_s' + dPhi_s_2';
    Cb(:,1:Npanels) = dPhi_d' + dPhi_d_2';

   
    % Influence of the trailing edge panel
    dPhi_d_2 = zeros(1,Npanels);
    Cte = zeros(Npanels,Npanels);
    
    [~,dPhi_d,~] = Phi(1,0,Xc',Zc',xTE',zTE',vtTE',vnTE');
    
    if grd == 1
        [~,dPhi_d_2,~] = Phi(1,0,Xc',Zc',xTE',-zTE',[vtTE(1,1) -vtTE(1,2)]',[vnTE(1,1) -vnTE(1,2)]');
    end
  
    Cte(:,Npanels) = dPhi_d' + dPhi_d_2';
    Cte(:,1) = -dPhi_d' - dPhi_d_2';
       
    % Influence of the wake panels    
     if i_t > 1   
        wakeInd = min([(i_t-1) Nlump*Nstep]);
%         wakeInd = min([(find(muW == 0,1)-1) Nlump*Nstep]); 
        [~,dPhi_d,~] = Phi(ones(1,length(muW(1:wakeInd))),zeros(1,length(muW(1:wakeInd))),Xc',Zc',xw(1:wakeInd+1),zw(1:wakeInd + 1),vtw(1:wakeInd,:)',vnw(1:wakeInd,:)');
        
        dPhi_d_2 = 0*dPhi_d;
        
        if grd == 1
            [~,dPhi_d_2,~] = Phi(ones(1,length(muW(1:wakeInd))),zeros(1,length(muW(1:wakeInd))),Xc',Zc',xw(1:wakeInd+1),-zw(1:wakeInd+1),[vtw(1:wakeInd,1) -vtw(1:wakeInd,2)]',[vnw(1:wakeInd,1) -vnw(1:wakeInd,2)]');
        end
        
        Cw(:,1:wakeInd) = dPhi_d' + dPhi_d_2';
    end
  
    % Influence of the lumped wake panel    
     if i_t > 1   
        [~,dPhi_d,~] = Phi(1,0,Xc',Zc',[xw(Nlump*Nstep+1) xl(1)],[zw(Nlump*Nstep+1) zl(1)],vtlw',vnlw');
        
        dPhi_d_2 = 0*dPhi_d;
        
        if grd == 1
            [~,dPhi_d_2,~] = Phi(1,0,Xc',Zc',[xw(Nlump*Nstep+1) xl(1)],[-zw(Nlump*Nstep+1) -zl(1)],[vtlw(1,1) -vtlw(1,2)]',[vnlw(1,1) -vnlw(1,2)]');
        end
        
        Clw(:,1) = dPhi_d' + dPhi_d_2';
    end

    if LES == 1 && t > 0
        
        % Influence of the LEtop edge panel
        CLEt = zeros(Npanels,Npanels);
        if Indsep_t == Npanels + 1
        else
            [~,dPhi_d,~] = Phi(1,0,Xc',Zc',xpt_LES(:,Ncyc*Nstep + 2 - i_t)',zpt_LES(:,Ncyc*Nstep + 2 - i_t)',vtLEt(Ncyc*Nstep + 2 - i_t,:)',vnLEt(Ncyc*Nstep + 2 - i_t,:)');
            CLEt(:,Indsep_t-1) = CircRed*dPhi_d';
            CLEt(:,Indsep_t) = -CircRed*dPhi_d';
        end

        % Influence of the LEbot edge panel
        CLEb = zeros(Npanels,Npanels);
        if Indsep_b == 1
        else
            [~,dPhi_d,~] = Phi(1,0,Xc',Zc',xpb_LES(:,Ncyc*Nstep + 2 - i_t)',zpb_LES(:,Ncyc*Nstep + 2 - i_t)',vtLEb(Ncyc*Nstep + 2 - i_t,:)',vnLEb(Ncyc*Nstep + 2 - i_t,:)');
            CLEb(:,Indsep_b) = -CircRed*dPhi_d';
            CLEb(:,Indsep_b-1) = CircRed*dPhi_d';
        end
        
        if i_t > 2
            for j = i_t:-1:3
                % Influence of the LE top wake sheet panels
                [~,dPhi_d,~] = Phi(ones(1,length(muLEt(Ncyc*Nstep + 3 - j))),zeros(1,length(muLEt(Ncyc*Nstep + 3 - j))),Xc',Zc',xpt_LES(:,Ncyc*Nstep + 3 - j)',zpt_LES(:,Ncyc*Nstep + 3 - j)',vtLEt(Ncyc*Nstep + 3 - j,:)',vnLEt(Ncyc*Nstep + 3 - j,:)');
                CLEwt(:,Ncyc*Nstep + 3 - j) = dPhi_d';

                % Influence of the LE top wake sheet panels
                [~,dPhi_d,~] = Phi(ones(1,length(muLEb(Ncyc*Nstep + 3 - j))),zeros(1,length(muLEb(Ncyc*Nstep + 3 - j))),Xc',Zc',xpb_LES(:,Ncyc*Nstep + 3 - j)',zpb_LES(:,Ncyc*Nstep + 3 - j)',vtLEb(Ncyc*Nstep + 3 - j,:)',vnLEb(Ncyc*Nstep + 3 - j,:)');
                CLEwb(:,Ncyc*Nstep + 3 - j) = dPhi_d';
            end
        end
    end   

    %% Constructing full matrix    
    if LES == 1 && i_t > 1
        A = Cb + Cte + CLEt + CLEb;
    else
        A = Cb + Cte;
    end
   
    % Defining Right Hand Side (RHS)
    if LES == 1 && i_t > 2 
        RHS = -B*sigma - wakeinf*Cw*muW' - wakeinf*Clw*muLump(1) - wakeinf*CLEwt*muLEt' - wakeinf*CLEwb*muLEb';
    elseif i_t > 1
        RHS = -B*sigma - wakeinf*Cw*muW' - wakeinf*Clw*muLump(1);
    else
        RHS = -B*sigma;
    end

    % Solving for doublet strengths.
    mu(:,2:end) = mu(:,1:end-1);
    mu(:,1) = A\RHS;
    
    % Solving for TE panel strength.
    muTE(1,2) = muTE(1,1);
    muTE(1,1) = mu(Npanels,1) - mu(1,1);
                  
    % Calculating wake circulation.
    GammaW(2:end) = GammaW(1:end-1);
    if i_t == 1
        GammaW(1,1) = -muTE(1,1);
    else
        GammaW(1,1) = -(muTE(1,1) - muTE(1,2));
    end
    
    % Solving for LE panel strengths.
    if LES == 1 && t > 0
        if Indsep_b == 1
            muLEb(1,Ncyc*Nstep + 2 - i_t) = 0;
        else
            muLEb(1,Ncyc*Nstep + 2 - i_t) = CircRed*(mu(Indsep_b - 1,i_t) - mu(Indsep_b,i_t));
        end
        
        if Indsep_t == Npanels + 1
            muLEt(1,Ncyc*Nstep + 2 - i_t) = 0;
        else
            muLEt(1,Ncyc*Nstep + 2 - i_t) = CircRed*(mu(Indsep_t - 1,i_t) - mu(Indsep_t,i_t));
        end
    end
    

    %% Calculating on-body velocities and pressures. 
    % Pressures are calculated from the unsteady Bernoulli's equation.  
    % The flux through the body should be zero.  
    
    % This is the order of finite difference scheme.
    Order = 2; 
 
    % Local velocity over the body due to the perturbation potential. 
    Qp(1) = SecondOrderForwDiff(mu(:,1),dL,[1, 2, 3]);
    Qp(2) = SecondOrderCenForwDiff(mu(:,1),dL,[1, 2, 3, 4]);
    
%     Qp(3:end-2,i_t) = (mu(2:end-3,i_t)- mu(4:end-1,i_t))./(1/2*dL(2:end-3,i_t) + dL(3:end-2,i_t) + 1/2*dL(4:end-1,i_t));
    Pan = [1:Npanels-4; 2:Npanels-3; 3:Npanels-2; 4:Npanels-1; 5:Npanels];
    for j = 1:Npanels-4  
        Qp(j+2) = CentDiff(mu(:,1),dL,Pan(:,j),Order);
    end

    Qp(Npanels-1) = SecondOrderCenBackDiff(mu(:,1),dL,[(Npanels - 3), (Npanels - 2), (Npanels - 1), Npanels]);
    Qp(Npanels) = SecondOrderBackDiff(mu(:,1),dL,[(Npanels - 2), (Npanels - 1), Npanels]);

    
    if LES == 1 && t > 0
        % Panel velocity calculation for the panel upstream and the panel 
        % downstream of the top separation location 
        if Indsep_t < Npanels  
            Qp(Indsep_t,i_t) = (mu(Indsep_t,i_t) - mu(Indsep_t+1,i_t))/(dL(Indsep_t,i_t) + dL(Indsep_t+1,i_t)); 
            Qp(Indsep_t-1,i_t) = (mu(Indsep_t-2,i_t) - mu(Indsep_t-1,i_t))/(dL(Indsep_t-2,i_t) + dL(Indsep_t-1,i_t));
        end

        % Panel velocity calculation for the panel upstream and the panel 
        % downstream of the top separation location 
        if Indsep_b >  2
            Qp(Indsep_b,i_t) = (mu(Indsep_b,i_t) - mu(Indsep_b+1,i_t))/(dL(Indsep_b,i_t) + dL(Indsep_b+1,i_t)); 
            Qp(Indsep_b-1,i_t) = (mu(Indsep_b-2,i_t) - mu(Indsep_b-1,i_t))/(dL(Indsep_b-2,i_t) + dL(Indsep_b-1,i_t));
        end
    end
    
    % Calculating the total velocity (free-stream plus the perturbation)
    Qt = Qp + sum(-Vt.*vt,2);

    if LES == 1 && t > 0 
        if Indsep_t < Npanels + 1
            Qt(Indsep_t,i_t) = Qt(Indsep_t,i_t) - Uinf*vt(Indsep_t,1,i_t)/2 - Winf*vt(Indsep_t,2,i_t)/2;
            Qt(Indsep_t-1,i_t) = Qt(Indsep_t-1,i_t) - Uinf*vt(Indsep_t-1,1,i_t)/2 - Winf*vt(Indsep_t-1,2,i_t)/2;
        end
        if Indsep_b > 1
            Qt(Indsep_b,i_t) = Qt(Indsep_b,i_t) - Uinf*vt(Indsep_b,1,i_t)/2 - Winf*vt(Indsep_b,2,i_t)/2;
            Qt(Indsep_b-1,i_t) = Qt(Indsep_b-1,i_t) - Uinf*vt(Indsep_b-1,1,i_t)/2 - Winf*vt(Indsep_b-1,2,i_t)/2;
        end
    end
    % Calculating the pressure coefficient.
    Qinf = sqrt(Q0(:,i_t)'*Q0(:,i_t));
    Cp_s = 1 - Qt.^2/Qinf^2;
    
    %%    
    if i_t == 1
        Cp_us = 2/Qinf^2*mu(:,1)/delT;
    elseif i_t == 2
        Cp_us = 2/Qinf^2*(mu(:,1) - mu(:,2))/delT;
    else
        Cp_us = 2/Qinf^2*(3*mu(:,1) - 4*mu(:,2) + mu(:,3))/2/delT;
    end
    Cp = Cp_s + Cp_us;
    
%     % Calculating the pressure coefficient.
    P_s = -rho*Qp.^2/2;
    
    %%    
    if i_t == 1
        P_us = rho*mu(:,1)/delT;
    elseif i_t == 2
        P_us = rho*(mu(:,1) - mu(:,2))/delT;
    else
        P_us = rho*(3*mu(:,1) - 4*mu(:,2) + mu(:,3))/2/delT;
    end
    P = P_s + P_us;

    %% Boundary layer solver.  
    if ViscDrag
%         [dFi_pct,dL_pct,factor,xbc_sep,xtc_sep,Indsep_b,Indsep_t,Stagpt] = SkinFrictionSolver(xc,Qt,nu,rho,c,tmax,Cp,Qinf);   
        [dFi_pct,dL_pct,factor,xbc_sep,xtc_sep] = SkinFrictionSolver(xc,Qt,nu,rho,c,tmax,Cp,Qinf,0);  
        Indsep_b = 1;
        Indsep_t = 1;
        Stagpt = 75;
        dFi = dFi_pct.*kron(dL*b,ones(factor,1));
        dFshear = real(kron(eye(Npanels),ones(1,factor))*dFi);   
    end
    
    % Eliminating NaNs if present from drag calculation
    dF_Ind = [];
    dF_Ind = find(isnan(dFshear));
    dFshear(dF_Ind) = zeros(1,length(dF_Ind));

    % Setting separation points for first two timesteps.
    if i_t < 2
        Indsep_b = 1;
        Indsep_t = Npanels + 1;
    end
        
    % Restricting the separation point to the leading edge if it wraps
    % around the leading edge.
    
    if Indsep_t < Indsep_tmax
        Indsep_t = Indsep_tmax;
    end
    
    if Indsep_b > Indsep_bmax
        Indsep_b = Indsep_bmax;
    end

    
%% Calculating the forces
    % Calculating lift on the body and pitching moment about the leading
    % edge.
    delFp = -kron((Cp*1/2*rho*Qinf^2.*dL*b),[1 1]).*vn;
%     delFp = -kron((P.*dL*b),[1 1]).*vn;
    delFs = ViscDrag*kron(dFshear,[1 1]).*vt;
    delF = delFp + delFs;
    delP = sum(-delF.*Vc,2);  
    
    F = sum(delF,1);
    Fvisc = sum(delFs,1);
    D_visc = Fvisc(1)
    
    if BodDrag == 1
        if i_t == 1
            Fx = F(1) + 1/2*Cd_bod*rho*A_bod*Q0(1,i_t)^2;
        else
            Fx = F(1) + 1/2*Cd_bod*rho*A_bod*Q0(1,i_t-1)^2;
        end
%         if i_t == 1
%             Fx = F(1) + Cd_bod*rho*A_bod*abs(Q0(1,i_t))^(3/2)*(nu/BL)^(1/2);
%         else
%             Fx = F(1) + Cd_bod*rho*A_bod*abs(Q0(1,i_t-1))^(3/2)*(nu/BL)^(1/2);
%         end
        
        Fz = F(2);
    else
        Fx = F(1);
        Fz = F(2);
    end
    
    Pow = sum(delP,1);
    
    L = Fz*cos(alpha) - Fx*sin(alpha);
    T = -(Fz*sin(alpha) + Fx*cos(alpha));
      
    % M0 = -sum(delL*cos(alpha).*xp(2:end-1));

    % Calculating non-dimensional coefficients.
    % Cm = M0/(1/2*rho*Qinf^2*c^2)
    Gamma = mu(Npanels,1) - mu(1,1);

    Cf = norm(F)/(1/2*rho*Qinf^2*c*b);
    Cl = L/(1/2*rho*Qinf^2*c*b);
    Ct = T/(1/2*rho*Qinf^2*c*b);
    Cpow = Pow/(1/2*rho*Qinf^3*c*b);
    
    if i_t >= (Ncyc-1)*Nstep + 1
        Cl_avg = Cl_avg + Cl/Nstep;
        Ct_avg = Ct_avg + Ct/Nstep;
        Tnet_avg = Tnet_avg + T/Nstep;
        D_avg = D_avg + D_visc/Nstep;
        Pow_avg = Pow_avg + Pow/Nstep;
        Cpow_avg = Cpow_avg + Cpow/Nstep;
    end

    %% Decomposing steady and unsteady force components
    delFp_s = -kron((Cp_s*1/2*rho*Qinf^2.*dL*b),[1 1]).*vn;
    delFp_us = -kron((Cp_us*1/2*rho*Qinf^2.*dL*b),[1 1]).*vn;
    delF_s = delFp_s;
    delF_us = delFp_us;
    delP_s = sum(-delF_s.*Vc,2);
    delP_us = sum(-delF_us.*Vc,2);
    
    F_s = sum(delF_s,1);
    Pow_s = sum(delP_s,1);
    Fx_s = F_s(1);
    Fz_s = F_s(2);
    
    F_us = sum(delF_us,1);
    Pow_us = sum(delP_us,1);
    Fx_us = F_us(1);
    Fz_us = F_us(2);
    
    L_s = Fz_s*cos(alpha) - Fx_s*sin(alpha);
    T_s = -(Fz_s*sin(alpha) + Fx_s*cos(alpha));
    
    Cl_s = L_s/(1/2*rho*Qinf^2*c*b);
    Ct_s = T_s/(1/2*rho*Qinf^2*c*b);
    Cpow_s = Pow_s/(1/2*rho*Qinf^3*c*b);
    
    L_us = Fz_us*cos(alpha) - Fx_us*sin(alpha);
    T_us = -(Fz_us*sin(alpha) + Fx_us*cos(alpha));
    
    Cl_us = L_us/(1/2*rho*Qinf^2*c*b);
    Ct_us = T_us/(1/2*rho*Qinf^2*c*b);
    Cpow_us = Pow_us/(1/2*rho*Qinf^3*c*b);

    %% Wake rollup calculations
    if Rollup == 1 && i_t > 1
        
        if i_t > Nlump*Nstep + 1
            lump = 1;
        else
            lump = [];
        end

        [u,w] = WakeRollupVelocity2D([xw(2:wakeInd+1) lump*xl(1)],[zw(2:wakeInd+1) lump*zl(1)],mu(:,1)',muTE(1),[muW(1:wakeInd) lump*muLump(1)],muLEt,muLEb,sigma',xp',zp',vt',vn',xTE',zTE',vtTE',vnTE',[xw(1:wakeInd+1) lump*xl(1)],[zw(1:wakeInd+1) lump*zl(1)],[vtw(1:wakeInd,:); lump*vtlw(1) lump*vtlw(2)]',[vnw(1:wakeInd,:); lump*vnlw(1) lump*vnlw(2)]',xpt_LES,zpt_LES,vtLEt,vnLEt,xpb_LES,zpb_LES,vtLEb,vnLEb,i_t,Ncyc,Nstep,LES,ep,epBod,SC);
                
        u_2 = 0*u;
        w_2 = 0*w;
        
        if grd == 1
            [u_2,w_2] = WakeRollupVelocity2D([xw(2:wakeInd+1) lump*xl(1)],[zw(2:wakeInd+1) lump*zl(1)],mu(:,1)',muTE(1),[muW(1:wakeInd) lump*muLump(1)],muLEt,muLEb,sigma',xp',-zp',[vt(:,1) -vt(:,2)]',[vn(:,1) -vn(:,2)]',xTE',-zTE',[vtTE(1,1) -vtTE(1,2)]',[vnTE(1,1) -vnTE(1,2)]',[xw(1:wakeInd+1) lump*xl(1)],[-zw(1:wakeInd+1) -lump*zl(1)],[[vtw(1:wakeInd,1) -vtw(1:wakeInd,2)]; lump*vtlw(1) -lump*vtlw(2)]',[[vnw(1:wakeInd,1) -vnw(1:wakeInd,2)]; lump*vnlw(1) -lump*vnlw(2)]',xpt_LES,zpt_LES,vtLEt,vnLEt,xpb_LES,zpb_LES,vtLEb,vnLEb,i_t,Ncyc,Nstep,LES,ep,epBod,SC);
        end
        
        % Applying fencing scheme
        [xstar_w,zstar_w,fence] = fencing(c,Vp,[xw(2:wakeInd+1) lump*xl(1)],[zw(2:wakeInd+1) lump*zl(1)],u + u_2,w + w_2,xp,zp,vn,x_b(2),z_b(2),delT,1/100*epBod);
        xstar_w = xstar_w.*fence + (xstar_w + (u + u_2)*delT).*(~fence);
        zstar_w = zstar_w.*fence + (zstar_w + (w + w_2)*delT).*(~fence);
        
        % Wake rollup
        if i_t > Nlump*Nstep + 1
            xw(2:wakeInd+1) = xstar_w(1:end-1);
            zw(2:wakeInd+1) = zstar_w(1:end-1);

            vttemp = [(xw(1,2:end) - xw(1,1:end-1))' (zw(1,2:end) - zw(1,1:end-1))'];
            vtw = diag(1./sqrt(vttemp(:,1).^2 + vttemp(:,2).^2))*vttemp;
            vnw = [-vtw(:,2) vtw(:,1)]; 

            xl(1) = xstar_w(end);
            zl(1) = zstar_w(end);

            vttemp = [(xl(1) - xw(Nlump*Nstep+1))' (zl(1) - zw(Nlump*Nstep+1))'];
            vtlw(1,:) = diag(1./sqrt(vttemp(:,1).^2 + vttemp(:,2).^2))*vttemp;
            vnlw(1,:) = [-vtlw(1,2) vtlw(1,1)]; 
        else
            xw(2:wakeInd+1) = xstar_w;
            zw(2:wakeInd+1) = zstar_w;

            vttemp = [(xw(1,2:end) - xw(1,1:end-1))' (zw(1,2:end) - zw(1,1:end-1))'];
            vtw = diag(1./sqrt(vttemp(:,1).^2 + vttemp(:,2).^2))*vttemp;
            vnw = [-vtw(:,2) vtw(:,1)];
        end
            
        % LES rollup
        if LES == 1 && i_t > 2
            
            % Calculating wake rollup for top LE sheet.
            if Sep_t == 1
                
                % Calculating advection velocity on LEStop coordinates
                [uLES1,wLES1] = WakeRollupVelocity2D(xpt_LES(1,Ncyc*Nstep + 3 - i_t:end),zpt_LES(1,Ncyc*Nstep + 3 - i_t:end),mu(:,i_t)',muTE(i_t),muW(Ncyc*Nstep + 2 - i_t:end),muLEt,muLEb,sigma(:,i_t)',xp(:,i_t)',zp(:,i_t)',vt(:,:,i_t)',vn(:,:,i_t)',xTE(:,i_t)',zTE(:,i_t)',vtTE(1,:,i_t)',vnTE(1,:,i_t)',xw(Ncyc*Nstep + 2 - i_t:end),zw(Ncyc*Nstep + 2 - i_t:end),vtw(Ncyc*Nstep + 2 - i_t:end,:)',vnw(Ncyc*Nstep + 2 - i_t:end,:)',xpt_LES,zpt_LES,vtLEt,vnLEt,xpb_LES,zpb_LES,vtLEb,vnLEb,i_t,Ncyc,Nstep,LES,ep,epBod,SC);
                [uLES2,wLES2] = WakeRollupVelocity2D(xpt_LES(2,Ncyc*Nstep + 3 - i_t:end),zpt_LES(2,Ncyc*Nstep + 3 - i_t:end),mu(:,i_t)',muTE(i_t),muW(Ncyc*Nstep + 2 - i_t:end),muLEt,muLEb,sigma(:,i_t)',xp(:,i_t)',zp(:,i_t)',vt(:,:,i_t)',vn(:,:,i_t)',xTE(:,i_t)',zTE(:,i_t)',vtTE(1,:,i_t)',vnTE(1,:,i_t)',xw(Ncyc*Nstep + 2 - i_t:end),zw(Ncyc*Nstep + 2 - i_t:end),vtw(Ncyc*Nstep + 2 - i_t:end,:)',vnw(Ncyc*Nstep + 2 - i_t:end,:)',xpt_LES,zpt_LES,vtLEt,vnLEt,xpb_LES,zpb_LES,vtLEb,vnLEb,i_t,Ncyc,Nstep,LES,ep,epBod,SC);
            
                % Applying fencing scheme
                [xstar_t1,zstar_t1,fence1] = fencing(c,Vp(:,:,i_t),xpt_LES(1,Ncyc*Nstep + 3 - i_t:end),zpt_LES(1,Ncyc*Nstep + 3 - i_t:end),uLES1,wLES1,xp(:,i_t),zp(:,i_t),vn(:,:,i_t),x_b(1,i_t),z_b(1,i_t),delT,epB);           
                [xstar_t2,zstar_t2,fence2] = fencing(c,Vp(:,:,i_t),xpt_LES(2,Ncyc*Nstep + 3 - i_t:end),zpt_LES(2,Ncyc*Nstep + 3 - i_t:end),uLES2,wLES2,xp(:,i_t),zp(:,i_t),vn(:,:,i_t),x_b(1,i_t),z_b(1,i_t),delT,epB);           
                
                xstar_t = [xstar_t1; xstar_t2];
                zstar_t = [zstar_t1; zstar_t2];
                fence = [fence1; fence2];
                uLES = [uLES1; uLES2];
                wLES = [wLES1; wLES2];
                
                xstar_t = xstar_t.*fence + (xstar_t + uLES*delT).*(~fence);
                zstar_t = zstar_t.*fence + (zstar_t + wLES*delT).*(~fence);
            else
                % Calculating advection velocity on LEStop coordinates
                [uLES1,wLES1] = WakeRollupVelocity2D(xpt_LES(1,Ncyc*Nstep + 4 - i_t:end),zpt_LES(1,Ncyc*Nstep + 4 - i_t:end),mu(:,i_t)',muTE(i_t),muW(Ncyc*Nstep + 2 - i_t:end),muLEt,muLEb,sigma(:,i_t)',xp(:,i_t)',zp(:,i_t)',vt(:,:,i_t)',vn(:,:,i_t)',xTE(:,i_t)',zTE(:,i_t)',vtTE(1,:,i_t)',vnTE(1,:,i_t)',xw(Ncyc*Nstep + 2 - i_t:end),zw(Ncyc*Nstep + 2 - i_t:end),vtw(Ncyc*Nstep + 2 - i_t:end,:)',vnw(Ncyc*Nstep + 2 - i_t:end,:)',xpt_LES,zpt_LES,vtLEt,vnLEt,xpb_LES,zpb_LES,vtLEb,vnLEb,i_t,Ncyc,Nstep,LES,ep,epBod,SC);
                [uLES2,wLES2] = WakeRollupVelocity2D(xpt_LES(2,Ncyc*Nstep + 3 - i_t:end),zpt_LES(2,Ncyc*Nstep + 3 - i_t:end),mu(:,i_t)',muTE(i_t),muW(Ncyc*Nstep + 2 - i_t:end),muLEt,muLEb,sigma(:,i_t)',xp(:,i_t)',zp(:,i_t)',vt(:,:,i_t)',vn(:,:,i_t)',xTE(:,i_t)',zTE(:,i_t)',vtTE(1,:,i_t)',vnTE(1,:,i_t)',xw(Ncyc*Nstep + 2 - i_t:end),zw(Ncyc*Nstep + 2 - i_t:end),vtw(Ncyc*Nstep + 2 - i_t:end,:)',vnw(Ncyc*Nstep + 2 - i_t:end,:)',xpt_LES,zpt_LES,vtLEt,vnLEt,xpb_LES,zpb_LES,vtLEb,vnLEb,i_t,Ncyc,Nstep,LES,ep,epBod,SC);
                
                % Applying fencing scheme
                [xstar_t1,zstar_t1,fence1] = fencing(c,Vp(:,:,i_t),xpt_LES(1,Ncyc*Nstep + 4 - i_t:end),zpt_LES(1,Ncyc*Nstep + 4 - i_t:end),uLES1,wLES1,xp(:,i_t),zp(:,i_t),vn(:,:,i_t),x_b(1,i_t),z_b(1,i_t),delT,epB);           
                [xstar_t2,zstar_t2,fence2] = fencing(c,Vp(:,:,i_t),xpt_LES(2,Ncyc*Nstep + 3 - i_t:end),zpt_LES(2,Ncyc*Nstep + 3 - i_t:end),uLES2,wLES2,xp(:,i_t),zp(:,i_t),vn(:,:,i_t),x_b(1,i_t),z_b(1,i_t),delT,epB);           
                                
                xstar_t = [xpt_LES(1,Ncyc*Nstep + 3 - i_t) xstar_t1; xstar_t2];
                zstar_t = [zpt_LES(1,Ncyc*Nstep + 3 - i_t) zstar_t1; zstar_t2];
                fence = [1 fence1; fence2];
                uLES = [0 uLES1; uLES2];
                wLES = [0 wLES1; wLES2];
                
                xstar_t = xstar_t.*fence + (xstar_t + uLES*delT).*(~fence);
                zstar_t = zstar_t.*fence + (zstar_t + wLES*delT).*(~fence);
            end
           
            % Calculating wake rollup for bottom LE sheet.
            
            if Sep_b == 1
                % Calculating advection velocity on LESbot coordinates
                [uLES1,wLES1] = WakeRollupVelocity2D(xpb_LES(1,Ncyc*Nstep + 3 - i_t:end),zpb_LES(1,Ncyc*Nstep + 3 - i_t:end),mu(:,i_t)',muTE(i_t),muW(Ncyc*Nstep + 2 - i_t:end),muLEt,muLEb,sigma(:,i_t)',xp(:,i_t)',zp(:,i_t)',vt(:,:,i_t)',vn(:,:,i_t)',xTE(:,i_t)',zTE(:,i_t)',vtTE(1,:,i_t)',vnTE(1,:,i_t)',xw(Ncyc*Nstep + 2 - i_t:end),zw(Ncyc*Nstep + 2 - i_t:end),vtw(Ncyc*Nstep + 2 - i_t:end,:)',vnw(Ncyc*Nstep + 2 - i_t:end,:)',xpt_LES,zpt_LES,vtLEt,vnLEt,xpb_LES,zpb_LES,vtLEb,vnLEb,i_t,Ncyc,Nstep,LES,ep,epBod,SC);
                [uLES2,wLES2] = WakeRollupVelocity2D(xpb_LES(2,Ncyc*Nstep + 3 - i_t:end),zpb_LES(2,Ncyc*Nstep + 3 - i_t:end),mu(:,i_t)',muTE(i_t),muW(Ncyc*Nstep + 2 - i_t:end),muLEt,muLEb,sigma(:,i_t)',xp(:,i_t)',zp(:,i_t)',vt(:,:,i_t)',vn(:,:,i_t)',xTE(:,i_t)',zTE(:,i_t)',vtTE(1,:,i_t)',vnTE(1,:,i_t)',xw(Ncyc*Nstep + 2 - i_t:end),zw(Ncyc*Nstep + 2 - i_t:end),vtw(Ncyc*Nstep + 2 - i_t:end,:)',vnw(Ncyc*Nstep + 2 - i_t:end,:)',xpt_LES,zpt_LES,vtLEt,vnLEt,xpb_LES,zpb_LES,vtLEb,vnLEb,i_t,Ncyc,Nstep,LES,ep,epBod,SC);
            
                % Applying fencing scheme
                [xstar_b1,zstar_b1,fence1] = fencing(c,Vp(:,:,i_t),xpb_LES(1,Ncyc*Nstep + 3 - i_t:end),zpb_LES(1,Ncyc*Nstep + 3 - i_t:end),uLES1,wLES1,xp(:,i_t),zp(:,i_t),vn(:,:,i_t),x_b(1,i_t),z_b(1,i_t),delT,epB);           
                [xstar_b2,zstar_b2,fence2] = fencing(c,Vp(:,:,i_t),xpb_LES(2,Ncyc*Nstep + 3 - i_t:end),zpb_LES(2,Ncyc*Nstep + 3 - i_t:end),uLES2,wLES2,xp(:,i_t),zp(:,i_t),vn(:,:,i_t),x_b(1,i_t),z_b(1,i_t),delT,epB);           
                
                xstar_b = [xstar_b1; xstar_b2];
                zstar_b = [zstar_b1; zstar_b2];
                fence = [fence1; fence2];
                uLES = [uLES1; uLES2];
                wLES = [wLES1; wLES2];
                
                xstar_b = xstar_b.*fence + (xstar_b + uLES*delT).*(~fence);
                zstar_b = zstar_b.*fence + (zstar_b + wLES*delT).*(~fence);
            else
                % Calculating advection velocity on LESbot coordinates
                [uLES1,wLES1] = WakeRollupVelocity2D(xpb_LES(1,Ncyc*Nstep + 4 - i_t:end),zpb_LES(1,Ncyc*Nstep + 4 - i_t:end),mu(:,i_t)',muTE(i_t),muW(Ncyc*Nstep + 2 - i_t:end),muLEt,muLEb,sigma(:,i_t)',xp(:,i_t)',zp(:,i_t)',vt(:,:,i_t)',vn(:,:,i_t)',xTE(:,i_t)',zTE(:,i_t)',vtTE(1,:,i_t)',vnTE(1,:,i_t)',xw(Ncyc*Nstep + 2 - i_t:end),zw(Ncyc*Nstep + 2 - i_t:end),vtw(Ncyc*Nstep + 2 - i_t:end,:)',vnw(Ncyc*Nstep + 2 - i_t:end,:)',xpt_LES,zpt_LES,vtLEt,vnLEt,xpb_LES,zpb_LES,vtLEb,vnLEb,i_t,Ncyc,Nstep,LES,ep,epBod,SC);
                [uLES2,wLES2] = WakeRollupVelocity2D(xpb_LES(2,Ncyc*Nstep + 3 - i_t:end),zpb_LES(2,Ncyc*Nstep + 3 - i_t:end),mu(:,i_t)',muTE(i_t),muW(Ncyc*Nstep + 2 - i_t:end),muLEt,muLEb,sigma(:,i_t)',xp(:,i_t)',zp(:,i_t)',vt(:,:,i_t)',vn(:,:,i_t)',xTE(:,i_t)',zTE(:,i_t)',vtTE(1,:,i_t)',vnTE(1,:,i_t)',xw(Ncyc*Nstep + 2 - i_t:end),zw(Ncyc*Nstep + 2 - i_t:end),vtw(Ncyc*Nstep + 2 - i_t:end,:)',vnw(Ncyc*Nstep + 2 - i_t:end,:)',xpt_LES,zpt_LES,vtLEt,vnLEt,xpb_LES,zpb_LES,vtLEb,vnLEb,i_t,Ncyc,Nstep,LES,ep,epBod,SC);
            
                % Applying fencing scheme
                [xstar_b1,zstar_b1,fence1] = fencing(c,Vp(:,:,i_t),xpb_LES(1,Ncyc*Nstep + 4 - i_t:end),zpb_LES(1,Ncyc*Nstep + 4 - i_t:end),uLES1,wLES1,xp(:,i_t),zp(:,i_t),vn(:,:,i_t),x_b(1,i_t),z_b(1,i_t),delT,epB);           
                [xstar_b2,zstar_b2,fence2] = fencing(c,Vp(:,:,i_t),xpb_LES(2,Ncyc*Nstep + 3 - i_t:end),zpb_LES(2,Ncyc*Nstep + 3 - i_t:end),uLES2,wLES2,xp(:,i_t),zp(:,i_t),vn(:,:,i_t),x_b(1,i_t),z_b(1,i_t),delT,epB);           
                
                xstar_b = [xpb_LES(1,Ncyc*Nstep + 3 - i_t) xstar_b1; xstar_b2];
                zstar_b = [zpb_LES(1,Ncyc*Nstep + 3 - i_t) zstar_b1; zstar_b2];
                fence = [1 fence1; fence2];
                uLES = [0 uLES1; uLES2];
                wLES = [0 wLES1; wLES2];
                
                xstar_b = xstar_b.*fence + (xstar_b + uLES*delT).*(~fence);
                zstar_b = zstar_b.*fence + (zstar_b + wLES*delT).*(~fence);
            end
            
            % LES panel rollup.
            xpt_LES(:,Ncyc*Nstep + 3 - i_t:end) = xstar_t;
            zpt_LES(:,Ncyc*Nstep + 3 - i_t:end) = zstar_t;

            xpb_LES(:,Ncyc*Nstep + 3 - i_t:end) = xstar_b;
            zpb_LES(:,Ncyc*Nstep + 3 - i_t:end) = zstar_b;
            
            % Calculating normal and tangential vectors for LE panels.
            vttemp = [(xpt_LES(2,Ncyc*Nstep + 3 - i_t:end) - xpt_LES(1,Ncyc*Nstep + 3 - i_t:end))' (zpt_LES(2,Ncyc*Nstep + 3 - i_t:end) - zpt_LES(1,Ncyc*Nstep + 3 - i_t:end))'];
            vtLEt(Ncyc*Nstep + 3 - i_t:end,:) = diag(1./sqrt(vttemp(:,1).^2 + vttemp(:,2).^2))*vttemp;
            vnLEt(Ncyc*Nstep + 3 - i_t:end,:) = [-vtLEt(Ncyc*Nstep + 3 - i_t:end,2) vtLEt(Ncyc*Nstep + 3 - i_t:end,1)]; 

            vttemp = [(xpb_LES(2,Ncyc*Nstep + 3 - i_t:end) - xpb_LES(1,Ncyc*Nstep + 3 - i_t:end))' (zpb_LES(2,Ncyc*Nstep + 3 - i_t:end) - zpb_LES(1,Ncyc*Nstep + 3 - i_t:end))'];
            vtLEb(Ncyc*Nstep + 3 - i_t:end,:) = diag(1./sqrt(vttemp(:,1).^2 + vttemp(:,2).^2))*vttemp;
            vnLEb(Ncyc*Nstep + 3 - i_t:end,:) = [-vtLEb(Ncyc*Nstep + 3 - i_t:end,2) vtLEb(Ncyc*Nstep + 3 - i_t:end,1)]; 
        end
        
        
    end
    
    
    %% Storing the separation index.
    if LES == 1
        Sep_t_Store(1,i_t) = Indsep_t;
        Sep_b_Store(1,i_t) = Indsep_b;
    end

    if i_t > 1 && LES == 1
        % Determining when the separation location changes from the last 
        % timestep.
        if Sep_t_Store(i_t-1) ~= Sep_t_Store(i_t)
            Sep_t = 1;
        else
            Sep_t = 0;
        end

        if Sep_b_Store(i_t-1) ~= Sep_b_Store(i_t)
            Sep_b = 1;
        else
            Sep_b = 0;
        end
        
    end
    
    
    if LES == 1 
        % Normal vector at separation point
        if Indsep_t == Npanels + 1
            vnLESt = vn(1,:,i_t)/2 + vn(end,:,i_t)/2;
            vtLESt = vt(1,:,i_t)/2 + vt(end,:,i_t)/2;
        else
            vnLESt = vn(Indsep_t,:,i_t)/2 + vn(Indsep_t-1,:,i_t)/2;
            vtLESt = vt(Indsep_t,:,i_t)/2 + vt(Indsep_t-1,:,i_t)/2;
        end
        
        if Indsep_b == 1
            vnLESb = vn(1,:,i_t)/2 + vn(end,:,i_t)/2;
            vtLESb = vt(1,:,i_t)/2 + vt(end,:,i_t)/2;
        else
            vnLESb = vn(Indsep_b,:,i_t)/2 + vn(Indsep_b-1,:,i_t)/2;
            vtLESb = vt(Indsep_b,:,i_t)/2 + vt(Indsep_b-1,:,i_t)/2;
        end
        
        % Calculating a point away from the separation point in the
        % positive tangential vector direction.
        xplus_t = xp(Indsep_t,i_t) + 1/10*epB*(vtLESt(1,1));
        xplus_b = xp(Indsep_b,i_t) + 1/10*epB*(vtLESb(1,1));
        zplus_t = zp(Indsep_t,i_t) + 1/10*epB*(vtLESt(1,2));
        zplus_b = zp(Indsep_b,i_t) + 1/10*epB*(vtLESb(1,2));
        
        % Calculating a point away from the separation point in the
        % negative tangential vector direction.
        xminus_t = xp(Indsep_t,i_t) - 1/10*epB*(vtLESt(1,1));
        xminus_b = xp(Indsep_b,i_t) - 1/10*epB*(vtLESb(1,1));
        zminus_t = zp(Indsep_t,i_t) - 1/10*epB*(vtLESt(1,2));
        zminus_b = zp(Indsep_b,i_t) - 1/10*epB*(vtLESb(1,2));
        
        if t > 0
            % Calculating the velocity upstream and downstream of the
            % separation points.
            [Uplus_t,Wplus_t] = WakeRollupVelocity2D(xplus_t,zplus_t,mu(:,i_t)',muTE(i_t),muW(Ncyc*Nstep + 2 - i_t:end),muLEt,muLEb,sigma(:,i_t)',xp(:,i_t)',zp(:,i_t)',vt(:,:,i_t)',vn(:,:,i_t)',xTE(:,i_t)',zTE(:,i_t)',vtTE(1,:,i_t)',vnTE(1,:,i_t)',xw(Ncyc*Nstep + 2 - i_t:end),zw(Ncyc*Nstep + 2 - i_t:end),vtw(Ncyc*Nstep + 2 - i_t:end,:)',vnw(Ncyc*Nstep + 2 - i_t:end,:)',xpt_LES,zpt_LES,vtLEt,vnLEt,xpb_LES,zpb_LES,vtLEb,vnLEb,i_t,Ncyc,Nstep,LES,ep,epBod,SC);
            [Uplus_b,Wplus_b] = WakeRollupVelocity2D(xplus_b,zplus_b,mu(:,i_t)',muTE(i_t),muW(Ncyc*Nstep + 2 - i_t:end),muLEt,muLEb,sigma(:,i_t)',xp(:,i_t)',zp(:,i_t)',vt(:,:,i_t)',vn(:,:,i_t)',xTE(:,i_t)',zTE(:,i_t)',vtTE(1,:,i_t)',vnTE(1,:,i_t)',xw(Ncyc*Nstep + 2 - i_t:end),zw(Ncyc*Nstep + 2 - i_t:end),vtw(Ncyc*Nstep + 2 - i_t:end,:)',vnw(Ncyc*Nstep + 2 - i_t:end,:)',xpt_LES,zpt_LES,vtLEt,vnLEt,xpb_LES,zpb_LES,vtLEb,vnLEb,i_t,Ncyc,Nstep,LES,ep,epBod,SC);

            [Uminus_t,Wminus_t] = WakeRollupVelocity2D(xminus_t,zminus_t,mu(:,i_t)',muTE(i_t),muW(Ncyc*Nstep + 2 - i_t:end),muLEt,muLEb,sigma(:,i_t)',xp(:,i_t)',zp(:,i_t)',vt(:,:,i_t)',vn(:,:,i_t)',xTE(:,i_t)',zTE(:,i_t)',vtTE(1,:,i_t)',vnTE(1,:,i_t)',xw(Ncyc*Nstep + 2 - i_t:end),zw(Ncyc*Nstep + 2 - i_t:end),vtw(Ncyc*Nstep + 2 - i_t:end,:)',vnw(Ncyc*Nstep + 2 - i_t:end,:)',xpt_LES,zpt_LES,vtLEt,vnLEt,xpb_LES,zpb_LES,vtLEb,vnLEb,i_t,Ncyc,Nstep,LES,ep,epBod,SC);
            [Uminus_b,Wminus_b] = WakeRollupVelocity2D(xminus_b,zminus_b,mu(:,i_t)',muTE(i_t),muW(Ncyc*Nstep + 2 - i_t:end),muLEt,muLEb,sigma(:,i_t)',xp(:,i_t)',zp(:,i_t)',vt(:,:,i_t)',vn(:,:,i_t)',xTE(:,i_t)',zTE(:,i_t)',vtTE(1,:,i_t)',vnTE(1,:,i_t)',xw(Ncyc*Nstep + 2 - i_t:end),zw(Ncyc*Nstep + 2 - i_t:end),vtw(Ncyc*Nstep + 2 - i_t:end,:)',vnw(Ncyc*Nstep + 2 - i_t:end,:)',xpt_LES,zpt_LES,vtLEt,vnLEt,xpb_LES,zpb_LES,vtLEb,vnLEb,i_t,Ncyc,Nstep,LES,ep,epBod,SC);
        else
             % Calculating the velocity upstream and downstream of the
             % separation points.
            [Uplus_t,Wplus_t] = WakeRollupVelocity2D(xplus_t,zplus_t,mu(:,i_t)',muTE(i_t),muW(Ncyc*Nstep + 1 - i_t:end),muLEt,muLEb,sigma(:,i_t)',xp(:,i_t)',zp(:,i_t)',vt(:,:,i_t)',vn(:,:,i_t)',xTE(:,i_t)',zTE(:,i_t)',vtTE(1,:,i_t)',vnTE(1,:,i_t)',xw(Ncyc*Nstep + 1 - i_t:end),zw(Ncyc*Nstep + 1 - i_t:end),vtw(Ncyc*Nstep + 1 - i_t:end,:)',vnw(Ncyc*Nstep + 1 - i_t:end,:)',xpt_LES,zpt_LES,vtLEt,vnLEt,xpb_LES,zpb_LES,vtLEb,vnLEb,i_t,Ncyc,Nstep,LES,ep,epBod,SC);
            [Uplus_b,Wplus_b] = WakeRollupVelocity2D(xplus_b,zplus_b,mu(:,i_t)',muTE(i_t),muW(Ncyc*Nstep + 1 - i_t:end),muLEt,muLEb,sigma(:,i_t)',xp(:,i_t)',zp(:,i_t)',vt(:,:,i_t)',vn(:,:,i_t)',xTE(:,i_t)',zTE(:,i_t)',vtTE(1,:,i_t)',vnTE(1,:,i_t)',xw(Ncyc*Nstep + 1 - i_t:end),zw(Ncyc*Nstep + 1 - i_t:end),vtw(Ncyc*Nstep + 1 - i_t:end,:)',vnw(Ncyc*Nstep + 1 - i_t:end,:)',xpt_LES,zpt_LES,vtLEt,vnLEt,xpb_LES,zpb_LES,vtLEb,vnLEb,i_t,Ncyc,Nstep,LES,ep,epBod,SC);

            [Uminus_t,Wminus_t] = WakeRollupVelocity2D(xminus_t,zminus_t,mu(:,i_t)',muTE(i_t),muW(Ncyc*Nstep + 1 - i_t:end),muLEt,muLEb,sigma(:,i_t)',xp(:,i_t)',zp(:,i_t)',vt(:,:,i_t)',vn(:,:,i_t)',xTE(:,i_t)',zTE(:,i_t)',vtTE(1,:,i_t)',vnTE(1,:,i_t)',xw(Ncyc*Nstep + 1 - i_t:end),zw(Ncyc*Nstep + 1 - i_t:end),vtw(Ncyc*Nstep + 1 - i_t:end,:)',vnw(Ncyc*Nstep + 1 - i_t:end,:)',xpt_LES,zpt_LES,vtLEt,vnLEt,xpb_LES,zpb_LES,vtLEb,vnLEb,i_t,Ncyc,Nstep,LES,ep,epBod,SC);
            [Uminus_b,Wminus_b] = WakeRollupVelocity2D(xminus_b,zminus_b,mu(:,i_t)',muTE(i_t),muW(Ncyc*Nstep + 1 - i_t:end),muLEt,muLEb,sigma(:,i_t)',xp(:,i_t)',zp(:,i_t)',vt(:,:,i_t)',vn(:,:,i_t)',xTE(:,i_t)',zTE(:,i_t)',vtTE(1,:,i_t)',vnTE(1,:,i_t)',xw(Ncyc*Nstep + 1 - i_t:end),zw(Ncyc*Nstep + 1 - i_t:end),vtw(Ncyc*Nstep + 1 - i_t:end,:)',vnw(Ncyc*Nstep + 1 - i_t:end,:)',xpt_LES,zpt_LES,vtLEt,vnLEt,xpb_LES,zpb_LES,vtLEb,vnLEb,i_t,Ncyc,Nstep,LES,ep,epBod,SC);
        end
        
        % Averaging the upstream and downstream velocities for the finding
        % the shear velocity
        Vshear_t(:,i_t + 1) = [1/2*Uplus_t + 1/2*Uminus_t + Uinf; 1/2*Wplus_t + 1/2*Wminus_t + Winf];
        Vshear_b(:,i_t + 1) = [1/2*Uplus_b + 1/2*Uminus_b + Uinf; 1/2*Wplus_b + 1/2*Wminus_b + Winf];
    end
       
%% Plotting

   
if i_t > 1 && PlotTimeStepFig == 1 
   
%     pause(0.25)
    % Closing previous figure 
    if i_t > 3
        close(fighand(i_t-2))
    end
    
    FontSizeAx = 24;
    FontSizeLb = 32;
    afFigurePosition = [0 15 20 10];
    axespos = [0.15 0.2 0.84 0.8];
    ylabelpos = [-0.12 0.5];
    xlabelpos = [0.5 -0.25];


    
    % Plotting airfoil with LE vortex sheets and the TE vortex sheet.
    fighand(i_t) = figure;
    set(fighand(i_t), 'Units', 'centimeters','PaperPositionMode', 'auto','Position', afFigurePosition);
    set(fighand(i_t),'DefaultAxesFontSize',FontSizeAx,'DefaultAxesFontName','TimesNewRoman','DefaultAxesGridLineStyle','-.','DefaultAxesLineWidth',2,'DefaultAxesFontWeight','Normal')

    hold on
    axis equal
    plot(xp,zp,'-k','linewidth',2)
    plot(xTE,zTE,'.-b','linewidth',2)
%     val = 1/10;
    
    if grd == 1
        plot(xp,-zp,'-k','linewidth',2)
        plot(xTE,-zTE,'.-b','linewidth',2)
        plot([min(xp)-10*c; min(xp)+10*c],[0 0],'-k','linewidth',2)
    end
  
    if LES == 1
        plot(xpt_LES(1,Ncyc*Nstep + 2 - i_t),zpt_LES(1,Ncyc*Nstep + 2 - i_t),'og','linewidth',2)
        plot(xpb_LES(1,Ncyc*Nstep + 2 - i_t),zpb_LES(1,Ncyc*Nstep + 2 - i_t),'or','linewidth',2)

        for j = i_t:-1:2
            if muLEt(Ncyc*Nstep + 2 - j) == 0
            else
                plot(xpt_LES(:,Ncyc*Nstep + 2 - j),zpt_LES(:,Ncyc*Nstep + 2 - j),'.-g','linewidth',2,'markersize',14)
%                 quiver(xpt_LES(1,Ncyc*Nstep + 1 - j),zpt_LES(1,Ncyc*Nstep + 1 - j),val*vtLEt(Ncyc*Nstep + 1 - j,1),val*vtLEt(Ncyc*Nstep + 1 - j,2),'b')
%                 quiver(xpt_LES(1,Ncyc*Nstep + 1 - j),zpt_LES(1,Ncyc*Nstep + 1 - j),val*vnLEt(Ncyc*Nstep + 1 - j,1),val*vnLEt(Ncyc*Nstep + 1 - j,2),'k')

            end
        end

        for j = i_t:-1:2
            if muLEb(Ncyc*Nstep + 2 - j) == 0  
            else
                plot(xpb_LES(:,Ncyc*Nstep + 2 - j),zpb_LES(:,Ncyc*Nstep + 2 - j),'.-r','linewidth',2,'markersize',14)
%                 quiver(xpb_LES(1,Ncyc*Nstep + 1 - j),zpb_LES(1,Ncyc*Nstep + 1 - j),val*vtLEb(Ncyc*Nstep + 1 - j,1),val*vtLEb(Ncyc*Nstep + 1 - j,2),'b')
%                 quiver(xpb_LES(1,Ncyc*Nstep + 1 - j),zpb_LES(1,Ncyc*Nstep + 1 - j),val*vnLEb(Ncyc*Nstep + 1 - j,1),val*vnLEb(Ncyc*Nstep + 1 - j,2),'k')

            end
        end
    end
%         plot(xp(Stagpt,i_t),zp(Stagpt,i_t),'xk','linewidth',2)

    bluevec = [linspace(0,1,100)' linspace(0,1,100)' ones(100,1)];
    redvec = [ones(100,1) linspace(1,0,100)' linspace(1,0,100)'];
    vec = [bluevec; redvec(2:end,:)];

    WakeCirc = GammaW'/max(abs(GammaW))/(1/1);
    WakeCirc(WakeCirc > 1) = 1;  
    WakeCirc(WakeCirc < -1) = -1;
    WakeCirc = WakeCirc + 1;
    WakeCirc = round(WakeCirc*100) + 1;
    WakeCirc(WakeCirc > 199) = 199;

    WakeCirc_w = -GammaW'/max(abs(GammaW))/(1/1);
    WakeCirc_w(WakeCirc_w > 1) = 1;  
    WakeCirc_w(WakeCirc_w < -1) = -1;
    WakeCirc_w = WakeCirc_w + 1;
    WakeCirc_w = round(WakeCirc_w*100) + 1;
    WakeCirc_w(WakeCirc_w > 199) = 199;
    
    WakeCirc_l = -muLump(1)'/max(abs(GammaW))/(1/1);
    WakeCirc_l(WakeCirc_l > 1) = 1;  
    WakeCirc_l(WakeCirc_l < -1) = -1;
    WakeCirc_l = WakeCirc_l + 1;
    WakeCirc_l = round(WakeCirc_l*100) + 1;
    WakeCirc_l(WakeCirc_l > 199) = 199;
    
    WakeCirc_lw = muLump(1)'/max(abs(GammaW))/(1/1);
    WakeCirc_lw(WakeCirc_lw > 1) = 1;  
    WakeCirc_lw(WakeCirc_lw < -1) = -1;
    WakeCirc_lw = WakeCirc_lw + 1;
    WakeCirc_lw = round(WakeCirc_lw*100) + 1;
    WakeCirc_lw(WakeCirc_lw > 199) = 199;
    
    for i_w = 1:wakeInd+1
        plot(xw(i_w),zw(i_w),'.-','color',vec(WakeCirc(i_w),:),'linewidth',2,'markersize',14)
    end
    if grd == 1 && i_t > 1
        for i_w = 1:wakeInd+1
            plot(xw(i_w),-zw(i_w),'.','color',vec(WakeCirc_w(i_w),:),'linewidth',2,'markersize',14)
        end
    end
    if i_t > Nlump*Nstep + 1
        plot(xl(1),zl(1),'.','color',vec(WakeCirc_l,:),'linewidth',2,'markersize',20)
        if grd == 1
            plot(xl(1),-zl(1),'.','color',vec(WakeCirc_lw,:),'linewidth',2,'markersize',20)
        end
    end
   
%     plot((xp(Npanels/2 + 1) + upfence*c)*[1 1],2*c*[-1 1],'--k','linewidth',1)

    hold off
    axis([-c/2 + x_b(2) 6*c + x_b(2) z_b(2) - 1.5*c z_b(2) + 1.5*c])

    set(gca, 'FontName', 'TimesNewRoman', 'FontSize', FontSizeAx)
    set(gca, 'Units', 'normalized', 'Position', axespos);
    % legend('d/c = \infty','d/c = 1/2','d/c = 1/4','Location','NorthWest')
    xlabel('$$x$$ (m)','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 0,'Units', 'Normalize','Position',xlabelpos,'Interpreter', 'LaTeX');
    ylabel('$$z$$ (m)','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 90,'Units', 'Normalize','Position',ylabelpos,'Interpreter', 'LaTeX');

%     print('-depsc','-r300',['Free_eps/2D_Pitch_Free_',num2str(i_t)]);
%     print('-dpng','-r300',['Grd_png_free_lump4/2D_Pitch_Grd_',num2str(i_t)]);
%     pause(0.01);
 end
    
%   
%     
% %     if i_t > 2
% %         close(fighand)
% %     end
% %     
% %     % Plotting discretized airfoil with normals and lift vectors.
% %     fighand = figure;
% %     hold on
% %     axis equal
% %     plot(xp(:,i_t),zp(:,i_t),'-k')
% %     plot(xp(Stagpt,i_t),zp(Stagpt,i_t),'ok','linewidth',2)
% %     
% % %     plot(xpt_LES(1,Ncyc*Nstep + 2 - i_t),zpt_LES(1,Ncyc*Nstep + 2 - i_t),'og','linewidth',2)
% % %     plot(xpt_LES(1,Ncyc*Nstep + 2 - i_t:end-1),zpt_LES(1,Ncyc*Nstep + 2 - i_t:end-1),'.-g','linewidth',1)
% % %     plot(xpb_LES(1,Ncyc*Nstep + 2 - i_t),zpb_LES(1,Ncyc*Nstep + 2 - i_t),'or','linewidth',2)
% % %     plot(xpb_LES(1,Ncyc*Nstep + 2 - i_t:end-1),zpb_LES(1,Ncyc*Nstep + 2 - i_t:end-1),'.-r','linewidth',1)
% %     
% %     plot(xTE(:,i_t),zTE(:,i_t),'.-b')
% % %     quiver(xc(:,i_t),zc(:,i_t),Vc(:,1,i_t),Vc(:,2,i_t),'r')
% % 
% %     if t > 0
% %         plot(xw(Ncyc*Nstep + 2 - i_t:end),zw(Ncyc*Nstep + 2 - i_t:end),'.-b')
% % 
% % %         quiver(xw(Ncyc*Nstep + 2 - i_t:end-1)',zw(Ncyc*Nstep + 2 - i_t:end-1)',vnw(Ncyc*Nstep + 2 - i_t:end,1),vnw(Ncyc*Nstep + 2 - i_t:end,2))
% % %         quiver(xw(Ncyc*Nstep + 2 - i_t:end-1)',zw(Ncyc*Nstep + 2 - i_t:end-1)',vtw(Ncyc*Nstep + 2 - i_t:end,1),vtw(Ncyc*Nstep + 2 - i_t:end,2))
% %     end
% %     
% %     xlabel('$$x$$','interpreter','latex','fontsize',12,'fontname','TimesNewRoman')
% %     ylabel('$$z$$','interpreter','latex','fontsize',12,'fontname','TimesNewRoman')
% % % 
% %     end
    if Flowfield == 1 && t > 0
        % Calculating the flowfield around the wing.
        nx = 51;
        nz = 51;

        U = Uinf*ones(nz,nx);
        W = Winf*ones(nz,nx);
       
        % Flow field for ground effect calculations
        xf = linspace(-c/4 + x_b(2),2*c + x_b(2),nx)';
        zf = linspace(-c/4,c/4,nz)';

        [Xf,Zf] = meshgrid(xf,zf);

        X = zeros(1,nx*nz);
        Z = zeros(1,nx*nz);

        for i = 1:nz
            vec = (i-1)*nx + 1:i*nx;

            X(1,vec) = Xf(i,:);
            Z(1,vec) = Zf(i,:);
        end

        % Calculating flow field

        % Body contribution
        [u_st,w_st,u_bdt,w_bdt] = DubSorV(mu(:,1)',sigma',X,Z,xp',zp',vt',vn',epSC,SC);

        u_st_2 = 0*u_st;
        w_st_2 = 0*w_st;
        u_bdt_2 = 0*u_bdt;
        w_bdt_2 = 0*w_bdt;
        
        if grd == 1
            [u_st_2,w_st_2,u_bdt_2,w_bdt_2] = DubSorV(mu(:,1)',sigma',X,Z,xp_2',zp_2',vt_2',vn_2',epSC,SC);
        end
        
        u_s = zeros(nz,nx);
        w_s = zeros(nz,nx);
        u_bd = zeros(nz,nx);
        w_bd = zeros(nz,nx);
        
        for i = 1:nz
                vec = (i-1)*nx + 1:i*nx;

                u_s(i,:) = u_st(vec) +  u_st_2(vec);
                w_s(i,:) = w_st(vec) + w_st_2(vec);
                u_bd(i,:) = u_bdt(vec) + u_bdt_2(vec);
                w_bd(i,:) = w_bdt(vec) + w_bdt_2(vec);
        end


        % TE contribution  
        [~,~,u_TEdt,w_TEdt] = DubSorV(muTE(1),0*muTE(1),X,Z,xTE',zTE',vtTE',vnTE',epSC,SC);

        u_TEdt_2 = 0*u_TEdt;
        w_TEdt_2 = 0*w_TEdt;
        
        if grd == 1
            [~,~,u_TEdt_2,w_TEdt_2] = DubSorV(muTE(1),0*muTE(1),X,Z,xTE_2',zTE_2',vtTE_2',vnTE_2',epSC,SC);
        end
                    
        u_TEd = zeros(nz,nx);
        w_TEd = zeros(nz,nx);
        for i = 1:nz
                vec = (i-1)*nx + 1:i*nx;

                u_TEd(i,:) = u_TEdt(vec) + u_TEdt_2(vec);
                w_TEd(i,:) = w_TEdt(vec) + w_TEdt_2(vec);
        end

        % Wake contribution
        [~,~,u_wdt,w_wdt] = DubSorV(muW(1:wakeInd),0*muW(1:wakeInd),X,Z,xw(1:wakeInd+1),zw(1:wakeInd+1),vtw(1:wakeInd,:)',vnw(1:wakeInd,:)',epSC,SC);

        u_wdt_2 = 0*u_wdt;
        w_wdt_2 = 0*w_wdt;
        
        if grd == 1
            [~,~,u_wdt_2,w_wdt_2] = DubSorV(muW(1:wakeInd),0*muW(1:wakeInd),X,Z,xw_2(1:wakeInd+1),zw_2(1:wakeInd+1),vtw_2(1:wakeInd,:)',vnw_2(1:wakeInd,:)',epSC,SC);
        end
        
        u_wd = zeros(nz,nx);
        w_wd = zeros(nz,nx);
        for i = 1:nz
                vec = (i-1)*nx + 1:i*nx;

                u_wd(i,:) = u_wdt(vec) + u_wdt_2(vec);
                w_wd(i,:) = w_wdt(vec) + w_wdt_2(vec);
        end

        if LES == 1
            vec = i_t:-1:2;
            u_LEtdt = zeros(1,length(X),length(vec));
            u_LEbdt = zeros(1,length(X),length(vec));
            w_LEtdt = zeros(1,length(X),length(vec));
            w_LEbdt = zeros(1,length(X),length(vec));

            for j = vec
                % LE sheet top contribution
                [~,~,u_LEtdt(1,:,j-1),w_LEtdt(1,:,j-1)] = DubSorV(muLEt(Ncyc*Nstep + 2 - j),0*muLEt(Ncyc*Nstep + 2 - j),X,Z,xpt_LES(:,Ncyc*Nstep + 2 - j)',zpt_LES(:,Ncyc*Nstep + 2 - j),vtLEt(Ncyc*Nstep + 2 - j,:)',vnLEt(Ncyc*Nstep + 2 - j,:)',epSC,SC);

                % LE sheet bot contribution
                [~,~,u_LEbdt(1,:,j-1),w_LEbdt(1,:,j-1)] = DubSorV(muLEb(Ncyc*Nstep + 2 - j),0*muLEb(Ncyc*Nstep + 2 - j),X,Z,xpb_LES(:,Ncyc*Nstep + 2 - j)',zpb_LES(:,Ncyc*Nstep + 2 - j),vtLEb(Ncyc*Nstep + 2 - j,:)',vnLEb(Ncyc*Nstep + 2 - j,:)',epSC,SC);
            end

            u_LEtdt = sum(u_LEtdt,3);
            u_LEbdt = sum(u_LEbdt,3);
            w_LEtdt = sum(w_LEtdt,3);
            w_LEbdt = sum(w_LEbdt,3);

            u_LEtd = zeros(nz,nx);
            w_LEtd = zeros(nz,nx);
            u_LEbd = zeros(nz,nx);
            w_LEbd = zeros(nz,nx);
            for i = 1:nz
                    vec = (i-1)*nx + 1:i*nx;

                    u_LEtd(i,:) = u_LEtdt(vec);
                    w_LEtd(i,:) = w_LEtdt(vec);
                    u_LEbd(i,:) = u_LEbdt(vec);
                    w_LEbd(i,:) = w_LEbdt(vec);
            end
        end

        if LES == 1
            u_d = u_bd + u_TEd + u_wd + u_LEtd + u_LEbd;
            w_d = w_bd + w_TEd + w_wd + w_LEtd + w_LEbd;
        else
            u_d = u_bd + u_TEd + u_wd;
            w_d = w_bd + w_TEd + w_wd;
        end

        u_p = u_d + u_s;
        w_p = w_d + w_s;

        Ut = U + u_p;
        Wt = W + w_p;
        
        omega_y = (u_p(3:end,2:end-1) - u_p(1:end-2,2:end-1))./(Zf(3:end,2:end-1) - Zf(1:end-2,2:end-1)) - (w_p(2:end-1,3:end) - w_p(2:end-1,1:end-2))./(Xf(2:end-1,3:end) - Xf(2:end-1,1:end-2));
        Xstar = (Xf(2:end-1,3:end) + Xf(2:end-1,1:end-2))/2;
        Zstar = (Zf(3:end,2:end-1) + Zf(1:end-2,2:end-1))/2;

        if i_t > 2
            close(flowfig)
        end
        
        
         % Plotting airfoil with LE vortex sheets and the TE vortex sheet.
        flowfig = figure;
        FontSizeAx = 24;
        afFigurePosition = [15 7 25 15];
        
    
        set(flowfig, 'Units', 'centimeters','PaperPositionMode', 'auto','Position', afFigurePosition);
        set(flowfig,'DefaultAxesFontSize',FontSizeAx,'DefaultAxesFontName','TimesNewRoman','DefaultAxesGridLineStyle','-.','DefaultAxesLineWidth',2,'DefaultAxesFontWeight','Normal')

        hold on
        axis equal
        
        if grd == 1
           plot([min(xp)-10*c; min(xp)+10*c],[0 0],'-k','linewidth',1)
        end
        
        bluevec = [linspace(0,1,100)' linspace(0,1,100)' ones(100,1)];
        redvec = [ones(100,1) linspace(1,0,100)' linspace(1,0,100)'];
        vec = [bluevec; redvec(2:end,:)];
        colormap(vec)
        
        pcolor(Xstar,Zstar,-omega_y)
        quiver(Xf,Zf,u_p,w_p,'k')
        plot(xp,zp,'-k','linewidth',2)
        
%         if grd == 1
%             plot([min(xp(:,i_t))-10*c; min(xp(:,i_t))+10*c],[0 0],'-k','linewidth',4)
%         end
        
        shading interp
        
        if LES == 1
            plot(xplus_t,zplus_t,'xk')
            plot(xminus_t,zminus_t,'xk')
            plot(xplus_b,zplus_b,'xk')
            plot(xminus_b,zminus_b,'xk')
            quiver(xp(Indsep_t,i_t),zp(Indsep_t,i_t),Vshear_t(1,i_t+1)*delT,Vshear_t(2,i_t+1)*delT,'r')
            quiver(xp(Indsep_b,i_t),zp(Indsep_b,i_t),Vshear_b(1,i_t+1)*delT,Vshear_b(2,i_t+1)*delT,'r')
        end
%         plot(xTE(:,i_t),zTE(:,i_t),'.-b','linewidth',2)
        
        val = 1/10;

        if LES == 1 && t > 0
            plot(xpt_LES(1,Ncyc*Nstep + 2 - i_t),zpt_LES(1,Ncyc*Nstep + 2 - i_t),'og','linewidth',2)
            plot(xpb_LES(1,Ncyc*Nstep + 2 - i_t),zpb_LES(1,Ncyc*Nstep + 2 - i_t),'or','linewidth',2)

            for j = i_t:-1:2
                if muLEt(Ncyc*Nstep + 2 - j) == 0
                else
                    plot(xpt_LES(:,Ncyc*Nstep + 2 - j),zpt_LES(:,Ncyc*Nstep + 2 - j),'.-g','linewidth',1)
    %                 quiver(xpt_LES(1,Ncyc*Nstep + 1 - j),zpt_LES(1,Ncyc*Nstep + 1 - j),val*vtLEt(Ncyc*Nstep + 1 - j,1),val*vtLEt(Ncyc*Nstep + 1 - j,2),'b')
    %                 quiver(xpt_LES(1,Ncyc*Nstep + 1 - j),zpt_LES(1,Ncyc*Nstep + 1 - j),val*vnLEt(Ncyc*Nstep + 1 - j,1),val*vnLEt(Ncyc*Nstep + 1 - j,2),'k')

                end
            end

            for j = i_t:-1:2
                if muLEb(Ncyc*Nstep + 2 - j) == 0   
                else
                    plot(xpb_LES(:,Ncyc*Nstep + 2 - j),zpb_LES(:,Ncyc*Nstep + 2 - j),'.-r','linewidth',1)
    %                 quiver(xpb_LES(1,Ncyc*Nstep + 1 - j),zpb_LES(1,Ncyc*Nstep + 1 - j),val*vtLEb(Ncyc*Nstep + 1 - j,1),val*vtLEb(Ncyc*Nstep + 1 - j,2),'b')
    %                 quiver(xpb_LES(1,Ncyc*Nstep + 1 - j),zpb_LES(1,Ncyc*Nstep + 1 - j),val*vnLEb(Ncyc*Nstep + 1 - j,1),val*vnLEb(Ncyc*Nstep + 1 - j,2),'k')

                end
            end

        end
        
        if t > 0
%             plot(xp(Stagpt,i_t),zp(Stagpt,i_t),'xk','linewidth',2)
%             plot(xw(Ncyc*Nstep + 2 - i_t:end),zw(Ncyc*Nstep + 2 - i_t:end),'.-b')
        end
            
        

%         colorbar
        caxis(6*[-1 1])
        colorbar
        axis([-c/4 + x_b(2) 2*c + x_b(2) -c c])
%         axis([-c/2 + x_b(i_t) 3*c + x_b(i_t) -3*c/4 + z_b(i_t) 3*c/4 + z_b(i_t)])
%         xlabel('$$x$$','interpreter','latex','fontsize',12,'fontname','TimesNewRoman')
%         ylabel('$$z$$','interpreter','latex','fontsize',12,'fontname','TimesNewRoman')
        
%         VortMov(i_t-1) = getframe;
%         print('-dpng','-r300',['Grd_St_25_dc_16_Ac_17_',num2str(i_t)]);

%         print('-dpng','-r300',['2D_NearBody_',num2str(i_t)]);
    end

    
    
    
    %% Free-swimming calculations
    
    F_t(i_t) = Fx;   
    if free == 1 && i_t > 1       
        a_b = Fx/M;       
        Q0(1,i_t+1) = a_b*delT + Q0(1,i_t);
        
        x_b(1) = (Q0(1,i_t + 1) + Q0(1,i_t))*delT/2 + x_b(2); 

        
        if PlotVelFig == 1 && (rem(i_t,5) == 0)
            if i_t > 10
                close(velh(i_t-10))
            end

            FontSizeAx = 24;
            FontSizeLb = 32;
            afFigurePosition = [0 0 20 10];
            axespos = [0.15 0.15 0.84 0.75];
            ylabelpos = [-0.12 0.5];
            xlabelpos = [0.5 -0.25];

            velh(i_t) = figure;
            set(gcf, 'Units', 'centimeters','PaperPositionMode', 'auto','Position', afFigurePosition);
            set(gcf,'DefaultAxesFontSize',FontSizeAx,'DefaultAxesFontName','TimesNewRoman','DefaultAxesGridLineStyle','-.','DefaultAxesLineWidth',2,'DefaultAxesFontWeight','Normal')

            hold on
                plot([0 Ncyc/f],U_swim*[1 1],'-k','linewidth',1)
                plot(N_a/f*[1 1],[0 1.5*max(-Q0(1,1:i_t))],'--k','linewidth',1)
                plot(((1:i_t) - 1)*delT,-Q0(1,1:i_t)','-b','linewidth',2)
            hold off

            set(gca, 'FontName', 'TimesNewRoman', 'FontSize', FontSizeAx)
            set(gca, 'Units', 'normalized', 'Position', axespos);
            xlabel('$$t/T$$, sec','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 0,'Units', 'Normalize','Position',xlabelpos,'Interpreter', 'LaTeX');
            ylabel('$$U(t)$$','Interpreter', 'LaTeX','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 90,'Units', 'Normalize','Position',ylabelpos);
            axis([0 Ncyc/f 0 1.5*max(-Q0(1,1:i_t))])
        end
%         subplot(2,1,2)
%         plot((1:Ncyc*Nstep+1)*delT,D_visc,'-r','linewidth',2)
%         xlabel('$$t$$, sec','interpreter','latex','fontname','TimesNewRoman','fontsize',12)
%         ylabel('Viscous Drag','interpreter','latex','fontname','TimesNewRoman','fontsize',12)
%         axis([0 Ncyc/f 0 1.5*max(D_visc)])
    else
        x_b(1) = Q0(1,i_t)*delT + x_b(2);
    end
    z_b(1) = Q0(2,i_t)*delT + z_b(2);
    
        
%     close(wbar);
    delTime(i_t) = toc;
    EstT = delTime(i_t)*(Ncyc*Nstep+1 - i_t);
    EstThour = floor(EstT/3600);
    EstTmin = floor((EstT - EstThour*3600)/60);
    EstTsec = round(EstT - EstThour*3600 - EstTmin*60);
%     sprintf(['Calculating...Time Remaining: ',num2str(EstThour),' hrs ',num2str(EstTmin),' mins ',num2str(EstTsec),' secs']);

%     wbar = waitbar(i_t/(Ncyc*Nstep+1),['Calculating...Time Remaining: ',num2str(EstThour),' hrs ',num2str(EstTmin),' mins ',num2str(EstTsec),' secs'],'Position',[(scrsz(3) - scrsz(4)/3) 0 scrsz(4)/3 1/16*scrsz(4)]);
    
    if SaveData == 1           
        Data = [wakeInd,Q0(1,i_t),Q0(2,i_t),x_b(2),z_b(2),a_b,Fx,Fz,Pow,L,T,D_visc,Gamma,Cf,Cl,Ct,Cpow,...
            Fx_s,Fx_us,Fz_s,Fz_us,Pow_s,Pow_us,L_s,L_us,T_s,T_us,Cl_s,Cl_us,...
            Ct_s,Ct_us,Cpow_s,Cpow_us];
        fprintf(fid_Data,'%16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f\n',Data');

        PanelProp = [xp_0(1:end-1)'; zp_0(1:end-1)';xp(1:end-1)';zp(1:end-1)';...
            Vc(:,1)';Vc(:,2)';xc';zc';vt(:,1)';vt(:,2)';vn(:,1)';vn(:,2)';...
            dL';Xc';Zc';sigma';mu(:,1)';Qp';Qt';Cp_s';Cp_us';Cp';dFshear']';
        fprintf(fid_PanelProp,'%16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f\n',PanelProp');

        WakeProp = [xw xl(1); zw zl(1);vtTE' vtw' vtlw(1,:)';...
            vnTE' vnw' vnlw(1,:)'; muTE(1) muW muLump(1); 
            GammaW -muLump(1)]';
        fprintf(fid_WakeProp,'%16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f\n',WakeProp');
    end
    
    sprintf([num2str(round((i_t/(Ncyc*Nstep+1))*100)),'%% complete'],'FontName','TimesNewRoman')
end

% matlabpool close force local
% close(wbar);
if SaveData == 1
    fclose(fid_Data);
    fclose(fid_PanelProp);
    fclose(fid_WakeProp);
end

% End timer.
Time = sum(delTime);
Thour = floor(Time/3600);
Tmin = floor((Time - Thour*3600)/60);
Tsec = round(Time - Thour*3600 - Tmin*60);
sprintf(['Calculation Time: ',num2str(Thour),' hrs ',num2str(Tmin),' mins ',num2str(Tsec),' secs'],'FontName','TimesNewRoman')


% Time averaged forces
Ucruise = mean(-Q0(1,(Ncyc-1)*Nstep + 1:end))
Tnet_avg
Ct_avg
D_avg
Cl_avg
Ct_avg
Cpow_avg
Pow_avg
% T_avg = Ct_avg*(1/2*rho*Qinf^2*c*b);
np = Tnet_avg*Ucruise/Pow_avg
Ec = Ucruise/Pow_avg                        %km/kJ
V_t = -Q0(1,:)';

% profile viewer

