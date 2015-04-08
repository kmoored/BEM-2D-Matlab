% Written by Keith Moored, 2/22/12

% clear
% clc

function [Ct_avg,Cl_avg,Cpow_avg,np,Cl_end,Ct_end] = PanelMethod2DGrdSweep(d_c,St,Npanels,Nstep,val)
% d_c = 0.5
% St = 0.45
% Npanels = 210;
% Nstep = 80
% val = 1e-3;

%% Parameters: 

% Fluid parameters
rho = 1000;                 % kg/m^3
nu = 1.004e-6;              % m^2/s
% Re = 2.4*10^5;
Qinf = 0.06;

% Geometry parameters
c = 0.12;                   % meters
b = 0.28;                   % meters
tmax = 0.05;                % Percent of the chord
% d_c = 0.5;                    % Distance from the ground in chords
D = d_c*c;                  % Distance from the ground in m

% Paneling parameters
Ncyc = 3;
% Nstep = 100;
% Npanels = 100;
TEfac = 0.25;
LEfac = 1;
kappa = 8;
CircRed = 1;

% 2P Kin: f = 1-2, 2S Kin: f = 2.5, h_c = 1 c = 1/10, phi = -pi/2,
% alpha_max = 30.

% Kinematic parameters
alpha = 0*(pi/180);        % Enter in degrees, converts to radians.
h_c = 0;
% f = 0;
% St = 0.25;
% f_red = 2.5;
alpha_max = asin(0.125);     % Max pitch angle
% alpha_max = 7.5*(pi/180);     % Max pitch angle
% alpha_max = 0*(10 - 0.2316)*(pi/180);     % Max pitch angle
phi = -pi/2;                 % Phase angle between heaving and pitching
phase = pi/2;
M = 100;

% Settings
LES = 0;                                        % 0 - off, 1 - on. 
SC = 1;                                         % 0 - off, 1 - on. 
Rollup = 1;                                     % 0 - off, 1 - on.
Flowfield = 0;                                  % 0 - off, 1 - on.
wakeinf = 1;                                    % 0 - off, 1 - on.
bod = 0;                                        % 0 - off, 1 - on.
grd = 0;                                        % 0 - off, 1 - on.
free = 0;                                       % 0 - off, 1 - on.
SurfaceColor = 0;                               % 0 - shear stress, 1 - pressure
ViscDrag = 1;                                   % 0 - off, 1 - on.

% Cuttoff radii
% val = 1e-3;
ep = val*c;
epSC = ep;
epBod = ep;
epB = val*c

%% Calculations

A_TE = c*sin(alpha_max)
A_c = 2*A_TE/c
% h_c = 2*pi*St/f_red                     % Heave to chord ratio
% St = 
% A_TE = h_c*c;
% Qinf = Re*nu/c                       % m/s
Re = c*Qinf/nu
% f = 1;
f = St*Qinf/2/A_TE
% f = Qinf*f_red/2/pi/c                      % Hz
AoA_max = abs(atan2(pi*f*A_TE,Qinf) + alpha_max)*180/pi
% St = f*2*A_TE/Qinf
delT = 1/f/Nstep;           % s
% delT = 20/Nstep;
% ep = max([kappa*h_c*c*f^2*delT^2  2*kappa*A_TE*f^2*delT^2 1/5*Qinf*delT])
% ep = 1/5*Qinf*delT;
% ep = c*1e-3;
% epSC = ep;
% epBod = c*1e-6;
% epB = 0.01*c


NonDimTimeStep = delT*Qinf/c


%% Initializing Matrices
xTE = zeros(2,Ncyc*Nstep+1);
zTE = zeros(2,Ncyc*Nstep+1);
xw = zeros(1,Ncyc*Nstep+1);
zw = zeros(1,Ncyc*Nstep+1);
xpbod = zeros(Npanels + 1,Ncyc*Nstep + 1);
zpbod = zeros(Npanels + 1,Ncyc*Nstep + 1);

xpt_LES = zeros(2,Ncyc*Nstep);
zpt_LES = zeros(2,Ncyc*Nstep);
xpb_LES = zeros(2,Ncyc*Nstep);
zpb_LES = zeros(2,Ncyc*Nstep);
Sep_t_Store = ones(1,Ncyc*Nstep+1);
Sep_b_Store = ones(1,Ncyc*Nstep+1);
Sep_t = 0;
Sep_b = 0;


vt = zeros(Npanels,2,Ncyc*Nstep+1);
vn = zeros(Npanels,2,Ncyc*Nstep+1);
vtTE = zeros(1,2,Ncyc*Nstep+1);
vnTE = zeros(1,2,Ncyc*Nstep+1);
vtw = zeros(Ncyc*Nstep,2);
vnw = zeros(Ncyc*Nstep,2);
vtLEt = zeros(Ncyc*Nstep,2);
vnLEt = zeros(Ncyc*Nstep,2);
vtLEb = zeros(Ncyc*Nstep,2);
vnLEb = zeros(Ncyc*Nstep,2);

dL = zeros(Npanels,Ncyc*Nstep+1);

mu = zeros(Npanels,Ncyc*Nstep+1);
muTE = zeros(1,Ncyc*Nstep+1);
muW = zeros(1,Ncyc*Nstep);
muLEt = zeros(1,Ncyc*Nstep);
muLEb = zeros(1,Ncyc*Nstep);
sigma = zeros(Npanels,Ncyc*Nstep+1);

Vp = zeros(Npanels+1,2,Ncyc*Nstep+1);
Vbod = zeros(Npanels+1,2,Ncyc*Nstep+1);
Vc = zeros(Npanels,2,Ncyc*Nstep+1);
Vt = zeros(Npanels,2,Ncyc*Nstep+1);
Vshear_t = zeros(2,Ncyc*Nstep+1);
Vshear_b = zeros(2,Ncyc*Nstep+1);
Qp = zeros(Npanels,Ncyc*Nstep+1);
Qt = zeros(Npanels,Ncyc*Nstep+1);

Cw = zeros(Npanels,Ncyc*Nstep);
dPhi_wd = zeros(Ncyc*Nstep,Npanels);
dPhi_wd_2 = zeros(Ncyc*Nstep,Npanels);
CLEwt = zeros(Npanels,Ncyc*Nstep-1);
CLEwb = zeros(Npanels,Ncyc*Nstep-1);


Cp_s = zeros(Npanels,Ncyc*Nstep+1);
Cp_us = zeros(Npanels,Ncyc*Nstep+1);
Cp = zeros(Npanels,Ncyc*Nstep+1);

dFshear = zeros(Npanels,Ncyc*Nstep+1);
delF = zeros(Npanels,2,Ncyc*Nstep+1);
delFp = zeros(Npanels,2,Ncyc*Nstep+1);
delFs = zeros(Npanels,2,Ncyc*Nstep+1);
delP = zeros(Npanels,Ncyc*Nstep+1);
F = zeros(1,2,Ncyc*Nstep+1);
Fx = zeros(1,Ncyc*Nstep+1);
Fz = zeros(1,Ncyc*Nstep+1);
L = zeros(1,Ncyc*Nstep+1);
T = zeros(1,Ncyc*Nstep+1);
Pow = zeros(1,Ncyc*Nstep+1);

Gamma = zeros(1,Ncyc*Nstep+1);
GammaW = zeros(1,Ncyc*Nstep+1);

Cf = zeros(1,Ncyc*Nstep+1);
Cl = zeros(1,Ncyc*Nstep+1);
Ct = zeros(1,Ncyc*Nstep+1);
Cpow = zeros(1,Ncyc*Nstep+1);

% Initializing matrices for the ground image.
if grd == 1
    xTE_2 = zeros(2,Ncyc*Nstep+1);
    zTE_2 = zeros(2,Ncyc*Nstep+1);
    xw_2 = zeros(1,Ncyc*Nstep+1);
    zw_2 = zeros(1,Ncyc*Nstep+1);

    vt_2 = zeros(Npanels,2,Ncyc*Nstep+1);
    vn_2 = zeros(Npanels,2,Ncyc*Nstep+1);
    vtTE_2 = zeros(1,2,Ncyc*Nstep+1);
    vnTE_2 = zeros(1,2,Ncyc*Nstep+1);
    vtw_2 = zeros(Ncyc*Nstep,2);
    vnw_2 = zeros(Ncyc*Nstep,2);
end


%% Geometry
% Specifying the surface shape function (Van de Vooren airfoil) and panel 
% corner points.
% [xpt,zpt,xpb,zpb,~,~,~,~,~,~,~,~] = VanDeVooren(c,Npanels+1,tmax/2);
[xpt,zpt,xpb,zpb,~,~,~,~,~,~,~,~] = NACA(c,Npanels+1,tmax);
xp = [xpb;xpt(2:end)];
zp = [zpb(1:end-1)-zpb(1);0;zpt(2:end)-zpt(end)];

[~,Indsep_tmax] = min(abs(xp(Npanels/2 + 1:end)  - 0.05*c)); 
Indsep_tmax = Indsep_tmax + Npanels/2;
[~,Indsep_bmax] = min(abs(xp(1:Npanels/2 + 1) - 0.05*c));
    


% Locations of collocation points (mid-panel on the surface).
xc = (xp(2:end,1) + xp(1:end-1,1))/2;
zc = (zp(2:end,1) + zp(1:end-1,1))/2;
zcval = zc(Npanels,1);

% Copying panel points to each timestep.
xp = kron(xp,ones(1,Ncyc*Nstep+1));
zp = kron(zp,ones(1,Ncyc*Nstep+1));
xc = kron(xc,ones(1,Ncyc*Nstep+1));
zc = kron(zc,ones(1,Ncyc*Nstep+1));


% Normal vectors for the panels at collocation points.  Each row is a new
% collocation point.
vttemp = [(xp(2:end,1) - xp(1:end-1,1)) (zp(2:end,1) - zp(1:end-1,1))];
vt(:,:,1) = diag(1./sqrt(vttemp(:,1).^2 + vttemp(:,2).^2))*vttemp;
vn(:,:,1) = [-vt(:,2,1) vt(:,1,1)];


% Calculating ground image panels.
if grd == 1
    xc_2 = xc;
    zc_2  = -zc;

    xp_2  = xp;
    zp_2  = -zp;

    vt_2(:,:,1) = [vt(:,1,1) -vt(:,2,1)];
    vn_2(:,:,1) = [vn(:,1,1) -vn(:,2,1)];
end

%% Flight trajectory: Translation
% Free-stream velocity
Uinf = Qinf*cos(alpha);
Winf = Qinf*sin(alpha);

if free == 1
    Q0 = kron(-[Uinf 0]',ones(1,Ncyc*Nstep+1));
else
    Q0 = kron(-[Uinf Winf]',ones(1,Ncyc*Nstep+1));
end
x_b = zeros(1,Ncyc*Nstep+1);
z_b = zeros(1,Ncyc*Nstep+1);
a_b = zeros(1,Ncyc*Nstep+1);

%% Time-stepping solution
scrsz = get(0,'ScreenSize');
wbar = waitbar(0,'Calculating...');
delTime = zeros(Ncyc*Nstep+1,1);

for i_t = 1:Ncyc*Nstep+1
    
    % Start timer.
    tic

    t = (i_t - 1)*delT;
    
    if t > 0
        muW(1,Ncyc*Nstep + 2 - i_t) = muTE(1,i_t-1);
    end
    

      
    % Calculating shift in panel positions with the swimming velocity. 
    delxp = x_b(1,i_t)*ones(Npanels+1,1);
    delzp = (z_b(1,i_t) + D)*ones(Npanels+1,1);
    delxc = x_b(1,i_t)*ones(Npanels,1);
    delzc = (z_b(1,i_t) + D)*ones(Npanels,1);
    
    
    
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
    [xptemp,zptemp,Vptemp] = KinematicsHeavePitch2D(xp(:,i_t),zp(:,i_t),alpha_max,h_c,f,t,phi,phase); 
    Vp(:,:,i_t) = Vptemp';
    xpbod(:,i_t) = xptemp';
    zpbod(:,i_t) = zptemp';
    xp(:,i_t) = xptemp';
    zp(:,i_t) = zptemp';
    Vc(:,:,i_t) = 1/2*Vp(2:end,:,i_t) + 1/2*Vp(1:end-1,:,i_t);
    if i_t > 1
        Vbod(:,1,i_t) = (xpbod(:,i_t) - xpbod(:,i_t-1))/delT;
        Vbod(:,2,i_t) = (zpbod(:,i_t) - zpbod(:,i_t-1))/delT;
    else
        Vbod(:,:,i_t) = Vp(:,:,i_t);
    end
    
    % Calculating locations of collocation points (mid-panel on the surface).
    xc(:,i_t) = (xp(2:end,i_t) + xp(1:end-1,i_t))/2;
    zc(:,i_t) = (zp(2:end,i_t) + zp(1:end-1,i_t))/2;
    
    % Superposing the kinematics and swimming translations.
    xp(:,i_t) = xp(:,i_t) + delxp;
    zp(:,i_t) = zp(:,i_t) + delzp;
    xc(:,i_t) = xc(:,i_t) + delxc;
    zc(:,i_t) = zc(:,i_t) + delzc;
   
    % Normal vectors for the panels at collocation points.  Each row is a new
    % collocation point.
    vttemp = [(xp(2:end,i_t) - xp(1:end-1,i_t)) (zp(2:end,i_t) - zp(1:end-1,i_t))];
    vt(:,:,i_t) = diag(1./sqrt(vttemp(:,1).^2 + vttemp(:,2).^2))*vttemp;
    vn(:,:,i_t) = [-vt(:,2,i_t) vt(:,1,i_t)];
    
    
    % Calculating ground image panels.
    if grd == 1
        xc_2(:,i_t) = xc(:,i_t);
        zc_2(:,i_t)  = -zc(:,i_t);
        
        xp_2(:,i_t)  = xp(:,i_t);
        zp_2(:,i_t)  = -zp(:,i_t);
        
        vt_2(:,:,i_t) = [vt(:,1,i_t) -vt(:,2,i_t)];
        vn_2(:,:,i_t) = [vn(:,1,i_t) -vn(:,2,i_t)];
    end
    
    
    %% Applying fencing after updating the kinematics
    if LES == 1 && t > 0

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
   

    %% Calculating source strengths.  
    Vt(:,:,i_t) = kron(ones(Npanels,1),Q0(:,i_t)') + Vc(:,:,i_t);
    sigma(:,i_t) = sum(vn(:,:,i_t).*Vt(:,:,i_t),2);    
    
    %% Calculating the new trailing edge panels.
    xp_avg = xp(end-1,i_t)/2 + xp(2,i_t)/2;
    zp_avg = zp(end-1,i_t)/2 + zp(2,i_t)/2;
    TEvec = [xp(end,i_t) - xp_avg; zp(end,i_t) - zp_avg]/norm([xp(end,i_t) - xp_avg; zp(end,i_t) - zp_avg]);
    xTE(:,i_t) = [xp(end,i_t); xp(end,i_t) + TEvec(1)*Qinf*delT*TEfac];
    zTE(:,i_t) = [zp(end,i_t); zp(end,i_t) + TEvec(2)*Qinf*delT*TEfac];

%     TEvec = -([Q0(1,end); 0*Q0(2,end)] + Vp(end,:,i_t)')/norm(([Q0(1,end); 0*Q0(2,end)]+ Vp(end,:,i_t)'));
%     xTE(:,i_t) = [xp(end,i_t); xp(end,i_t) + TEvec(1)*Qinf*delT*TEfac];
%     zTE(:,i_t) = [zp(end,i_t); zp(end,i_t) + TEvec(2)*Qinf*delT*TEfac];
%     xTE(:,i_t) = [xp(end,i_t); xp(end,i_t) + TEvec(1)*c*TEfac];
%     zTE(:,i_t) = [zp(end,i_t); zp(end,i_t) + TEvec(2)*c*TEfac];
    
    vttemp = [(xTE(2,i_t) - xTE(1,i_t)) (zTE(2,i_t) - zTE(1,i_t))];
    vtTE(1,:,i_t) = diag(1./sqrt(vttemp(:,1).^2 + vttemp(:,2).^2))*vttemp;
    vnTE(1,:,i_t) = [-vtTE(1,2,i_t) vtTE(1,1,i_t)];

    % Calculating ground image TE panels.
    if grd == 1
        xTE_2(:,i_t) = xTE(:,i_t);
        zTE_2(:,i_t)  = -zTE(:,i_t);
        
        vtTE_2(:,:,i_t) = [vtTE(:,1,i_t) -vtTE(:,2,i_t)];
        vnTE_2(:,:,i_t) = [vnTE(:,1,i_t) -vnTE(:,2,i_t)];
    end
    
    % Shedding a wake panel after first time step.
    if t > 0
        
        % Calculating the downstream position of the newly shed wake panel
        xw(1,Ncyc*Nstep + 3 - i_t) = xTE(2,i_t-1);
        zw(1,Ncyc*Nstep + 3 - i_t) = zTE(2,i_t-1);
        
        % Calculating the upstream position of the newly shed wake panel
        xw(1,Ncyc*Nstep + 2 - i_t) = xTE(2,i_t);
        zw(1,Ncyc*Nstep + 2 - i_t) = zTE(2,i_t);
        
        vttemp = [(xw(1,Ncyc*Nstep + 3 - i_t:end) - xw(1,Ncyc*Nstep + 2 - i_t:end-1))' (zw(1,Ncyc*Nstep + 3 - i_t:end) - zw(1,Ncyc*Nstep + 2 - i_t:end-1))'];
        vtw(Ncyc*Nstep + 2 - i_t:end,:) = diag(1./sqrt(vttemp(:,1).^2 + vttemp(:,2).^2))*vttemp;
        vnw(Ncyc*Nstep + 2 - i_t:end,:) = [-vtw(Ncyc*Nstep + 2 - i_t:end,2) vtw(Ncyc*Nstep + 2 - i_t:end,1)]; 

        % Calculating ground image wake panels.
        if grd == 1
            xw_2(1,Ncyc*Nstep + 2 - i_t:Ncyc*Nstep + 3 - i_t) = xw(1,Ncyc*Nstep + 2 - i_t:Ncyc*Nstep + 3 - i_t);
            zw_2(1,Ncyc*Nstep + 2 - i_t:Ncyc*Nstep + 3 - i_t) = -zw(1,Ncyc*Nstep + 2 - i_t:Ncyc*Nstep + 3 - i_t);

            vtw_2 = [vtw(:,1) -vtw(:,2)];
            vnw_2 = [vnw(:,1) -vnw(:,2)];
        end
    end
       
    
    % Moving collocation points inward.
    val = 0.25;
    InVal = abs(val*zcval);
    Xc = xc(:,i_t) - InVal.*vn(:,1,i_t);
    Zc = zc(:,i_t) - InVal.*vn(:,2,i_t);
    
%     Xc = xc(:,i_t) - abs(val*zcval).*vn(:,1,i_t);
%     Zc = zc(:,i_t) - abs(val*zcval).*vn(:,2,i_t);


    %% Calculating influence coefficients 
    % (a = q*n, q = [u, w]') at each collocation point.  The row denotes the collocation
    % point while the column denotes the doublet element.  Also calculating RHS,
    % RHS = -Qinf*n.
    
    % Influence of the body panels
    dPhi_s_2 = zeros(Npanels,Npanels);
    dPhi_d_2 = zeros(Npanels,Npanels);
    
    [dPhi_s,dPhi_d,dLtemp] = Phi(ones(1,Npanels),ones(1,Npanels),Xc',Zc',xp(:,i_t)',zp(:,i_t)',vt(:,:,i_t)',vn(:,:,i_t)');
    dL(:,i_t) = dLtemp;
    
    if grd == 1
        [dPhi_s_2,dPhi_d_2,dLtemp2] = Phi(ones(1,Npanels),ones(1,Npanels),Xc',Zc',xp_2(:,i_t)',zp_2(:,i_t)',vt_2(:,:,i_t)',vn_2(:,:,i_t)');
    end
    
    B = dPhi_s' + dPhi_s_2';
    Cb = dPhi_d' + dPhi_d_2';
    
     
    % Influence of the trailing edge panel
    dPhi_d_2 = zeros(1,Npanels);
    Cte = zeros(Npanels,Npanels);
    
    [~,dPhi_d,~] = Phi(1,0,Xc',Zc',xTE(:,i_t)',zTE(:,i_t)',vtTE(1,:,i_t)',vnTE(1,:,i_t)');
    
    if grd == 1
        [~,dPhi_d_2,~] = Phi(1,0,Xc',Zc',xTE_2(:,i_t)',zTE_2(:,i_t)',vtTE_2(1,:,i_t)',vnTE_2(1,:,i_t)');
    end
    
    Cte(:,Npanels) = dPhi_d' + dPhi_d_2';
    Cte(:,1) = - dPhi_d' - dPhi_d_2';
    
    
    % Influence of the wake panels
    if t > 0        
        [~,dPhi_d,~] = Phi(ones(1,length(muW(Ncyc*Nstep + 2 - i_t:end))),zeros(1,length(muW(Ncyc*Nstep + 2 - i_t:end))),Xc',Zc',xw(Ncyc*Nstep + 2 - i_t:end),zw(Ncyc*Nstep + 2 - i_t:end),vtw(Ncyc*Nstep + 2 - i_t:end,:)',vnw(Ncyc*Nstep + 2 - i_t:end,:)');
        
        dPhi_d_2 = 0*dPhi_d;
        
        if grd == 1
            [~,dPhi_d_2,~] = Phi(ones(1,length(muW(Ncyc*Nstep + 2 - i_t:end))),zeros(1,length(muW(Ncyc*Nstep + 2 - i_t:end))),Xc',Zc',xw_2(Ncyc*Nstep + 2 - i_t:end),zw_2(Ncyc*Nstep + 2 - i_t:end),vtw_2(Ncyc*Nstep + 2 - i_t:end,:)',vnw_2(Ncyc*Nstep + 2 - i_t:end,:)');
        end
        
        Cw(:,Ncyc*Nstep + 2 - i_t:end) = dPhi_d' + dPhi_d_2';
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


    %% Applying the Kutta condition, solving.
    if LES == 1 && i_t > 1
        A = Cb + Cte + CLEt + CLEb;
    else
        A = Cb + Cte;
    end
    
    % Defining Right Hand Side (RHS)
    if LES == 1 && i_t > 2 
        RHS = -B*sigma(:,i_t) - wakeinf*Cw*muW' - wakeinf*CLEwt*muLEt' - wakeinf*CLEwb*muLEb';
    elseif i_t > 1
        RHS = -B*sigma(:,i_t) - wakeinf*Cw*muW';
    else
        RHS = -B*sigma(:,i_t);
    end

    % Solving for doublet strengths.
    mu(:,i_t) = A\RHS;
    
    % Solving for TE panel strength.
    muTE(1,i_t) = mu(Npanels,i_t) - mu(1,i_t);
    
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
    
    % Calculating wake circulation.
    if i_t == 1
        GammaW(1,Ncyc*Nstep + 1) = muTE(1,i_t);
    else
        GammaW(1,Ncyc*Nstep + 2 - i_t) = muTE(1,i_t) - muTE(1,i_t - 1);
    end

    %% Calculating on-body velocities and pressures. 
    % Pressures are calculated from the unsteady Bernoulli's equation.  
    % The flux through the body should be zero.  
 
    % Local velocity over the body due to the perturbation potential.
    Qp(2:end-1,i_t) = (mu(1:end-2,i_t) - mu(3:end,i_t))./(dL(3:end,i_t)/2 + dL(2:end-1,i_t) + dL(1:end-2,i_t)/2);
      
    % First and last panel velocity calculation (explicitly imposing the Kutta condition)
%     Qp(1,i_t) = (mu(1,i_t) - mu(2,i_t))/(dL(2,i_t)/2 + 3/2*dL(1,i_t)); 
%     Qp(Npanels,i_t) = (mu(Npanels - 1,i_t) - mu(Npanels,i_t))/(dL(Npanels - 1,i_t)/2 + 3/2*dL(Npanels,i_t));
    Qp(1,i_t) = (mu(1,i_t) - mu(2,i_t))/(dL(2,i_t) + dL(1,i_t)); 
    Qp(Npanels,i_t) = (mu(Npanels - 1,i_t) - mu(Npanels,i_t))/(dL(Npanels - 1,i_t) + dL(Npanels,i_t));
    
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
    Qt(:,i_t) = Qp(:,i_t) + sum(-Vt(:,:,i_t).*vt(:,:,i_t),2);
%     Qt(1,i_t) = Qt(1,i_t) - Uinf*vt(1,1,i_t)/2 - Winf*vt(1,2,i_t)/2;
%     Qt(Npanels,i_t) = Qt(Npanels,i_t) - Uinf*vt(Npanels,1,i_t)/2 - Winf*vt(Npanels,2,i_t)/2;
%     if LES == 1 && t > 0 
%         if Indsep_t < Npanels + 1
%             Qt(Indsep_t,i_t) = Qt(Indsep_t,i_t) - Uinf*vt(Indsep_t,1,i_t)/2 - Winf*vt(Indsep_t,2,i_t)/2;
%             Qt(Indsep_t-1,i_t) = Qt(Indsep_t-1,i_t) - Uinf*vt(Indsep_t-1,1,i_t)/2 - Winf*vt(Indsep_t-1,2,i_t)/2;
%         end
%         if Indsep_b > 1
%             Qt(Indsep_b,i_t) = Qt(Indsep_b,i_t) - Uinf*vt(Indsep_b,1,i_t)/2 - Winf*vt(Indsep_b,2,i_t)/2;
%             Qt(Indsep_b-1,i_t) = Qt(Indsep_b-1,i_t) - Uinf*vt(Indsep_b-1,1,i_t)/2 - Winf*vt(Indsep_b-1,2,i_t)/2;
%         end
%     end
    % Calculating the pressure coefficient.
    Cp_s(:,i_t) = 1 - Qt(:,i_t).^2/Qinf^2;
    
%     if t == 0
%         Cp_us(:,i_t) = 2/Qinf^2*mu(:,i_t)'/delT;
%     else
%         Cp_us(:,i_t) = 2/Qinf^2*(mu(:,i_t) - mu(:,i_t-1))'/delT;
%     end


%% Calculating Phi on the surface of each panel.

% % Moving collocation points outward.
%     val = 1e-5;
%     InVal = abs(val*dL(:,i_t));
%     Xc = xc(:,i_t) + InVal.*vn(:,1,i_t);
%     Zc = zc(:,i_t) + InVal.*vn(:,2,i_t);
%       
%     % Velocity Potential due to the body panels
%     dPhi_bs_2 = zeros(Npanels,Npanels);
%     dPhi_bd_2 = zeros(Npanels,Npanels);
%     
%     [dPhi_bs,dPhi_bd,~] = Phi(mu(:,i_t)',sigma(:,i_t)',Xc',Zc',xp(:,i_t)',zp(:,i_t)',vt(:,:,i_t)',vn(:,:,i_t)');
%     
%     if grd == 1
%         [dPhi_bs_2,dPhi_bd_2,~] = Phi(mu(:,i_t)',sigma(:,i_t)',Xc',Zc',xp_2(:,i_t)',zp_2(:,i_t)',vt_2(:,:,i_t)',vn_2(:,:,i_t)');
%     end
%     
%      
%     % Velocity Potential due to the trailing edge panel
%     dPhi_td_2 = zeros(1,Npanels);
%     
%     [~,dPhi_td,~] = Phi(muTE(i_t),0,Xc',Zc',xTE(:,i_t)',zTE(:,i_t)',vtTE(1,:,i_t)',vnTE(1,:,i_t)');
%     
%     if grd == 1
%         [~,dPhi_td_2,~] = Phi(muTE(i_t),0,Xc',Zc',xTE_2(:,i_t)',zTE_2(:,i_t)',vtTE_2(1,:,i_t)',vnTE_2(1,:,i_t)');
%     end
%     
%     
%     % Velocity Potential due to the wake panels
%     if t > 0        
%         [~,dPhi_wd,~] = Phi(muW(Ncyc*Nstep + 2 - i_t:end),zeros(1,length(muW(Ncyc*Nstep + 2 - i_t:end))),Xc',Zc',xw(Ncyc*Nstep + 2 - i_t:end),zw(Ncyc*Nstep + 2 - i_t:end),vtw(Ncyc*Nstep + 2 - i_t:end,:)',vnw(Ncyc*Nstep + 2 - i_t:end,:)');
%         
%         dPhi_wd_2 = 0*dPhi_wd;
%         
%         if grd == 1
%             [~,dPhi_wd_2,~] = Phi(muW(Ncyc*Nstep + 2 - i_t:end),zeros(1,length(muW(Ncyc*Nstep + 2 - i_t:end))),Xc',Zc',xw_2(Ncyc*Nstep + 2 - i_t:end),zw_2(Ncyc*Nstep + 2 - i_t:end),vtw_2(Ncyc*Nstep + 2 - i_t:end,:)',vnw_2(Ncyc*Nstep + 2 - i_t:end,:)');
%         end
%         
%     end
    
%     dPhi(:,i_t) = (sum(dPhi_bs' + dPhi_bd',2) + dPhi_td' + sum(dPhi_wd',2));
%     dPhi(:,i_t) = (sum(dPhi_bs',2) + sum(dPhi_bd',2) + dPhi_td' + sum(dPhi_wd',2) + sum(dPhi_bs_2',2) + sum(dPhi_bd_2',2) + dPhi_td_2' + sum(dPhi_wd_2',2));
%     dPhi(:,i_t) = 2*(-mu(:,i_t)/2 + sigma(:,i_t)/4/pi.*dL(:,i_t).*log((dL(:,i_t)/2).^2));

    
%     if LES == 1 && t > 0
%         
%         % Influence of the LEtop edge panel
%         CLEt = zeros(Npanels,Npanels);
%         if Indsep_t == Npanels + 1
%         else
%             [~,dPhi_d,~] = Phi(1,0,Xc',Zc',xpt_LES(:,Ncyc*Nstep + 2 - i_t)',zpt_LES(:,Ncyc*Nstep + 2 - i_t)',vtLEt(Ncyc*Nstep + 2 - i_t,:)',vnLEt(Ncyc*Nstep + 2 - i_t,:)');
%             CLEt(:,Indsep_t-1) = CircRed*dPhi_d';
%             CLEt(:,Indsep_t) = -CircRed*dPhi_d';
%         end
% 
%         % Influence of the LEbot edge panel
%         CLEb = zeros(Npanels,Npanels);
%         if Indsep_b == 1
%         else
%             [~,dPhi_d,~] = Phi(1,0,Xc',Zc',xpb_LES(:,Ncyc*Nstep + 2 - i_t)',zpb_LES(:,Ncyc*Nstep + 2 - i_t)',vtLEb(Ncyc*Nstep + 2 - i_t,:)',vnLEb(Ncyc*Nstep + 2 - i_t,:)');
%             CLEb(:,Indsep_b) = -CircRed*dPhi_d';
%             CLEb(:,Indsep_b-1) = CircRed*dPhi_d';
%         end
%         
%         if i_t > 2
%             for j = i_t:-1:3
%                 % Influence of the LE top wake sheet panels
%                 [~,dPhi_d,~] = Phi(ones(1,length(muLEt(Ncyc*Nstep + 3 - j))),zeros(1,length(muLEt(Ncyc*Nstep + 3 - j))),Xc',Zc',xpt_LES(:,Ncyc*Nstep + 3 - j)',zpt_LES(:,Ncyc*Nstep + 3 - j)',vtLEt(Ncyc*Nstep + 3 - j,:)',vnLEt(Ncyc*Nstep + 3 - j,:)');
%                 CLEwt(:,Ncyc*Nstep + 3 - j) = dPhi_d';
% 
%                 % Influence of the LE top wake sheet panels
%                 [~,dPhi_d,~] = Phi(ones(1,length(muLEb(Ncyc*Nstep + 3 - j))),zeros(1,length(muLEb(Ncyc*Nstep + 3 - j))),Xc',Zc',xpb_LES(:,Ncyc*Nstep + 3 - j)',zpb_LES(:,Ncyc*Nstep + 3 - j)',vtLEb(Ncyc*Nstep + 3 - j,:)',vnLEb(Ncyc*Nstep + 3 - j,:)');
%                 CLEwb(:,Ncyc*Nstep + 3 - j) = dPhi_d';
%             end
%         end
%     end
    
    %%    
    
%     if i_t == 1
%         Cp_us(:,i_t) = (-2/Qinf^2)*dPhi(:,i_t)'/delT;
%     elseif i_t == 2
%         Cp_us(:,i_t) = (-2/Qinf^2)*(dPhi(:,i_t) - dPhi(:,i_t-1))'/delT;
%     else
%         Cp_us(:,i_t) = (-2/Qinf^2)*(3*dPhi(:,i_t) - 4*dPhi(:,i_t-1) + dPhi(:,i_t-2))'/2/delT;
%     end
%     Cp(:,i_t) = Cp_s(:,i_t) + Cp_us(:,i_t);    
    
   
    if i_t == 1
        Cp_us(:,i_t) = 2/Qinf^2*mu(:,i_t)'/delT;
    elseif i_t == 2
        Cp_us(:,i_t) = 2/Qinf^2*(mu(:,i_t) - mu(:,i_t-1))'/delT;
    else
        Cp_us(:,i_t) = 2/Qinf^2*(3*mu(:,i_t) - 4*mu(:,i_t-1) + mu(:,i_t-2))'/2/delT;
    end
    Cp(:,i_t) = Cp_s(:,i_t) + Cp_us(:,i_t);
%     
%     if t == 0
%         Cp_us(:,i_t) = (-2/Qinf^2)*(1/2)*(mu(:,i_t)'/delT - sigma(:,i_t)'/delT);
%     else
%         Cp_us(:,i_t) = (-2/Qinf^2)*(1/2)*((mu(:,i_t) - mu(:,i_t-1))'/delT - (sigma(:,i_t) - sigma(:,i_t-1))'/delT);
%     end
%     Cp(:,i_t) = Cp_s(:,i_t) + Cp_us(:,i_t);
    
%     if t == 0
%         Cp_us(:,i_t) = -2/Qinf^2*Phi_VP(:,i_t)'/delT;
%     else
%         Cp_us(:,i_t) = -2/Qinf^2*(Phi_VP(:,i_t) - Phi_VP(:,i_t-1))'/delT;
%     end
%     Cp(:,i_t) = Cp_s(:,i_t) + Cp_us(:,i_t);


    % Calculating viscous drag correction.  
    if ViscDrag
        [dFi_pct,dL_pct,factor,xbc_sep,xtc_sep,Indsep_b,Indsep_t,Stagpt] = SkinFrictionSolver(xc(:,i_t),Qt(:,i_t),nu,rho,c,tmax,Cp(:,i_t),Qinf);   
        dFi = dFi_pct.*kron(dL(:,i_t)*b,ones(factor,1));
        dFshear(:,i_t) = kron(eye(Npanels),ones(1,factor))*dFi;
    else
        dFshear(:,i_t) = zeros(Npanels,1);
    end


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

%     if Indsep_t < Stagpt + 5
%         Indsep_t = Stagpt + 5;
%     end
%     
%     if Indsep_b > Stagpt - 5
%         Indsep_b = Stagpt - 5;
%     end
    
    % Calculating lift on the body and pitching moment about the leading
    % edge.
    delFp(:,:,i_t) = -kron((Cp(:,i_t)*1/2*rho*Qinf^2.*dL(:,i_t)*b),[1 1]).*vn(:,:,i_t);
    delFs(:,:,i_t) = ViscDrag*kron(dFshear(:,i_t),[1 1]).*vt(:,:,i_t);
    delF(:,:,i_t) = delFp(:,:,i_t) + delFs(:,:,i_t);
%     delP(:,i_t) = sum(-delF(:,:,i_t).*Vc(:,:,i_t),2);

    delP(:,i_t) = sum(-delF(:,:,i_t).*Vc(:,:,i_t),2);
    
    
    F(1,:,i_t) = sum(delF(:,:,i_t),1);
    Pow(1,i_t) = abs(sum(delP(:,i_t),1));
    Fx(1,i_t) = F(1,1,i_t);
    Fz(1,i_t) = F(1,2,i_t);
    
    L(1,i_t) = Fz(1,i_t)*cos(alpha) - Fx(1,i_t)*sin(alpha);
    T(1,i_t) = -(Fz(1,i_t)*sin(alpha) + Fx(1,i_t)*cos(alpha));
    
    % M0 = -sum(delL*cos(alpha).*xp(2:end-1));

    % Calculating non-dimensional coefficients.
    % Cm = M0/(1/2*rho*Qinf^2*c^2)
    Gamma(1,i_t) = mu(1,i_t) - mu(Npanels,i_t);

    Cf(1,i_t) = norm(F(1,:,i_t))/(1/2*rho*Qinf^2*c*b);
    Cl(1,i_t) = L(1,i_t)/(1/2*rho*Qinf^2*c*b);
    Ct(1,i_t) = T(1,i_t)/(1/2*rho*Qinf^2*c*b);
    Cpow(1,i_t) = Pow(1,i_t)/(1/2*rho*Qinf^3*c*b);
    
%     Gamma_ideal = pi*c*Qinf*alpha;
%     Cl_ideal = 2*pi*sin(alpha);
%     Cl_KJ = rho*Qinf*Gamma/(1/2*rho*Qinf^2*c);
    

    %% Wake rollup calculations
    if Rollup == 1 && t > 0
                
        [u,w] = WakeRollupVelocity2D(xw(Ncyc*Nstep + 3 - i_t:end),zw(Ncyc*Nstep + 3 - i_t:end),mu(:,i_t)',muTE(i_t),muW(Ncyc*Nstep + 2 - i_t:end),muLEt,muLEb,sigma(:,i_t)',xp(:,i_t)',zp(:,i_t)',vt(:,:,i_t)',vn(:,:,i_t)',xTE(:,i_t)',zTE(:,i_t)',vtTE(1,:,i_t)',vnTE(1,:,i_t)',xw(Ncyc*Nstep + 2 - i_t:end),zw(Ncyc*Nstep + 2 - i_t:end),vtw(Ncyc*Nstep + 2 - i_t:end,:)',vnw(Ncyc*Nstep + 2 - i_t:end,:)',xpt_LES,zpt_LES,vtLEt,vnLEt,xpb_LES,zpb_LES,vtLEb,vnLEb,i_t,Ncyc,Nstep,LES,ep,epBod,SC);
        
        u_2 = 0*u;
        w_2 = 0*w;
        
        if grd == 1
            [u_2,w_2] = WakeRollupVelocity2D(xw(Ncyc*Nstep + 3 - i_t:end),zw(Ncyc*Nstep + 3 - i_t:end),mu(:,i_t)',muTE(i_t),muW(Ncyc*Nstep + 2 - i_t:end),muLEt,muLEb,sigma(:,i_t)',xp_2(:,i_t)',zp_2(:,i_t)',vt_2(:,:,i_t)',vn_2(:,:,i_t)',xTE_2(:,i_t)',zTE_2(:,i_t)',vtTE_2(1,:,i_t)',vnTE_2(1,:,i_t)',xw_2(Ncyc*Nstep + 2 - i_t:end),zw_2(Ncyc*Nstep + 2 - i_t:end),vtw_2(Ncyc*Nstep + 2 - i_t:end,:)',vnw_2(Ncyc*Nstep + 2 - i_t:end,:)',xpt_LES,zpt_LES,vtLEt,vnLEt,xpb_LES,zpb_LES,vtLEb,vnLEb,i_t,Ncyc,Nstep,LES,ep,epBod,SC);
        end
        
        % Applying fencing scheme
        [xstar_w,zstar_w,fence] = fencing(c,Vc(:,:,i_t),xw(Ncyc*Nstep + 3 - i_t:end),zw(Ncyc*Nstep + 3 - i_t:end),u + u_2,w + w_2,xp(:,i_t),zp(:,i_t),vn(:,:,i_t),x_b(1,i_t),z_b(1,i_t),delT,epB);
        xstar_w = xstar_w.*fence + (xstar_w + (u + u_2)*delT).*(~fence);
        zstar_w = zstar_w.*fence + (zstar_w + (w + w_2)*delT).*(~fence);
        
        % Wake rollup
        xw(Ncyc*Nstep + 3 - i_t:end) = xstar_w;
        zw(Ncyc*Nstep + 3 - i_t:end) = zstar_w;
       
        vttemp = [(xw(1,Ncyc*Nstep + 3 - i_t:end) - xw(1,Ncyc*Nstep + 2 - i_t:end-1))' (zw(1,Ncyc*Nstep + 3 - i_t:end) - zw(1,Ncyc*Nstep + 2 - i_t:end-1))'];
        vtw(Ncyc*Nstep + 2 - i_t:end,:) = diag(1./sqrt(vttemp(:,1).^2 + vttemp(:,2).^2))*vttemp;
        vnw(Ncyc*Nstep + 2 - i_t:end,:) = [-vtw(Ncyc*Nstep + 2 - i_t:end,2) vtw(Ncyc*Nstep + 2 - i_t:end,1)]; 

        if grd == 1
            % Wake rollup
            xw_2(Ncyc*Nstep + 3 - i_t:end) = xw(Ncyc*Nstep + 3 - i_t:end);
            zw_2(Ncyc*Nstep + 3 - i_t:end) = -zw(Ncyc*Nstep + 3 - i_t:end);
        
            vtw_2(Ncyc*Nstep + 2 - i_t:end,:) = [vtw(Ncyc*Nstep + 2 - i_t:end,1)  -vtw(Ncyc*Nstep + 2 - i_t:end,2)];
            vnw_2(Ncyc*Nstep + 2 - i_t:end,:) = [vnw(Ncyc*Nstep + 2 - i_t:end,1)  -vnw(Ncyc*Nstep + 2 - i_t:end,2)];
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
    Sep_t_Store(1,i_t) = Indsep_t;
    Sep_b_Store(1,i_t) = Indsep_b;

    if t > 0
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
    
if t > 0
    %     % Closing previous figure 
    if i_t > 2
        close(fighand)
    end
    
    FontSizeAx = 24;
    FontSizeLb = 32;
    afFigurePosition = [0 2 30 10];
    axespos = [0.15 0.2 0.84 0.8];
    ylabelpos = [-0.12 0.5];
    xlabelpos = [0.5 -0.25];


    % Plotting airfoil with LE vortex sheets and the TE vortex sheet.
    fighand = figure;
    set(gcf, 'Units', 'centimeters','PaperPositionMode', 'auto','Position', afFigurePosition);
    set(gcf,'DefaultAxesFontSize',FontSizeAx,'DefaultAxesFontName','TimesNewRoman','DefaultAxesGridLineStyle','-.','DefaultAxesLineWidth',2,'DefaultAxesFontWeight','Normal')

    hold on
    axis equal
    plot(xp(:,i_t),zp(:,i_t),'-k','linewidth',2)
    plot(xTE(:,i_t),zTE(:,i_t),'.-b','linewidth',2)
    val = 1/10;
    
    if grd == 1
%         plot(xp_2(:,i_t),zp_2(:,i_t),'-k','linewidth',2)
%         plot(xTE_2(:,i_t),zTE_2(:,i_t),'.-b','linewidth',2)
        plot([min(xp(:,i_t))-10*c; min(xp(:,i_t))+10*c],[0 0],'-k','linewidth',4)
%         if t > 0
%             plot(xw_2(Ncyc*Nstep + 2 - i_t:end),zw_2(Ncyc*Nstep + 2 - i_t:end),'.-b','linewidth',2,'markersize',14)
%         end
    end

    if t > 0
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
        vec = vec(end:-1:1,:);
        
        WakeCirc = GammaW'/max(GammaW)/(1/2);
        WakeCirc(WakeCirc > 1) = 1;  
        WakeCirc(WakeCirc < -1) = -1;
        WakeCirc = WakeCirc + 1;
        WakeCirc = round(WakeCirc*100) + 1;
        WakeCirc(WakeCirc > 199) = 199;
        
        for i_w = Ncyc*Nstep + 2 - i_t:Ncyc*Nstep+1
            plot(xw(i_w),zw(i_w),'.','color',vec(WakeCirc(i_w),:),'linewidth',2,'markersize',14)
        end
    end

    axis equal
%     axis([-c/2 + x_b(i_t) 2*c + x_b(i_t) -1/2*c + z_b(i_t) 1/2*c + z_b(i_t)])
    axis([-c/2 + x_b(i_t) 4*c + x_b(i_t) -c/10 c])

    set(gca, 'FontName', 'TimesNewRoman', 'FontSize', FontSizeAx)
    set(gca, 'Units', 'normalized', 'Position', axespos);
    % legend('d/c = \infty','d/c = 1/2','d/c = 1/4','Location','NorthWest')
    xlabel('$$x$$ (m)','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 0,'Units', 'Normalize','Position',xlabelpos,'Interpreter', 'LaTeX');
    ylabel('$$z$$ (m)','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 90,'Units', 'Normalize','Position',ylabelpos,'Interpreter', 'LaTeX');


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
        nx = 50;
        nz = 50;

        U = Uinf*ones(nz,nx);
        W = Winf*ones(nz,nx);

%         xf = linspace(-c/2 + x_b(i_t),3*c + x_b(i_t),nx)';
%         zf = linspace(-3*c/4 + z_b(i_t),3*c/4 + z_b(i_t),nz)';
        
        % Flow field for ground effect calculations
        xf = linspace(-c/4 + x_b(i_t),2*c + x_b(i_t),nx)';
        zf = linspace(0,c,nz)';

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
        [u_st,w_st,u_bdt,w_bdt] = DubSorV(mu(:,i_t)',sigma(:,i_t)',X,Z,xp(:,i_t)',zp(:,i_t)',vt(:,:,i_t)',vn(:,:,i_t)',epSC,SC);

        u_st_2 = 0*u_st;
        w_st_2 = 0*w_st;
        u_bdt_2 = 0*u_bdt;
        w_bdt_2 = 0*w_bdt;
        
        if grd == 1
            [u_st_2,w_st_2,u_bdt_2,w_bdt_2] = DubSorV(mu(:,i_t)',sigma(:,i_t)',X,Z,xp_2(:,i_t)',zp_2(:,i_t)',vt_2(:,:,i_t)',vn_2(:,:,i_t)',epSC,SC);
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
        [~,~,u_TEdt,w_TEdt] = DubSorV(muTE(i_t),0*muTE(i_t),X,Z,xTE(:,i_t)',zTE(:,i_t)',vtTE(1,:,i_t)',vnTE(1,:,i_t)',epSC,SC);

        u_TEdt_2 = 0*u_TEdt;
        w_TEdt_2 = 0*w_TEdt;
        
        if grd == 1
            [~,~,u_TEdt_2,w_TEdt_2] = DubSorV(muTE(i_t),0*muTE(i_t),X,Z,xTE_2(:,i_t)',zTE_2(:,i_t)',vtTE_2(1,:,i_t)',vnTE_2(1,:,i_t)',epSC,SC);
        end
                    
        u_TEd = zeros(nz,nx);
        w_TEd = zeros(nz,nx);
        for i = 1:nz
                vec = (i-1)*nx + 1:i*nx;

                u_TEd(i,:) = u_TEdt(vec) + u_TEdt_2(vec);
                w_TEd(i,:) = w_TEdt(vec) + w_TEdt_2(vec);
        end

        % Wake contribution
        [~,~,u_wdt,w_wdt] = DubSorV(muW(Ncyc*Nstep + 2 - i_t:end),0*muW(Ncyc*Nstep + 2 - i_t:end),X,Z,xw(Ncyc*Nstep + 2 - i_t:end),zw(Ncyc*Nstep + 2 - i_t:end),vtw(Ncyc*Nstep + 2 - i_t:end,:)',vnw(Ncyc*Nstep + 2 - i_t:end,:)',epSC,SC);

        u_wdt_2 = 0*u_wdt;
        w_wdt_2 = 0*w_wdt;
        
        if grd == 1
            [~,~,u_wdt_2,w_wdt_2] = DubSorV(muW(Ncyc*Nstep + 2 - i_t:end),0*muW(Ncyc*Nstep + 2 - i_t:end),X,Z,xw_2(Ncyc*Nstep + 2 - i_t:end),zw_2(Ncyc*Nstep + 2 - i_t:end),vtw_2(Ncyc*Nstep + 2 - i_t:end,:)',vnw_2(Ncyc*Nstep + 2 - i_t:end,:)',epSC,SC);
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
        set(gcf,'DefaultAxesfontsize',20,'DefaultAxesfontname','TimesNewRoman','DefaultAxesGridLineStyle','-.')

        hold on
        axis equal
        
        if grd == 1
           plot([min(xp(:,i_t))-10*c; min(xp(:,i_t))+10*c],[0 0],'-k','linewidth',4)
        end
        
        bluevec = [linspace(0,1,100)' linspace(0,1,100)' ones(100,1)];
        redvec = [ones(100,1) linspace(1,0,100)' linspace(1,0,100)'];
        vec = [bluevec; redvec(2:end,:)];
        colormap(vec)
        
        pcolor(Xstar,Zstar,-omega_y)
        quiver(Xf,Zf,u_p,w_p,'k')
        plot(xp(:,i_t),zp(:,i_t),'-k','linewidth',2)
        
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
        axis([-c/4 + x_b(i_t) 2*c + x_b(i_t) -c/10 c])
%         axis([-c/2 + x_b(i_t) 3*c + x_b(i_t) -3*c/4 + z_b(i_t) 3*c/4 + z_b(i_t)])
%         xlabel('$$x$$','interpreter','latex','fontsize',12,'fontname','TimesNewRoman')
%         ylabel('$$z$$','interpreter','latex','fontsize',12,'fontname','TimesNewRoman')
        
        VortMov(i_t-1) = getframe;
%         print('-dpng','-r300',['Grd_St_25_dc_16_Ac_17_',num2str(i_t)]);

%         print('-dpng','-r300',['2D_NearBody_',num2str(i_t)]);
    end

    
    
    
    %% Free-swimming calculations
    
        
    if free == 1 && i_t > 1
        if i_t > 2
            close(velh)
        end
        
        a_b(1,i_t) = Fx(1,i_t)/M;
        Q0(1,i_t + 1) = a_b(1,i_t)*delT + Q0(1,i_t);
        if i_t == 1
            x_b(1,i_t + 1) = a_b(1,i_t)*delT^2 + 2*x_b(1,i_t);
        else
            x_b(1,i_t + 1) = a_b(1,i_t)*delT^2 + 2*x_b(1,i_t) - x_b(1,i_t - 1); 
        end
        velh = figure('Position',[0 0 scrsz(3)/2 scrsz(4)/3]);
    
        hold on
%         subplot(2,1,1)
        plot(((1:i_t) - 1)*delT,-Q0(1,1:i_t)','-b','linewidth',2)
        xlabel('$$t$$, sec','interpreter','latex','fontname','TimesNewRoman','fontsize',12)
        ylabel('Velocity','interpreter','latex','fontname','TimesNewRoman','fontsize',12)
        axis([0 Ncyc/f 0 1.5*max(-Q0(1,1:i_t))])
        
%         subplot(2,1,2)
%         plot((1:Ncyc*Nstep+1)*delT,D_visc,'-r','linewidth',2)
%         xlabel('$$t$$, sec','interpreter','latex','fontname','TimesNewRoman','fontsize',12)
%         ylabel('Viscous Drag','interpreter','latex','fontname','TimesNewRoman','fontsize',12)
%         axis([0 Ncyc/f 0 1.5*max(D_visc)])
    else
        x_b(1,i_t + 1) = Q0(1,i_t)*delT + x_b(1,i_t);
    end
    z_b(1,i_t + 1) = Q0(2,i_t)*delT + z_b(1,i_t);
    
    
    
    close(wbar);
    delTime(i_t) = toc;
    EstT = delTime(i_t)*(Ncyc*Nstep+1 - i_t);
    EstThour = floor(EstT/3600);
    EstTmin = floor((EstT - EstThour*3600)/60);
    EstTsec = round(EstT - EstThour*3600 - EstTmin*60);
    wbar = waitbar(i_t/(Ncyc*Nstep+1),['Calculating...Time Remaining: ',num2str(EstThour),' hrs ',num2str(EstTmin),' mins ',num2str(EstTsec),' secs'],'Position',[(scrsz(3) - scrsz(4)/3) 0 scrsz(4)/3 1/16*scrsz(4)]);
end

close(wbar);

% End timer.
Time = sum(delTime);
Thour = floor(Time/3600);
Tmin = floor((Time - Thour*3600)/60);
Tsec = round(Time - Thour*3600 - Tmin*60);
sprintf(['Calculation Time: ',num2str(Thour),' hrs ',num2str(Tmin),' mins ',num2str(Tsec),' secs'],'FontName','TimesNewRoman')


% Time averaged forces
L_avg = mean(L(1,(Ncyc-1)*Nstep+1:end))
T_avg = mean(T(1,(Ncyc-1)*Nstep+1:end))
Pow_avg = mean(Pow(1,(Ncyc-1)*Nstep+1:end))
Cl_avg = mean(Cl(1,(Ncyc-1)*Nstep+1:end))
Ct_avg = mean(Ct(1,(Ncyc-1)*Nstep+1:end))
Cpow_avg = mean(Cpow(1,(Ncyc-1)*Nstep+1:end))
np = real(Ct_avg)/real(Cpow_avg)
Ec = mean(-Q0(1,(Ncyc-1)*Nstep+1:end))/Pow_avg

Cl_end = Cl(end);
Ct_end = Ct(end);

thickness = (max(zp) - min(zp))/c;


% save(['Pitch_St0',num2str(St*10),'_A_c0',num2str(2*A_TE/c*100),'_d_cInf_N210_Nstep80_Ep1e-5.mat'],'-v7.3')



% % Analytical solution.
% 
% [~,~,~,~,a_vdv,k,epsilon,theta1,theta2,theta,r1,r2] = VanDeVooren(c,200,tmax/2);
% % 
% % xp_a = [xpb;xpt(2:end)];
% % % zp = [zpb;zpt(2:end)];
% % zp_a = [zpb(1:end-1)-zpb(1);0;zpt(2:end)-zpt(end)];
% 
% Aa = cos((k-1)*theta1).*cos(k*theta2) + sin((k-1)*theta1).*sin(k*theta2);
% Ba = sin((k-1)*theta1).*cos(k*theta2) - cos((k-1)*theta1).*sin(k*theta2);
% D0 = a_vdv*(1 - k + k*epsilon);
% D1 = Aa.*(a_vdv*cos(theta) - D0) - Ba.*a_vdv.*sin(theta);
% D2 = Aa.*a_vdv.*sin(theta) + Ba.*(a_vdv.*cos(theta) - D0);
% 
% u_anltc = 2*Qinf*r2.^k./r1.^(k-1).*(sin(alpha) - sin(alpha - theta))./(D1.^2 + D2.^2).*(D1.*sin(theta) + D2.*cos(theta));
% w_anltc = -2*Qinf*r2.^k./r1.^(k-1).*(sin(alpha) - sin(alpha - theta))./(D1.^2 + D2.^2).*(D1.*cos(theta) - D2.*sin(theta));
% 
% Cp_anltc_t = 1 - (u_anltc.^2 + w_anltc.^2)/Qinf^2;
% 
% alpha_b = -alpha;
% [xpt,zpt,xpb,zpb,a_vdv,k,epsilon,theta1,theta2,theta,r1,r2] = VanDeVooren(c,200,tmax/2);
% 
% xp_a = [xpb;xpt(2:end)];
% % zp = [zpb;zpt(2:end)];
% zp_a = [zpb(1:end-1)-zpb(1);0;zpt(2:end)-zpt(end)];
% 
% Aa = cos((k-1)*theta1).*cos(k*theta2) + sin((k-1)*theta1).*sin(k*theta2);
% Ba = sin((k-1)*theta1).*cos(k*theta2) - cos((k-1)*theta1).*sin(k*theta2);
% D0 = a_vdv*(1 - k + k*epsilon);
% D1 = Aa.*(a_vdv*cos(theta) - D0) - Ba.*a_vdv.*sin(theta);
% D2 = Aa.*a_vdv.*sin(theta) + Ba.*(a_vdv.*cos(theta) - D0);
% 
% u_anltc = 2*Qinf*r2.^k./r1.^(k-1).*(sin(alpha_b) - sin(alpha_b - theta))./(D1.^2 + D2.^2).*(D1.*sin(theta) + D2.*cos(theta));
% w_anltc = -2*Qinf*r2.^k./r1.^(k-1).*(sin(alpha_b) - sin(alpha_b - theta))./(D1.^2 + D2.^2).*(D1.*cos(theta) - D2.*sin(theta));
% 
% Cp_anltc_b = 1 - (u_anltc.^2 + w_anltc.^2)/Qinf^2;


% % Plotting discretized airfoil with normals and lift vectors.
% figure
% hold on
% plot(xp(:,end),zp(:,end),'.-k')
% plot(xTE(:,end),zTE(:,end),'.-g')
% plot(xw,zw,'.-b')
% if grd == 1
%     plot([min(xp(:,i_t))-10*c; min(xp(:,i_t))+10*c],[0 0],'-k','linewidth',4)
% end
% axis equal
% xlabel('$$x$$','interpreter','latex','fontsize',12,'fontname','TimesNewRoman')
% ylabel('$$z$$','interpreter','latex','fontsize',12,'fontname','TimesNewRoman')
% axis([-c/2 + x_b(i_t) 4*c + x_b(i_t) -c/10 c])


% figure
% hold on
% plot(((1:i_t) - 1)*delT,x_b(1:end-1),'-b','linewidth',2)
% xlabel('$$t$$, sec','interpreter','latex','fontname','TimesNewRoman','fontsize',12)
% ylabel('Position','interpreter','latex','fontname','TimesNewRoman','fontsize',12)
      

% % Plotting pressure coefficient both numerical and analytical solutions.
% figure
% hold on
% % plot(xp_a(round(199/2):end),-Cp_anltc_t,'-k')
% % plot(xp_a(round(199/2):end),-Cp_anltc_b,'-k')
% plot((xc(1:Npanels/2,end) - x_b(end-1))/c,-Cp(1:Npanels/2,end),'sk')
% plot((xc(Npanels/2 + 1:end,end) - x_b(end-1))/c,-Cp(Npanels/2 + 1:end,end),'^k')
% 
% % plot(xp(2:round((Npanels + 1)/2)),-Cp(1:round((Npanels+1)/2)-1),'sk')
% % plot(xp(round((Npanels + 1)/2):end-1),-Cp(round((Npanels+1)/2)-1:end),'^k')
% xlabel('$$\frac{x}{c}$$','interpreter','latex','fontname','TimesNewRoman','fontsize',12)
% ylabel('$$-C_p$$','interpreter','latex','fontname','TimesNewRoman','fontsize',12)
% % axis([0 c -1 2])
% % axis([-1.6 -1.3 -1 2])
% legend('Bottom','Top','Location','NorthEast')
% 



%% Plotting coefficient of pressure.
% FontSizeAx = 24;
% FontSizeLb = 32;
% afFigurePosition = [1 1 20 15];
% ylabelpos = [-0.1 0.5];
% xlabelpos = [0.5 -0.15];
% axespos = [0.15 0.2 0.8 0.75];
% Panel = 186;
% 
% figure
% set(gcf, 'Units', 'centimeters','PaperPositionMode', 'auto','Position', afFigurePosition);
% set(gcf,'DefaultAxesFontSize',FontSizeAx,'DefaultAxesFontName','TimesNewRoman','DefaultAxesGridLineStyle','-.','DefaultAxesLineWidth',2,'DefaultAxesFontWeight','Normal')
% 
% t = ([1:Ncyc*Nstep+1] - 1)*delT;  
% delt = (t((Ncyc-1)*Nstep+1:end) - t((Ncyc-1)*Nstep+1))*f;
% Press = Cp(Panel,(Ncyc-1)*Nstep+1:end)*(1/2*rho*Qinf^2);
% Press_avg = mean(Press);
% Press_var = Press - Press_avg;
% 
% hold on
% %     plot(delt,Cp_s(Panel,(Ncyc-1)*Nstep+1:end)*(1/2*rho*Qinf^2),'--b','linewidth',3)
% %     plot(delt,Cp_us(Panel,(Ncyc-1)*Nstep+1:end)*(1/2*rho*Qinf^2),'--r','linewidth',3)
%     plot(delt,Press_var,'-k','linewidth',3)
% hold off
% 
% axis([0 1 1.2*min(Press_var) 1.2*max(Press_var)])
% set(gca, 'FontName', 'TimesNewRoman', 'FontSize', FontSizeAx)
% set(gca, 'Units', 'normalized', 'Position', axespos);
% xlabel('$$t/T$$','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 0,'Units', 'Normalize','Position', xlabelpos,'Interpreter', 'LaTeX');
% ylabel('$$P$$ (Pa)','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 90,'Units', 'Normalize','Position', ylabelpos,'Interpreter', 'LaTeX');
% 
% 
% 
% 
% 
% % Plotting of the forces acting on the body.
% figure
% set(gcf, 'Units', 'centimeters','PaperPositionMode', 'auto','Position', afFigurePosition);
% set(gcf,'DefaultAxesFontSize',FontSizeAx,'DefaultAxesFontName','TimesNewRoman','DefaultAxesGridLineStyle','-.','DefaultAxesLineWidth',2,'DefaultAxesFontWeight','Normal')
% 
% hold on
%     plot(delt,L((Ncyc-1)*Nstep+1:end),'-r','linewidth',3)
%     plot(delt,T((Ncyc-1)*Nstep+1:end),'-b','linewidth',3)
% hold off
% 
% % axis([0 1 1.2*min(Press_var) 1.2*max(Press_var)])
% set(gca, 'FontName', 'TimesNewRoman', 'FontSize', FontSizeAx)
% set(gca, 'Units', 'normalized', 'Position', axespos);
% xlabel('$$t/T$$','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 0,'Units', 'Normalize','Position', xlabelpos,'Interpreter', 'LaTeX');
% ylabel('$$F$$ (N)','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 90,'Units', 'Normalize','Position', ylabelpos,'Interpreter', 'LaTeX');
% 
% 
% % Calculating the moment from each panel about the leading edge.
% % Leading edge origin:
% x0 = kron(xp(Npanels/2 + 1,:),ones(Npanels,1));
% z0 = kron(zp(Npanels/2 + 1,:),ones(Npanels,1));
% r0(:,1,:) = x0;
% r0(:,2,:) = 0*z0;
% r0(:,3,:) = z0;
% r(:,1,:) = xc;
% r(:,2,:) = 0*zc;
% r(:,3,:) = zc;
% delF3(:,1,:) = delF(:,1,:);
% delF3(:,2,:) = 0*delF(:,1,:);
% delF3(:,3,:) = delF(:,2,:);
% 
% dM = cross(r-r0,delF3);
% dM = dM(:,2,:);
% Torque(1,:) = -sum(dM,1);
% theta_dot = -2*pi*f*alpha_max*cos(2*pi*f*t);
% Pow_T = abs(Torque.*theta_dot);
% Pow_T_avg = mean(Pow_T((Ncyc-1)*Nstep+1:end))
% 
% C_pow_T = Pow_T_avg/(1/2*rho*Qinf^3*c*b)
% 
% % Plotting of the forces acting on the body.
% FontSizeAx = 24;
% FontSizeLb = 32;
% afFigurePosition = [1 1 20 15];
% ylabelpos = [-0.125 0.5];
% xlabelpos = [0.5 -0.15];
% axespos = [0.175 0.2 0.8 0.75];
% 
% figure
% set(gcf, 'Units', 'centimeters','PaperPositionMode', 'auto','Position', afFigurePosition);
% set(gcf,'DefaultAxesFontSize',FontSizeAx,'DefaultAxesFontName','TimesNewRoman','DefaultAxesGridLineStyle','-.','DefaultAxesLineWidth',2,'DefaultAxesFontWeight','Normal')
% 
% hold on
%     plot(delt,Pow((Ncyc-1)*Nstep+1:end),'-r','linewidth',3)
%      plot([0 1],[Pow_avg Pow_avg],'-r','linewidth',1)
%     plot(delt,Pow_T((Ncyc-1)*Nstep+1:end),'-b','linewidth',3)
%     plot([0 1],[Pow_T_avg Pow_T_avg],'-b','linewidth',1)
% hold off
% 
% % axis([0 1 1.2*min(Press_var) 1.2*max(Press_var)])
% set(gca, 'FontName', 'TimesNewRoman', 'FontSize', FontSizeAx)
% set(gca, 'Units', 'normalized', 'Position', axespos);
% xlabel('$$t/T$$','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 0,'Units', 'Normalize','Position', xlabelpos,'Interpreter', 'LaTeX');
% ylabel('$$Pow$$ (W)','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 90,'Units', 'Normalize','Position', ylabelpos,'Interpreter', 'LaTeX');

% figure
% set(gcf, 'Units', 'centimeters','PaperPositionMode', 'auto','Position', afFigurePosition);
% set(gcf,'DefaultAxesFontSize',FontSizeAx,'DefaultAxesFontName','TimesNewRoman','DefaultAxesGridLineStyle','-.','DefaultAxesLineWidth',2,'DefaultAxesFontWeight','Normal')
% 
% a_TE = [(zp(end,(Ncyc-1)*Nstep+1:end) - d_c*c)/A_TE]';
% hold on
%     plot(delt,a_TE,'-k','linewidth',3)
% hold off
% 
% axis([0 1 1.2*min(a_TE) 1.2*max(a_TE)])
% set(gca, 'FontName', 'TimesNewRoman', 'FontSize', FontSizeAx)
% set(gca, 'Units', 'normalized', 'Position', axespos);
% xlabel('$$t/T$$','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 0,'Units', 'Normalize','Position', xlabelpos,'Interpreter', 'LaTeX');
% ylabel('$$z_{TE}/A_{TE}$$','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 90,'Units', 'Normalize','Position', ylabelpos,'Interpreter', 'LaTeX');


%% Calculating the flowfield around the wing.
% nx = 75;
% nz = 75;
% 
% U = Uinf*ones(nz,nx);
% W = Winf*ones(nz,nx);
% 
% %         xf = linspace(-c/2 + x_b(i_t),3*c + x_b(i_t),nx)';
% %         zf = linspace(-3*c/4 + z_b(i_t),3*c/4 + z_b(i_t),nz)';
% 
% % Flow field for ground effect calculations
% xf = linspace(-c/2 + x_b(i_t),5*c + x_b(i_t),nx)';
% zf = linspace(0,c,nz)';
% 
% [Xf,Zf] = meshgrid(xf,zf);
% 
% X = zeros(1,nx*nz);
% Z = zeros(1,nx*nz);
% 
% for i = 1:nz
%     vec = (i-1)*nx + 1:i*nx;
% 
%     X(1,vec) = Xf(i,:);
%     Z(1,vec) = Zf(i,:);
% end
% 
% % Calculating flow field
% 
% % Body contribution
% [u_st,w_st,u_bdt,w_bdt] = DubSorV(mu(:,i_t)',sigma(:,i_t)',X,Z,xp(:,i_t)',zp(:,i_t)',vt(:,:,i_t)',vn(:,:,i_t)',epSC,SC);
% 
% u_st_2 = 0*u_st;
% w_st_2 = 0*w_st;
% u_bdt_2 = 0*u_bdt;
% w_bdt_2 = 0*w_bdt;
% 
% if grd == 1
%     [u_st_2,w_st_2,u_bdt_2,w_bdt_2] = DubSorV(mu(:,i_t)',sigma(:,i_t)',X,Z,xp_2(:,i_t)',zp_2(:,i_t)',vt_2(:,:,i_t)',vn_2(:,:,i_t)',epSC,SC);
% end
% 
% u_s = zeros(nz,nx);
% w_s = zeros(nz,nx);
% u_bd = zeros(nz,nx);
% w_bd = zeros(nz,nx);
% 
% for i = 1:nz
%         vec = (i-1)*nx + 1:i*nx;
% 
%         u_s(i,:) = u_st(vec) +  u_st_2(vec);
%         w_s(i,:) = w_st(vec) + w_st_2(vec);
%         u_bd(i,:) = u_bdt(vec) + u_bdt_2(vec);
%         w_bd(i,:) = w_bdt(vec) + w_bdt_2(vec);
% end
% 
% 
% % TE contribution  
% [~,~,u_TEdt,w_TEdt] = DubSorV(muTE(i_t),0*muTE(i_t),X,Z,xTE(:,i_t)',zTE(:,i_t)',vtTE(1,:,i_t)',vnTE(1,:,i_t)',epSC,SC);
% 
% u_TEdt_2 = 0*u_TEdt;
% w_TEdt_2 = 0*w_TEdt;
% 
% if grd == 1
%     [~,~,u_TEdt_2,w_TEdt_2] = DubSorV(muTE(i_t),0*muTE(i_t),X,Z,xTE_2(:,i_t)',zTE_2(:,i_t)',vtTE_2(1,:,i_t)',vnTE_2(1,:,i_t)',epSC,SC);
% end
% 
% u_TEd = zeros(nz,nx);
% w_TEd = zeros(nz,nx);
% for i = 1:nz
%         vec = (i-1)*nx + 1:i*nx;
% 
%         u_TEd(i,:) = u_TEdt(vec) + u_TEdt_2(vec);
%         w_TEd(i,:) = w_TEdt(vec) + w_TEdt_2(vec);
% end
% 
% % Wake contribution
% [~,~,u_wdt,w_wdt] = DubSorV(muW(Ncyc*Nstep + 2 - i_t:end),0*muW(Ncyc*Nstep + 2 - i_t:end),X,Z,xw(Ncyc*Nstep + 2 - i_t:end),zw(Ncyc*Nstep + 2 - i_t:end),vtw(Ncyc*Nstep + 2 - i_t:end,:)',vnw(Ncyc*Nstep + 2 - i_t:end,:)',epSC,SC);
% 
% u_wdt_2 = 0*u_wdt;
% w_wdt_2 = 0*w_wdt;
% 
% if grd == 1
%     [~,~,u_wdt_2,w_wdt_2] = DubSorV(muW(Ncyc*Nstep + 2 - i_t:end),0*muW(Ncyc*Nstep + 2 - i_t:end),X,Z,xw_2(Ncyc*Nstep + 2 - i_t:end),zw_2(Ncyc*Nstep + 2 - i_t:end),vtw_2(Ncyc*Nstep + 2 - i_t:end,:)',vnw_2(Ncyc*Nstep + 2 - i_t:end,:)',epSC,SC);
% end
% 
% u_wd = zeros(nz,nx);
% w_wd = zeros(nz,nx);
% for i = 1:nz
%         vec = (i-1)*nx + 1:i*nx;
% 
%         u_wd(i,:) = u_wdt(vec) + u_wdt_2(vec);
%         w_wd(i,:) = w_wdt(vec) + w_wdt_2(vec);
% end
% 
% if LES == 1
%     vec = i_t:-1:2;
%     u_LEtdt = zeros(1,length(X),length(vec));
%     u_LEbdt = zeros(1,length(X),length(vec));
%     w_LEtdt = zeros(1,length(X),length(vec));
%     w_LEbdt = zeros(1,length(X),length(vec));
% 
%     for j = vec
%         % LE sheet top contribution
%         [~,~,u_LEtdt(1,:,j-1),w_LEtdt(1,:,j-1)] = DubSorV(muLEt(Ncyc*Nstep + 2 - j),0*muLEt(Ncyc*Nstep + 2 - j),X,Z,xpt_LES(:,Ncyc*Nstep + 2 - j)',zpt_LES(:,Ncyc*Nstep + 2 - j),vtLEt(Ncyc*Nstep + 2 - j,:)',vnLEt(Ncyc*Nstep + 2 - j,:)',epSC,SC);
% 
%         % LE sheet bot contribution
%         [~,~,u_LEbdt(1,:,j-1),w_LEbdt(1,:,j-1)] = DubSorV(muLEb(Ncyc*Nstep + 2 - j),0*muLEb(Ncyc*Nstep + 2 - j),X,Z,xpb_LES(:,Ncyc*Nstep + 2 - j)',zpb_LES(:,Ncyc*Nstep + 2 - j),vtLEb(Ncyc*Nstep + 2 - j,:)',vnLEb(Ncyc*Nstep + 2 - j,:)',epSC,SC);
%     end
% 
%     u_LEtdt = sum(u_LEtdt,3);
%     u_LEbdt = sum(u_LEbdt,3);
%     w_LEtdt = sum(w_LEtdt,3);
%     w_LEbdt = sum(w_LEbdt,3);
% 
%     u_LEtd = zeros(nz,nx);
%     w_LEtd = zeros(nz,nx);
%     u_LEbd = zeros(nz,nx);
%     w_LEbd = zeros(nz,nx);
%     for i = 1:nz
%             vec = (i-1)*nx + 1:i*nx;
% 
%             u_LEtd(i,:) = u_LEtdt(vec);
%             w_LEtd(i,:) = w_LEtdt(vec);
%             u_LEbd(i,:) = u_LEbdt(vec);
%             w_LEbd(i,:) = w_LEbdt(vec);
%     end
% end
% 
% if LES == 1
%     u_d = u_bd + u_TEd + u_wd + u_LEtd + u_LEbd;
%     w_d = w_bd + w_TEd + w_wd + w_LEtd + w_LEbd;
% else
%     u_d = u_bd + u_TEd + u_wd;
%     w_d = w_bd + w_TEd + w_wd;
% end
% 
% 
% u_p = u_d + u_s;
% w_p = w_d + w_s;
% 
% U = U + u_p;
% W = W + w_p;
% 
% omega_y = (u_p(3:end,2:end-1) - u_p(1:end-2,2:end-1))./(Zf(3:end,2:end-1) - Zf(1:end-2,2:end-1)) - (w_p(2:end-1,3:end) - w_p(2:end-1,1:end-2))./(Xf(2:end-1,3:end) - Xf(2:end-1,1:end-2));
% Xstar = (Xf(2:end-1,3:end) + Xf(2:end-1,1:end-2))/2;
% Zstar = (Zf(3:end,2:end-1) + Zf(1:end-2,2:end-1))/2;
% 
% 
% 
% % % Plotting discretized airfoil with normals and lift vectors.
% % figure
% % set(gcf,'DefaultAxesfontsize',20,'DefaultAxesfontname','TimesNewRoman','DefaultAxesGridLineStyle','-.')
% % 
% % hold on
% % plot(xp(:,end),zp(:,end),'.-k','linewidth',2)
% % plot(xTE(:,end),zTE(:,end),'.-g')
% % plot(xw,zw,'.-g')
% % plot(Xc,Zc,'xb')
% % if grd == 1
% %     plot([min(xp(:,i_t))-10*c; min(xp(:,i_t))+10*c],[0 0],'-k','linewidth',4)
% % end
% % streamline(Xf,Zf,U,W,Xf(:,1),Zf(:,1))
% % % quiver(Xf,Zf,U,W)
% % axis equal
% % axis([-c/2 + x_b(end-1) 2*c + x_b(end-1) -c/10 c])
% % % axis([-c/2 + x_b(end-1) 2*c + x_b(end-1) -c/2 + z_b(end-1) c/2 + z_b(end-1)])
% % % title('2D Solution','interpreter','latex','fontsize',14,'fontname','TimesNewRoman')
% % % xlabel('$$x$$','interpreter','latex','fontsize',12,'fontname','TimesNewRoman')
% % % ylabel('$$z$$','interpreter','latex','fontsize',12,'fontname','TimesNewRoman')
% % % print('-dpng','-r300',['2D_Steady_Streamline_alpha_0']);
% 
% % Plotting discretized airfoil with normals and lift vectors.
% figure
% hold on
% 
% c_range = max(max(abs(omega_y)));
% plot(xp(:,end),zp(:,end),'.-k')
% plot(xTE(:,end),zTE(:,end),'.-g')
% if grd == 1
%     plot([min(xp(:,i_t))-10*c; min(xp(:,i_t))+10*c],[0 0],'-k','linewidth',4)
% end
% % plot(xw,zw,'.-g')
% % plot(Xc,Zc,'xb')
% % streamline(Xf,Zf,U,W,Xf(:,1),Zf(:,1))
% 
% bluevec = [linspace(0,1,100)' linspace(0,1,100)' ones(100,1)];
% redvec = [ones(100,1) linspace(1,0,100)' linspace(1,0,100)'];
% vec = [bluevec; redvec(2:end,:)];
% colormap(vec)
%     
% pcolor(Xstar,Zstar,-omega_y)
% quiver(Xf,Zf,U,W,'k')
% shading interp
% caxis(0.7*[-c_range c_range])
% axis equal
% axis([-c/2 + x_b(end-1) 5*c + x_b(end-1) -c/10 c])

% axis([-c/2 + x_b(end-1) 5*c + x_b(end-1) -c + z_b(end-1) c + z_b(end-1)])
% title('2D Solution','interpreter','latex','fontsize',14,'fontname','TimesNewRoman')
% xlabel('$$x$$','interpreter','latex','fontsize',12,'fontname','TimesNewRoman')
% ylabel('$$z$$','interpreter','latex','fontsize',12,'fontname','TimesNewRoman')
% print('-dpng','-r300',['2D_Wake_VF']);




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

