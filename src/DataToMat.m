function out = DataToMat(f,A_c,d_c,Qinf)

savefilename = ['_PitchGrd_f',num2str(f),...
            '_A_c',num2str(A_c),...
            '_d_c',num2str(d_c)];
        
%% Loading data 
load(['FlowfieldData/Figure9/Parameters',savefilename,'.mat']);
Data = load(['FlowfieldData/Figure9/Data',savefilename,'.txt']);
PanelProp = load(['FlowfieldData/Figure9/PanelProp',savefilename,'.txt']);
WakeProp = load(['FlowfieldData/Figure9/WakeProp',savefilename,'.txt']);

%% Defining variables
Qinf
alpha = 0;

% Free-stream velocity
Uinf = Qinf*cos(alpha);
Winf = Qinf*sin(alpha);

wakeInd = Data(:,1)';
Q0(1,:) = Data(:,2)';
Q0(2,:) = Data(:,3)';
x_b = Data(:,4)';
z_b = Data(:,5)';
a_b = Data(:,6)';
Fx = Data(:,7)';
Fz = Data(:,8)';
Pow = Data(:,9)';
L = Data(:,10)';
T = Data(:,11)';
D_visc = Data(:,12)';
Gamma = Data(:,13)';
Cf = Data(:,14)';
Cl = Data(:,15)';
Ct = Data(:,16)';
Cpow = Data(:,17)';
Fx_s = Data(:,18)';
Fx_us = Data(:,19)';
Fz_s = Data(:,20)';
Fz_us = Data(:,21)';
Pow_s = Data(:,22)';
Pow_us = Data(:,23)';
L_s = Data(:,24)';
L_us = Data(:,25)';
T_s = Data(:,26)';
T_us = Data(:,27)';
Cl_s = Data(:,28)';
Cl_us = Data(:,29)';
Ct_s = Data(:,30)';
Ct_us = Data(:,31)';
Cpow_s = Data(:,32)';
Cpow_us = Data(:,33)';

clear Data

xp_0 = reshape(PanelProp(:,1),[Npanels (Ncyc*Nstep + 1)]);
xp_0 = [xp_0;xp_0(1,:)];
zp_0 = reshape(PanelProp(:,2),[Npanels (Ncyc*Nstep + 1)]);
zp_0 = [zp_0;zp_0(1,:)];
xp = reshape(PanelProp(:,3),[Npanels (Ncyc*Nstep + 1)]);
xp = [xp;xp(1,:)];
zp = reshape(PanelProp(:,4),[Npanels (Ncyc*Nstep + 1)]);
zp = [zp;zp(1,:)];
Vc(:,:,1)  = reshape(PanelProp(:,5),[Npanels (Ncyc*Nstep + 1)]);
Vc(:,:,2) = reshape(PanelProp(:,6),[Npanels (Ncyc*Nstep + 1)]);
xc = reshape(PanelProp(:,7),[Npanels (Ncyc*Nstep + 1)]);
zc = reshape(PanelProp(:,8),[Npanels (Ncyc*Nstep + 1)]);
vt(:,:,1) = reshape(PanelProp(:,9),[Npanels (Ncyc*Nstep + 1)]);
vt(:,:,2) = reshape(PanelProp(:,10),[Npanels (Ncyc*Nstep + 1)]);
vn(:,:,1) = reshape(PanelProp(:,11),[Npanels (Ncyc*Nstep + 1)]);
vn(:,:,2) = reshape(PanelProp(:,12),[Npanels (Ncyc*Nstep + 1)]);
dL = reshape(PanelProp(:,13),[Npanels (Ncyc*Nstep + 1)]);
Xc = reshape(PanelProp(:,14),[Npanels (Ncyc*Nstep + 1)]);
Zc = reshape(PanelProp(:,15),[Npanels (Ncyc*Nstep + 1)]);
sigma = reshape(PanelProp(:,16),[Npanels (Ncyc*Nstep + 1)]);
mu = reshape(PanelProp(:,17),[Npanels (Ncyc*Nstep + 1)]);
Qp = reshape(PanelProp(:,18),[Npanels (Ncyc*Nstep + 1)]);
Qt = reshape(PanelProp(:,19),[Npanels (Ncyc*Nstep + 1)]);
Cp_s = reshape(PanelProp(:,20),[Npanels (Ncyc*Nstep + 1)]);
Cp_us = reshape(PanelProp(:,21),[Npanels (Ncyc*Nstep + 1)]);
Cp = reshape(PanelProp(:,22),[Npanels (Ncyc*Nstep + 1)]);
dFshear = reshape(PanelProp(:,23),[Npanels (Ncyc*Nstep + 1)]);
      
clear PanelProp

xwtemp = reshape(WakeProp(:,1),[(Nlump*Nstep + 2) (Ncyc*Nstep + 1)]);
xw = xwtemp(1:end-1,:);
xl = xwtemp(end,:);
zwtemp = reshape(WakeProp(:,2),[(Nlump*Nstep + 2) (Ncyc*Nstep + 1)]);
zw = zwtemp(1:end-1,:);
zl = zwtemp(end,:);
vtwtemp_1 = reshape(WakeProp(:,3),[(Nlump*Nstep + 2) (Ncyc*Nstep + 1)]);
vtwtemp_2 = reshape(WakeProp(:,4),[(Nlump*Nstep + 2) (Ncyc*Nstep + 1)]);
vtTE(1,:,1) = vtwtemp_1(1,:);
vtTE(1,:,2) = vtwtemp_2(1,:);
vtw(:,:,1) = vtwtemp_1(2:end-1,:);
vtw(:,:,2) = vtwtemp_2(2:end-1,:);
vtlw(1,:,1) = vtwtemp_1(end,:);
vtlw(1,:,2) = vtwtemp_2(end,:);
vnwtemp_1 = reshape(WakeProp(:,5),[(Nlump*Nstep + 2) (Ncyc*Nstep + 1)]);
vnwtemp_2 = reshape(WakeProp(:,6),[(Nlump*Nstep + 2) (Ncyc*Nstep + 1)]);
vnTE(1,:,1) = vnwtemp_1(1,:);
vnTE(1,:,2) = vnwtemp_2(1,:);
vnw(:,:,1) = vnwtemp_1(2:end-1,:);
vnw(:,:,2) = vnwtemp_2(2:end-1,:);
vnlw(1,:,1) = vnwtemp_1(end,:);
vnlw(1,:,2) = vnwtemp_2(end,:);
muWtemp = reshape(WakeProp(:,7),[(Nlump*Nstep + 2) (Ncyc*Nstep + 1)]);
muTE(1,:) = muWtemp(1,:);
muW = muWtemp(2:end-1,:);
muLump(1,:) = muWtemp(end,:);
GammaWtemp = reshape(WakeProp(:,8),[(Nlump*Nstep + 2) (Ncyc*Nstep + 1)]);
GammaW = GammaWtemp(1:end-1,:);
GammaLump = GammaWtemp(end,:);

clear WakeProp xwtemp zwtemp vtwtemp_1 vtwtemp_2 vnwtemp_1 vnwtemp_2 muWtemp GammaWtemp

save(['FlowfieldData/Figure9/Processed_',savefilename,'.mat'],'-v7.3')   

out = 1