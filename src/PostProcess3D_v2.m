clear
close all
clc

% profile -memory on
% setpref('profiler','showJitLines',1);

% Parameters for identifying data file.
% Atip_b = exp(-0.1342*(0.4568)^2)*sin(0.4568);
fvec = [0.25 0.5 0.75 1 1.25 1.5]';
 M_star = 1/3;
    A_bp = 1;
    A_c = 0.25;
%     f = 1.5;
% 
% savefilename = ['_Pitch_Mstar',num2str(M_star),...
%         '_A_bp',num2str(A_bp),...
%         '_A_c',num2str(A_c),...
%         '_f',num2str(f)];
    
% load(['Pow_Vel_Data.mat']);

for j = 1:length(fvec)
    f = fvec(j);
    delT = 1/f/100;

   
    
% % savefilename = ['_MantaBot_Free_A_b',num2str(round((Atip_b)*1000))];
   savefilename = ['_2_Pitch_Mstar',num2str(M_star),...
        '_A_bp',num2str(A_bp),...
        '_A_c',num2str(A_c),...
        '_f',num2str(f)];

%% Loading data 
load(['FreeV_Visc/Parameters',savefilename,'.mat']);


% load(['MantabotData/Parameters',savefilename,'.mat']);
% Data = load(['MantabotData/Data',savefilename,'.txt']);
% PanelProp = load(['MantabotData/PanelProp',savefilename,'.txt']);
% WakeProp = load(['MantabotData/WakeProp',savefilename,'.txt']);

%% Defining variables

% Npanels = Ncpanels*Nspanels;

% Free-stream velocity
Uinf = Qinf*cos(alpha);
Winf = Qinf*sin(alpha);

Data = load(['FreeV_Visc/Data',savefilename,'.txt']);

U(:,j) = -Data(:,2)';
Pow(:,j) = Data(:,9)';
t(:,j) = [((1:Ncyc*Nstep+1) - 1)*delT]';

% clear Data
% PanelProp = load(['Scratch/PanelProp',savefilename,'.txt']);
% 
% pc = reshape(PanelProp(:,1:12)',[12 2*Npanels (Ncyc*Nstep + 1)]);
% pcI = reshape(PanelProp(:,13:24)',[12 2*Npanels (Ncyc*Nstep + 1)]);
% Vc  = reshape(PanelProp(:,25:27)',[3 2*Npanels (Ncyc*Nstep + 1)]);
% cpts  = reshape(PanelProp(:,28:30)',[3 2*Npanels (Ncyc*Nstep + 1)]);
% cptsI  = reshape(PanelProp(:,31:33)',[3 2*Npanels (Ncyc*Nstep + 1)]);
% vc  = reshape(PanelProp(:,34:36)',[3 2*Npanels (Ncyc*Nstep + 1)]);
% vcI  = reshape(PanelProp(:,37:39)',[3 2*Npanels (Ncyc*Nstep + 1)]);
% vn  = reshape(PanelProp(:,40:42)',[3 2*Npanels (Ncyc*Nstep + 1)]);
% vnI  = reshape(PanelProp(:,43:45)',[3 2*Npanels (Ncyc*Nstep + 1)]);
% vt  = reshape(PanelProp(:,46:48)',[3 2*Npanels (Ncyc*Nstep + 1)]);
% vtI  = reshape(PanelProp(:,49:51)',[3 2*Npanels (Ncyc*Nstep + 1)]);
% S  = reshape(PanelProp(:,52)',[1 2*Npanels (Ncyc*Nstep + 1)]);
% SI  = reshape(PanelProp(:,53)',[1 2*Npanels (Ncyc*Nstep + 1)]);
% dl  = reshape(PanelProp(:,54)',[1 2*Npanels (Ncyc*Nstep + 1)]);
% dlI  = reshape(PanelProp(:,55)',[1 2*Npanels (Ncyc*Nstep + 1)]);
% dm  = reshape(PanelProp(:,56)',[1 2*Npanels (Ncyc*Nstep + 1)]);
% dmI  = reshape(PanelProp(:,57)',[1 2*Npanels (Ncyc*Nstep + 1)]);
% dm_hat  = reshape(PanelProp(:,58:60)',[3 2*Npanels (Ncyc*Nstep + 1)]);
% tp  = reshape(PanelProp(:,61)',[1 2*Npanels (Ncyc*Nstep + 1)]);
% cptsP  = reshape(PanelProp(:,62:64)',[3 2*Npanels (Ncyc*Nstep + 1)]);
% sig  = reshape(PanelProp(:,65)',[1 2*Npanels (Ncyc*Nstep + 1)]);
% mu  = reshape(PanelProp(:,66)',[1 2*Npanels (Ncyc*Nstep + 1)]);
% Qt  = reshape(PanelProp(:,67:69)',[3 2*Npanels (Ncyc*Nstep + 1)]);
% Cp_s  = reshape(PanelProp(:,70)',[1 2*Npanels (Ncyc*Nstep + 1)]);
% Cp_us  = reshape(PanelProp(:,71)',[1 2*Npanels (Ncyc*Nstep + 1)]);
% Cp  = reshape(PanelProp(:,72)',[1 2*Npanels (Ncyc*Nstep + 1)]);
% dFshear  = reshape(PanelProp(:,73)',[1 2*Npanels (Ncyc*Nstep + 1)]);
% dtau  = reshape(PanelProp(:,74)',[1 2*Npanels (Ncyc*Nstep + 1)]);
%      
% clear PanelProp
% WakeProp = load(['Scratch/WakeProp',savefilename,'.txt']);
% 
% 
% pcwtemp = reshape(WakeProp(:,1:12)',[12 Nspanels*(Nlump*Nstep + 2) (Ncyc*Nstep + 1)]);
% pcTE = pcwtemp(:,1:Nspanels,:);
% pcw = pcwtemp(:,Nspanels+1:Nspanels*(Nlump*Nstep + 1),:);
% pcl = pcwtemp(:,Nspanels*(Nlump*Nstep + 1)+1:end,:);
% cptstemp = reshape(WakeProp(:,13:15)',[3 Nspanels*(Nlump*Nstep + 2) (Ncyc*Nstep + 1)]);
% cptsTE = cptstemp(:,1:Nspanels,:);
% cptsw = cptstemp(:,Nspanels+1:Nspanels*(Nlump*Nstep + 1),:);
% cptslw = cptstemp(:,Nspanels*(Nlump*Nstep + 1)+1:end,:);
% dltemp = reshape(WakeProp(:,16)',[1 Nspanels*(Nlump*Nstep + 2) (Ncyc*Nstep + 1)]);
% dlTE = dltemp(1,1:Nspanels,:);
% dlw = dltemp(1,Nspanels+1:Nspanels*(Nlump*Nstep + 1),:);
% dllw = dltemp(1,Nspanels*(Nlump*Nstep + 1)+1:end,:);
% dmtemp = reshape(WakeProp(:,17)',[1 Nspanels*(Nlump*Nstep + 2) (Ncyc*Nstep + 1)]);
% dmTE = dmtemp(1,1:Nspanels,:);
% dmw = dmtemp(1,Nspanels+1:Nspanels*(Nlump*Nstep + 1),:);
% dmlw = dmtemp(1,Nspanels*(Nlump*Nstep + 1)+1:end,:);
% vcTEtemp = reshape(WakeProp(:,18:20)',[3 Nspanels*(Nlump*Nstep + 2) (Ncyc*Nstep + 1)]);
% vcTE = vcTEtemp(:,1:Nspanels,:);
% vcw = vcTEtemp(:,Nspanels+1:Nspanels*(Nlump*Nstep + 1),:);
% vclw = vcTEtemp(:,Nspanels*(Nlump*Nstep + 1)+1:end,:);
% vnTEtemp = reshape(WakeProp(:,21:23)',[3 Nspanels*(Nlump*Nstep + 2) (Ncyc*Nstep + 1)]);
% vnTE = vnTEtemp(:,1:Nspanels,:);
% vnw = vnTEtemp(:,Nspanels+1:Nspanels*(Nlump*Nstep + 1),:);
% vnlw = vnTEtemp(:,Nspanels*(Nlump*Nstep + 1)+1:end,:);
% vtTEtemp = reshape(WakeProp(:,24:26)',[3 Nspanels*(Nlump*Nstep + 2) (Ncyc*Nstep + 1)]);
% vtTE = vtTEtemp(:,1:Nspanels,:);
% vtw = vtTEtemp(:,Nspanels+1:Nspanels*(Nlump*Nstep + 1),:);
% vtlw = vtTEtemp(:,Nspanels*(Nlump*Nstep + 1)+1:end,:);
% STEtemp = reshape(WakeProp(:,27)',[1 Nspanels*(Nlump*Nstep + 2) (Ncyc*Nstep + 1)]);
% STE = STEtemp(1,1:Nspanels,:);
% Sw = STEtemp(1,Nspanels+1:Nspanels*(Nlump*Nstep + 1),:);
% Slw = STEtemp(1,Nspanels*(Nlump*Nstep + 1)+1:end,:);
% muTEtemp = reshape(WakeProp(:,28)',[1 Nspanels*(Nlump*Nstep + 2) (Ncyc*Nstep + 1)]);
% muTE = muTEtemp(1,1:Nspanels,:);
% muW = muTEtemp(1,Nspanels+1:Nspanels*(Nlump*Nstep + 1),:);
% muLump = muTEtemp(1,Nspanels*(Nlump*Nstep + 1)+1:end,:);
% 
% clear WakeProp pcwtemp cptstemp dltemp dmtemp vcTEtemp vnTEtemp vtTEtemp STEtemp muTEtemp 
% WakePropI = load(['Scratch/WakePropI',savefilename,'.txt']);
% 
% pcwtemp = reshape(WakePropI(:,1:12)',[12 Nspanels*(Nlump*Nstep + 2) (Ncyc*Nstep + 1)]);
% pcTEI = pcwtemp(:,1:Nspanels,:);
% pcwI = pcwtemp(:,Nspanels+1:Nspanels*(Nlump*Nstep + 1),:);
% pclI = pcwtemp(:,Nspanels*(Nlump*Nstep + 1)+1:end,:);
% cptstemp = reshape(WakePropI(:,13:15)',[3 Nspanels*(Nlump*Nstep + 2) (Ncyc*Nstep + 1)]);
% cptsTEI = cptstemp(:,1:Nspanels,:);
% cptswI = cptstemp(:,Nspanels+1:Nspanels*(Nlump*Nstep + 1),:);
% cptslI = cptstemp(:,Nspanels*(Nlump*Nstep + 1)+1:end,:);
% dltemp = reshape(WakePropI(:,16)',[1 Nspanels*(Nlump*Nstep + 2) (Ncyc*Nstep + 1)]);
% dlTEI = dltemp(1,1:Nspanels,:);
% dlwI = dltemp(1,Nspanels+1:Nspanels*(Nlump*Nstep + 1),:);
% dllI = dltemp(1,Nspanels*(Nlump*Nstep + 1)+1:end,:);
% dmtemp = reshape(WakePropI(:,17)',[1 Nspanels*(Nlump*Nstep + 2) (Ncyc*Nstep + 1)]);
% dmTEI = dmtemp(1,1:Nspanels,:);
% dmwI = dmtemp(1,Nspanels+1:Nspanels*(Nlump*Nstep + 1),:);
% dmlI = dmtemp(1,Nspanels*(Nlump*Nstep + 1)+1:end,:);
% vcTEtemp = reshape(WakePropI(:,18:20)',[3 Nspanels*(Nlump*Nstep + 2) (Ncyc*Nstep + 1)]);
% vcTEI = vcTEtemp(:,1:Nspanels,:);
% vcwI = vcTEtemp(:,Nspanels+1:Nspanels*(Nlump*Nstep + 1),:);
% vclI = vcTEtemp(:,Nspanels*(Nlump*Nstep + 1)+1:end,:);
% vnTEtemp = reshape(WakePropI(:,21:23)',[3 Nspanels*(Nlump*Nstep + 2) (Ncyc*Nstep + 1)]);
% vnTEI = vnTEtemp(:,1:Nspanels,:);
% vnwI = vnTEtemp(:,Nspanels+1:Nspanels*(Nlump*Nstep + 1),:);
% vnlI = vnTEtemp(:,Nspanels*(Nlump*Nstep + 1)+1:end,:);
% vtTEtemp = reshape(WakePropI(:,24:26)',[3 Nspanels*(Nlump*Nstep + 2) (Ncyc*Nstep + 1)]);
% vtTEI = vtTEtemp(:,1:Nspanels,:);
% vtwI = vtTEtemp(:,Nspanels+1:Nspanels*(Nlump*Nstep + 1),:);
% vtlI = vtTEtemp(:,Nspanels*(Nlump*Nstep + 1)+1:end,:);
% STEtemp = reshape(WakePropI(:,27)',[1 Nspanels*(Nlump*Nstep + 2) (Ncyc*Nstep + 1)]);
% STEI = STEtemp(1,1:Nspanels,:);
% SwI = STEtemp(1,Nspanels+1:Nspanels*(Nlump*Nstep + 1),:);
% SlI = STEtemp(1,Nspanels*(Nlump*Nstep + 1)+1:end,:);
% 
% clear WakePropI pcwtemp cptstemp dltemp dmtemp vcTEtemp vnTEtemp vtTEtemp STEtemp muTEtemp 
% 
end 
save(['Pow_Vel_Data.mat'],'-v7.3','U','Pow','t','f')   
% end


% load(['Scratch/MR_free_f_vary_2.mat']);

%% Velocity Plotting
i_t = Ncyc*Nstep+1

figure
hold on
    plot(((1:i_t) - 1)*delT,-Q0(1,1:i_t)','-b','linewidth',2)
hold off

xlabel('$$t$$, sec','interpreter','latex','fontname','TimesNewRoman','fontsize',12)
ylabel('Velocity','interpreter','latex','fontname','TimesNewRoman','fontsize',12)
axis([0 min([5*Nstep*(i_t)*delT Ncyc/f]) 0 1.5*max(-Q0(1,1:i_t))])

U_avg = mean(-Q0(1,(Ncyc-1)*Nstep+1:end))

%% Body/Wake Plotting
Flowfield = 1;
BodyWake = 1;
VelFig = 0;

C1 = [-1  1  0  0;...
       0 -1  1  0;...
       0  0 -1  1;...
       1  0  0 -1];  
   
ind = 1;

for i_t = Ncyc*Nstep+1:Ncyc*Nstep+1
        
    if BodyWake == 1
        wakeInd = min([(i_t-1) Nlump*Nstep]);
        

%         pause(0.25)
        % Closing previous figure 
        if i_t > 2
%             close(fighand)
%             close(fighand(i_t-2))
        end

        scrsz = get(0,'ScreenSize');
        FontSizeAx = 24;
        FontSizeLb = 32;
        afFigurePosition = [0 15 20 10];
        axespos = [0.15 0.2 0.84 0.8];
        ylabelpos = [-0.1 0.5];
        xlabelpos = [0.5 -0.1];
        
        
        % Plotting the three-dimensional fin geometry.
        wvec = 1:wakeInd*Nspanels;
        fighand = figure('Position',[0 1/2*scrsz(4) scrsz(3)/2 1/2*scrsz(4)]);
%         FigHand(i_t) = figure('Position',[0 1/2*scrsz(4) scrsz(3)/2 1/2*scrsz(4)]);
        set(gcf,'DefaultAxesfontsize',20,'DefaultAxesfontname','TimesNewRoman','DefaultAxesGridLineStyle','-.')

        hold on
        PlotPanelsMatrix(i_t,pc(:,:,i_t),pcI(:,:,i_t),Vc(:,:,i_t),pcTE(:,:,i_t),pcTEI(:,:,i_t),pcw(:,wvec,i_t),pcwI(:,wvec,i_t),pcl(:,:,i_t),pclI(:,:,i_t),C1,cptsP(:,:,i_t),cptsI(:,:,i_t),cptsTE(:,:,i_t),cptsTEI(:,:,i_t),cptsw(:,wvec,i_t),cptswI(:,wvec,i_t),vc(:,:,i_t),vcI(:,:,i_t),vn(:,:,i_t),vnI(:,:,i_t),vt(:,:,i_t),vtI(:,:,i_t),vcTE(:,:,i_t),vcTEI(:,:,i_t),vnTE(:,:,i_t),vnTEI(:,:,i_t),vtTE(:,:,i_t),vtTEI(:,:,i_t),vcw(:,wvec,i_t),vcwI(:,wvec,i_t),vnw(:,wvec,i_t),vnwI(:,wvec,i_t),vtw(:,wvec,i_t),vtwI(:,wvec,i_t),Npanels,Nspanels,Nlump,Nstep,Cp(1,:,i_t),dtau(1,:,i_t),SurfaceColor)
        if grd == 1
            PlotPanelsMatrix(i_t,pc2,pcI2,Vc,pcTE2,pcTEI2,pcw2(:,wvec),pcwI2(:,wvec),pcl2(:,:,1),pclI2(:,:,1),C1,cpts2,cptsI2,cptsTE2,cptsTEI2,cptsw2(:,wvec),cptswI2(:,wvec),vc2,vcI2,vn2,vnI2,vt2,vtI2,vcTE2,vcTEI2,vnTE2,vnTEI2,vtTE2,vtTEI2,vcw2(:,wvec),vcwI2(:,wvec),vnw2(:,wvec),vnwI2(:,wvec),vtw2(:,wvec),vtwI2(:,wvec),Npanels,Nspanels,Nlump,Nstep,Cp(1,:),dtau(1,:),SurfaceColor)
            PlotGrd(Ncyc,Nstep,delT,Q0,c_b,b,c_r,C1);
        end
    %     xlabel('x, m','Interpreter','Latex','FontName','TimesNewRoman','FontSize',12)
    %     ylabel('y, m','Interpreter','Latex','FontName','TimesNewRoman','FontSize',12)
    %     zlabel('z, m','Interpreter','Latex','FontName','TimesNewRoman','FontSize',12)
    %     title('Fin Geometry','Interpreter','Latex','FontName','TimesNewRoman','FontSize',14)
        axis equal
    %     axis([(Ncyc*Nstep + 2)*delT*Q0(1) -  c_r/2 1.05*c_r -1.05*b 1.05*b Ncyc*Nstep*delT*Q0(3)-3*Atip - D 3*Atip + D])
    %     axis([x_b(2) - 1 2.3 -3 3 -2 2])
    %     axis([x_b(i_t) - 1/4 1/4 -1/5 1/5 -1/5 1/5])
        view(3)
    %     print('-depsc','-r300',['Free_eps/2D_Pitch_Free_',num2str(i_t)]);
        if ind < 10 
            print('-dpng','-r300',['MantaFreeData/PanelWake/MantaFree_f',num2str(round(f*1000)),'_PanelWake_00',num2str(ind)]);
        elseif ind < 100
            print('-dpng','-r300',['MantaFreeData/PanelWake/MantaFree_f',num2str(round(f*1000)),'_PanelWake_0',num2str(ind)]);
        else
            print('-dpng','-r300',['MantaFreeData/PanelWake/MantaFree_f',num2str(round(f*1000)),'_PanelWake_',num2str(ind)]);
        end
%         pause(0.001);
    end
    
    if VelFig == 1
        if i_t > 2
            close(velh)
        end

        FontSizeAx = 24;
        FontSizeLb = 32;
        afFigurePosition = [0 0 20 10];
        axespos = [0.2 0.2 0.76 0.78];
        ylabelpos = [-0.175 0.5];
        xlabelpos = [0.5 -0.1];
        
        velh = figure;
        set(gcf, 'Units', 'centimeters','PaperPositionMode', 'auto','Position', afFigurePosition);
        set(gcf,'DefaultAxesFontSize',FontSizeAx,'DefaultAxesFontName','TimesNewRoman','DefaultAxesGridLineStyle','-.','DefaultAxesLineWidth',2,'DefaultAxesFontWeight','Normal')

        hold on
%             plot([0 Ncyc/f],U_swim*[1 1],'-k','linewidth',1)
%             plot(N_a/f*[1 1],[0 1.5*max(-Q0(1,1:i_t))],'--k','linewidth',1)
            plot(((1:i_t) - 1)*delT,-Q0(1,1:i_t)','-b','linewidth',2)
        hold off

        set(gca, 'FontName', 'TimesNewRoman', 'FontSize', FontSizeAx)
        set(gca, 'Units', 'normalized', 'Position', axespos);
        xlabel('t/T','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 0,'Units', 'Normalize','Position',xlabelpos,'Interpreter', 'LaTeX');
        ylabel('U(t), m/s','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 90,'Units', 'Normalize','Position',ylabelpos,'Interpreter', 'LaTeX');
        axis([0 Ncyc/f 0 1.5*max(-Q0(1,1:i_t))])
%         if ind < 10 
%             print('-depsc','-r600',['FreeSwim/Vel_2/2D_Pitch_Free_00',num2str(ind)]);
%         elseif ind < 100
%             print('-depsc','-r600',['FreeSwim/Vel_2/2D_Pitch_Free_0',num2str(ind)]);
%         else
%             print('-depsc','-r600',['FreeSwim/Vel_2/2D_Pitch_Free_',num2str(ind)]);
%         end
%         pause(0.001);
    end
    ind = ind + 1;
end


%% Flowfield calculation
figInd = 1;
% epSC = 1e-3;
ep = c_r*1e-2;
epSC = ep;

% Calculating the flowfield around the wing.
nx = 20;
ny = 20;
nz = 20;
% Ut = zeros(nx,ny,nz,Nstep);
% Wt = zeros(nx,ny,nz,Nstep);

if Flowfield == 1
    for i_t = Ncyc*Nstep+1:Ncyc*Nstep+1
%             for i_t = (Ncyc-1)*Nstep+1:Ncyc*Nstep+1
        if i_t > Nlump*Nstep + 1
            lump = 1;
        else
            lump = [];
        end
        wakeInd = min([(i_t-1) Nlump*Nstep]);
        
%         U = Uinf*ones(nx,ny,nz);
%         W = Winf*ones(nx,ny,nz);

        % Flow field for ground effect calculations
        xf = linspace(-(c_r - c_b)/2 + x_b(1,i_t),6*c_r + x_b(1,i_t),nx)';
        yf = linspace(-3/2*b*(b_b + 1),3/2*b*(b_b + 1),ny)';
        zf = linspace(z_b(1,i_t) - 1.5*c_r,z_b(1,i_t) + 1.5*c_r,nz)';

        [Xf,Yf,Zf] = meshgrid(xf,yf,zf);
        
        X = reshape(Xf,1,nx*ny*nz);
        Y = reshape(Yf,1,nx*ny*nz);
        Z = reshape(Zf,1,nx*ny*nz);

        [u_p,v_p,w_p] = VelocityField(X,Y,Z,pc(:,:,i_t),pcTE(:,:,i_t),[pcw(:,1:Nspanels*wakeInd,i_t) kron(lump,pcl(:,:,i_t))],pcI(:,:,i_t),pcTEI(:,:,i_t),[pcwI(:,1:Nspanels*wakeInd,i_t) kron(lump,pclI(:,:,i_t))],cpts(:,:,i_t),cptsTE(:,:,i_t),[cptsw(:,1:Nspanels*wakeInd,i_t) kron(lump,cptslw(:,:,i_t))],cptsI(:,:,i_t),cptsTEI(:,:,i_t),[cptswI(:,1:Nspanels*wakeInd,i_t) kron(lump,cptslI(:,:,i_t))],vc(:,:,i_t),vcTE(:,:,i_t),[vcw(:,1:Nspanels*wakeInd,i_t) kron(lump,vclw(:,:,i_t))],vcI(:,:,i_t),vcTEI(:,:,i_t),[vcwI(:,1:Nspanels*wakeInd,i_t) kron(lump,vclI(:,:,i_t))],vn(:,:,i_t),vnTE(:,:,i_t),[vnw(:,1:Nspanels*wakeInd,i_t) kron(lump,vnlw(:,:,i_t))],vnI(:,:,i_t),vnTEI(:,:,i_t),[vnwI(:,1:Nspanels*wakeInd,i_t) kron(lump,vnlI(:,:,i_t))],vt(:,:,i_t),vtTE(:,:,i_t),[vtw(:,1:Nspanels*wakeInd,i_t) kron(lump,vtlw(:,:,i_t))],vtI(:,:,i_t),vtTEI(:,:,i_t),[vtwI(:,1:Nspanels*wakeInd,i_t) kron(lump,vtlI(:,:,i_t))],mu(1,:,i_t)',muTE(1,:,i_t)',[muW(1,1:Nspanels*wakeInd,i_t)'; kron(lump,muLump(1,:,i_t)')],sig(1,:,i_t)',Nspanels,dl(1,:,i_t)',dm(1,:,i_t)',S(1,:,i_t)',dlI(1,:,i_t)',dmI(1,:,i_t)',SI(1,:,i_t)',dlTE(1,:,i_t)',dmTE(1,:,i_t)',STE(1,:,i_t)',dlTEI(1,:,i_t)',dmTEI(1,:,i_t)',STEI(1,:,i_t)',[dlw(1,1:Nspanels*wakeInd,i_t) kron(lump,dllw(1,:,i_t))],[dmw(1,1:Nspanels*wakeInd,i_t) kron(lump,dmlw(1,:,i_t))],[Sw(1,1:Nspanels*wakeInd,i_t) kron(lump,Slw(1,:,i_t))],[dlwI(1,1:Nspanels*wakeInd,i_t) kron(lump,dllI(1,:,i_t))],[dmwI(1,1:Nspanels*wakeInd,i_t) kron(lump,dmlI(1,:,i_t))],[SwI(1,1:Nspanels*wakeInd,i_t) kron(lump,SlI(1,:,i_t))],ep,epSC,SCw);

        u_p = reshape(u_p,nx,ny,nz);
        v_p = reshape(v_p,nx,ny,nz);
        w_p = reshape(w_p,nx,ny,nz);
        
        Ut(:,:,:,figInd) = u_p;
        Vt(:,:,:,figInd) = v_p;
        Wt(:,:,:,figInd) = w_p;
        
        xfi = linspace(-(c_r - c_b)/2 + x_b(1,i_t),6*c_r + x_b(1,i_t),4*nx)';
        yfi = linspace(-3/2*b*(b_b + 1),3/2*b*(b_b + 1),4*ny)';
        zfi = linspace(z_b(1,i_t) - 1.5*c_r,z_b(1,i_t) + 1.5*c_r,4*nz)';
        [Xfi,Yfi,Zfi] = meshgrid(xfi,yfi,zfi);
        
        u_pi = interp3(Xf,Yf,Zf,u_p,Xfi,Yfi,Zfi);
        v_pi = interp3(Xf,Yf,Zf,v_p,Xfi,Yfi,Zfi);
        w_pi = interp3(Xf,Yf,Zf,w_p,Xfi,Yfi,Zfi);

        [omega_x,omega_y,omega_z,Lambda2,Xstar,Ystar,Zstar] = Lambda2Crit3D(Xfi,Yfi,Zfi,u_pi,v_pi,w_pi);

        if figInd > 1
            close(flowfig)
        end

%%
        % Plotting airfoil with LE vortex sheets and the TE vortex
        % sheet.
        flowfig = figure;
        FontSizeAx = 24;
        afFigurePosition = [15 7 25 15];
        wvec = 1:wakeInd*Nspanels;

        set(flowfig, 'Units', 'centimeters','PaperPositionMode', 'auto','Position', afFigurePosition);
        set(flowfig,'DefaultAxesFontSize',FontSizeAx,'DefaultAxesFontName','TimesNewRoman','DefaultAxesGridLineStyle','-.','DefaultAxesLineWidth',2,'DefaultAxesFontWeight','Normal')

        


        bluevec = [linspace(0,1,100)' linspace(0,1,100)' ones(100,1)];
        redvec = [ones(100,1) linspace(1,0,100)' linspace(1,0,100)'];
        vec = [bluevec; redvec(2:end,:)];
        colormap(vec)

%         pcolor(Xstar,Zstar,-Lambda2.*omega_y)
%         pcolor(Xstar,Zstar,-omega_y)
%         quiver3(Xf,Yf,Zf,u_p,v_p,w_p,'k')
    
        hold on
        
            surf_patch = patch(isosurface(Xstar,Ystar,Zstar,Lambda2,-0.25));
            isonormals(Xstar,Ystar,Zstar,Lambda2,surf_patch)
            set(surf_patch,'FaceColor','b','FaceAlpha',0.55,'FaceLighting','Phong','EdgeColor', 'none')

            surf_patch = patch(isosurface(Xstar,Ystar,Zstar,Lambda2,-0.15));
            isonormals(Xstar,Ystar,Zstar,Lambda2,surf_patch)
            set(surf_patch,'FaceColor','b','FaceAlpha',0.35,'FaceLighting','Phong','EdgeColor', 'none')

            surf_patch = patch(isosurface(Xstar,Ystar,Zstar,Lambda2,-0.05));
            isonormals(Xstar,Ystar,Zstar,Lambda2,surf_patch)
            set(surf_patch,'FaceColor','r','FaceAlpha',0.35,'FaceLighting','Phong','EdgeColor', 'none')

            PlotBodyPanelsMatrix(pc(:,:,i_t),pcI(:,:,i_t),Npanels,Cp(1,:,i_t),dtau(1,:,i_t),0)         
        hold off
        
        axis equal
        view(3)
        axis([-c_r/2 + x_b(1,i_t) 6*c_r + x_b(1,i_t) -3/2*b*(b_b + 1) 3/2*b*(b_b + 1) z_b(1,i_t) - 1.5*c_r z_b(1,i_t) + 1.5*c_r])

        %%%%%%%%%%%%%%%%%%%%%%%%
        
%         flowfig2 = figure;
%         set(flowfig2, 'Units', 'centimeters','PaperPositionMode', 'auto','Position', afFigurePosition);
%         set(flowfig2,'DefaultAxesFontSize',FontSizeAx,'DefaultAxesFontName','TimesNewRoman','DefaultAxesGridLineStyle','-.','DefaultAxesLineWidth',2,'DefaultAxesFontWeight','Normal')
% 
%         
%         Vmag = sqrt(u_pi.^2 + v_pi.^2 + w_pi.^2);
%         
%         hold on
%        
%             surf_patch = patch(isosurface(Xfi,Yfi,Zfi,Vmag,0.7*max(max(max(Vmag)))));
%             isonormals(Xfi,Yfi,Zfi,Vmag,surf_patch)
%             set(surf_patch,'FaceColor',[1 0 0],'FaceAlpha',0.35,'FaceLighting','Phong','EdgeColor', 'none')
% 
%             surf_patch = patch(isosurface(Xfi,Yfi,Zfi,Vmag,0.6*max(max(max(Vmag)))));
%             isonormals(Xfi,Yfi,Zfi,Vmag,surf_patch)
%             set(surf_patch,'FaceColor',[0.75 0 0.25],'FaceAlpha',0.35,'FaceLighting','Phong','EdgeColor', 'none')
% 
%             surf_patch = patch(isosurface(Xfi,Yfi,Zfi,Vmag,0.5*max(max(max(Vmag)))));
%             isonormals(Xfi,Yfi,Zfi,Vmag,surf_patch)
%             set(surf_patch,'FaceColor',[0 0 1],'FaceAlpha',0.35,'FaceLighting','Phong','EdgeColor', 'none')
%   
%             PlotBodyPanelsMatrix(pc(:,:,i_t),pcI(:,:,i_t),Npanels,Cp(1,:,i_t),dtau(1,:,i_t),1) 
%         hold off
%         %         camlight; 
% %         lighting phong
%         axis equal
%         view(3)
%         axis([-c_r/2 + x_b(1,i_t) 6*c_r + x_b(1,i_t) -3/2*b*(b_b + 1) 3/2*b*(b_b + 1) z_b(1,i_t) - 1.5*c_r z_b(1,i_t) + 1.5*c_r])
% 
% %         colorbar
% %         caxis(6*[-1 1])
% %         colorbar
% %         axis([-c/2 + x_b(i_t) 3*c + x_b(i_t) -3*c/4 + z_b(i_t) 3*c/4 + z_b(i_t)])
% %         xlabel('$$x$$','interpreter','latex','fontsize',12,'fontname','TimesNewRoman')
% %         ylabel('$$z$$','interpreter','latex','fontsize',12,'fontname','TimesNewRoman')
% 
%         VortMov(i_t-1) = getframe;
%         print('-dpng','-r300',['Grd_St_25_dc_16_Ac_17_',num2str(i_t)]);
% 
        if ind < 10 
            print('-dpng','-r300',['MantaFreeData/VortWake/MantaFree_f',num2str(round(f*1000)),'_VortWake_00',num2str(ind)]);
        elseif ind < 100
            print('-dpng','-r300',['MantaFreeData/VortWake/MantaFree_f',num2str(round(f*1000)),'_VortWake_0',num2str(ind)]);
        else
            print('-dpng','-r300',['MantaFreeData/VortWake/MantaFree_f',num2str(round(f*1000)),'_VortWake_',num2str(ind)]);
        end
       figInd = figInd + 1;    
    end   
end
    
 
save(['MantaFreeData/Flowfield',savefilename,'.mat'],'-v7.3','Ut','Vt','Wt','Lambda2','Xstar','Ystar','Zstar')    

% end

% profile viewer

%% Plotting time-averaged flowfield
% figure;
% FontSizeAx = 24;
% afFigurePosition = [15 7 25 15];
% 
% set(gcf, 'Units', 'centimeters','PaperPositionMode', 'auto','Position', afFigurePosition);
% set(gcf,'DefaultAxesFontSize',FontSizeAx,'DefaultAxesFontName','TimesNewRoman','DefaultAxesGridLineStyle','-.','DefaultAxesLineWidth',2,'DefaultAxesFontWeight','Normal')
% 
% hold on
% axis equal
% % % % 
% % % % % Flow field for ground effect calculations
% % % % xf = linspace(-c/2 + x_b(1,end),6*c + x_b(1,end),nx)';
% % % % zf = linspace(z_b(1,end) - 1.5*c,z_b(1,end) + 1.5*c,nz)';
% % % % 
% % % % [Xf,Zf] = meshgrid(xf,zf);
% % % % 
% % % % % if grd == 1
% % % % %    plot([min(xp)-10*c; min(xp)+10*c],[0 0],'-k','linewidth',1)
% % % % % end
% % % % % 
% % % % % bluevec = [linspace(0,1,100)' linspace(0,1,100)' ones(100,1)];
% % % % % redvec = [ones(100,1) linspace(1,0,100)' linspace(1,0,100)'];
% % % % % vec = [bluevec; redvec(2:end,:)];
% % % % % colormap(vec)
% % % % 
% % % % FFmean(:,:,1) = Xf;
% % % % FFmean(:,:,2) = Zf;
% % % % FFmean(:,:,3) = mean(Ut,3);
% % % % FFmean(:,:,4) = mean(Wt,3);
% % % % 
% % % % % quiver(Xf,Zf,mean(Ut-Uinf,3),mean(Wt,3),'k')
% % % % % % quiver(Xf,Zf,Uinf*ones(nz,nx),mean(Wt(:,:,1:151)-Winf,3),'k')
% % % % % plot(xp_0 + x_b(end) ,zp_0 + d_c*c,'-k','linewidth',2)
% % % % % hold off
% % % % 
% % % % %         if grd == 1
% % % % %             plot([min(xp(:,i_t))-10*c; min(xp(:,i_t))+10*c],[0 0],'-k','linewidth',4)
% % % % %         end
% % % % 
% % % % % shading interp
% % % % 
% % % % 
% % % % 
% % % % %         colorbar
% % % % % caxis([0.5 3])
% % % % % colorbar
% % % % % axis([-c/2 + x_b(1,end) 6*c + x_b(1,end) z_b(1,end) - 1.5*c z_b(1,end) + 1.5*c])
% % % % % print('-depsc','-r600',['FixedV/A_c_',num2str(A_c),'/TA_FF_St_',num2str(St*1000)]);
% % % % 
% % % % save(['FixedV/A_c_',num2str(A_c),'/FF_St_',num2str(St),'.mat'],'Ut','Wt','FFmean','f','Uinf','c','xp_0','zp_0','xp','zp','x_b','-v7.3')
% % % % % save(['FixedV/A_c_',num2str(A_c),'/Xp.mat'],'xp_0','zp_0','x_b','-v7.3')
% % % % 
% % % % %         axis([-c/2 + x_b(i_t) 3*c + x_b(i_t) -3*c/4 + z_b(i_t) 3*c/4 + z_b(i_t)])
% % % % %         xlabel('$$x$$','interpreter','latex','fontsize',12,'fontname','TimesNewRoman')
% % % % %         ylabel('$$z$$','interpreter','latex','fontsize',12,'fontname','TimesNewRoman')
% % % % 
% % % % %         VortMov(i_t-1) = getframe;
% % % % %         print('-dpng','-r300',['Grd_St_25_dc_16_Ac_17_',num2str(i_t)]);
% % % % 
% % % % %             print('-dpng','-r300',['GrdData/Vorticity/GrdEffect_St',num2str(St*1000),'_d_c',num2str(d_c*100),'_A_c',num2str(A_c),'_',num2str(figInd),'.png']);



%%%%%%%%%%%%%%%% TIME AVG
% Cp_avg = mean(Cp(1,:,(Ncyc-1)*Nstep + 1:end),3)    
%%%%%%%%%%%%%%%%

%%%%%%%%%%% TRANSFORM
%  L = Fz*cos(alpha) + Fx*sin(alpha);
%  T = -(-Fz*sin(alpha) + Fx*cos(alpha));
%%%%%%%%%%%








% %% Plotting the forces acting on the body
% FontSizeAx = 24;
% FontSizeLb = 32;
% afFigurePosition = [0 2 30 20];
% axespos = [0.175 0.15 0.775 0.8];
% ylabelpos = [-0.125 0.5];
% xlabelpos = [0.5 -0.1];
% 
% t = ([1:Ncyc*Nstep+1] - 1)*delT;  
% delt = (t((Ncyc-2)*Nstep+1:end) - t((Ncyc-2)*Nstep+1))*f;
% span = 0.05;
% 
% % Lift
% figure
% set(gcf, 'Units', 'centimeters','PaperPositionMode', 'auto','Position', afFigurePosition);
% set(gcf,'DefaultAxesFontSize',FontSizeAx,'DefaultAxesFontName','TimesNewRoman','DefaultAxesGridLineStyle','-.','DefaultAxesLineWidth',2,'DefaultAxesFontWeight','Normal')
% 
% hold on
%     plot(delt,Cl_s((Ncyc-2)*Nstep+1:end),'-b','linewidth',3)
%     plot(delt,smooth(Cl_us((Ncyc-2)*Nstep+1:end),span,'loess'),'-r','linewidth',3)
%     plot(delt,smooth(Cl((Ncyc-2)*Nstep+1:end),span,'loess'),'-k','linewidth',3)
% 
%     plot([min(delt) max(delt)],mean(Cl((Ncyc-2)*Nstep+1:end))*[1 1],':k','linewidth',2)
%     plot([min(delt) max(delt)],[0 0],'-k','linewidth',2)
% hold off
% 
% axis([0 2 1.2*min(Cl) 1.2*max(Cl)])
% set(gca, 'FontName', 'TimesNewRoman', 'FontSize', FontSizeAx)
% set(gca, 'Units', 'normalized', 'Position', axespos);
% legend('Steady','Unsteady','Total','Location','NorthWest')
% xlabel('$$t/T$$','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 0,'Units', 'Normalize','Position', xlabelpos,'Interpreter', 'LaTeX');
% ylabel('$$C_l$$','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 0,'Units', 'Normalize','Position', ylabelpos,'Interpreter', 'LaTeX');
% 
% % Thrust
% figure
% set(gcf, 'Units', 'centimeters','PaperPositionMode', 'auto','Position', afFigurePosition);
% set(gcf,'DefaultAxesFontSize',FontSizeAx,'DefaultAxesFontName','TimesNewRoman','DefaultAxesGridLineStyle','-.','DefaultAxesLineWidth',2,'DefaultAxesFontWeight','Normal')
% 
% hold on
%     plot(delt,Ct_s((Ncyc-2)*Nstep+1:end),'-b','linewidth',3)
%     plot(delt,smooth(Ct_us((Ncyc-2)*Nstep+1:end),span,'loess'),'-r','linewidth',3)
%     plot(delt,smooth(Ct((Ncyc-2)*Nstep+1:end),span,'loess'),'-k','linewidth',3)
%     
%     plot([min(delt) max(delt)],mean(Ct((Ncyc-2)*Nstep+1:end))*[1 1],':k','linewidth',2)
%     plot([min(delt) max(delt)],[0 0],'-k','linewidth',2)
% hold off
% 
% axis([0 2 1.2*min(Ct) 1.2*max(Ct)])
% set(gca, 'FontName', 'TimesNewRoman', 'FontSize', FontSizeAx)
% set(gca, 'Units', 'normalized', 'Position', axespos);
% legend('Steady','Unsteady','Total','Location','SouthEast')
% xlabel('$$t/T$$','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 0,'Units', 'Normalize','Position', xlabelpos,'Interpreter', 'LaTeX');
% ylabel('$$C_t$$','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 0,'Units', 'Normalize','Position', ylabelpos,'Interpreter', 'LaTeX');
% 
% % Circulation
% figure
% set(gcf, 'Units', 'centimeters','PaperPositionMode', 'auto','Position', afFigurePosition);
% set(gcf,'DefaultAxesFontSize',FontSizeAx,'DefaultAxesFontName','TimesNewRoman','DefaultAxesGridLineStyle','-.','DefaultAxesLineWidth',2,'DefaultAxesFontWeight','Normal')
% 
% hold on
%     plot(delt,(-2/c/Uinf)*Gamma((Ncyc-2)*Nstep+1:end),'-k','linewidth',3)
%     
%     plot([min(delt) max(delt)],mean((-2/c/Uinf)*Gamma((Ncyc-2)*Nstep+1:end))*[1 1],':k','linewidth',2)
%     plot([min(delt) max(delt)],mean(Cl_s((Ncyc-2)*Nstep+1:end))*[1 1],'--b','linewidth',2)
%     plot([min(delt) max(delt)],mean(Cl_us((Ncyc-2)*Nstep+1:end))*[1 1],'--g','linewidth',2)
%     plot([min(delt) max(delt)],[0 0],'-k','linewidth',2)
% hold off
% 
% axis([0 2 1.2*min((-2/c/Uinf)*Gamma) 1.2*max((-2/c/Uinf)*Gamma)])
% set(gca, 'FontName', 'TimesNewRoman', 'FontSize', FontSizeAx)
% set(gca, 'Units', 'normalized', 'Position', axespos);
% legend('Circulatory force','Time-avg circulatory force','Time-avg steady force','Time-avg unsteady force','Location','SouthEast')
% xlabel('$$t/T$$','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 0,'Units', 'Normalize','Position', xlabelpos,'Interpreter', 'LaTeX');
% ylabel('$$-\frac{2\Gamma}{cU_{\infty}}$$','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 0,'Units', 'Normalize','Position', ylabelpos,'Interpreter', 'LaTeX');
% 
% 
% % Pressure coefficient
% figure
% set(gcf, 'Units', 'centimeters','PaperPositionMode', 'auto','Position', afFigurePosition);
% set(gcf,'DefaultAxesFontSize',FontSizeAx,'DefaultAxesFontName','TimesNewRoman','DefaultAxesGridLineStyle','-.','DefaultAxesLineWidth',2,'DefaultAxesFontWeight','Normal')
% 
% xc_0 = (xp_0(1:end-1,1)/2 + xp_0(2:end,1)/2);
% hold on
%     plot(xc_0(1:Npanels/2)/c,mean(Cp(1:Npanels/2,(Ncyc-2)*Nstep+1:end),2),'ob','linewidth',3,'markersize',6)
%     plot(xc_0(Npanels/2+1:end)/c,mean(Cp(Npanels/2+1:end,(Ncyc-2)*Nstep+1:end),2),'sk','linewidth',3,'markersize',6)
% hold off
% 
% axis([0 1 1.2*min(mean(Cp(:,(Ncyc-2)*Nstep+1:end),2)) 1.2*max(mean(Cp(:,(Ncyc-2)*Nstep+1:end),2))])
% set(gca, 'FontName', 'TimesNewRoman', 'FontSize', FontSizeAx)
% set(gca, 'Units', 'normalized', 'Position', axespos);
% legend('Bottom','Top','Location','South')
% xlabel('$$x/c$$','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 0,'Units', 'Normalize','Position', xlabelpos,'Interpreter', 'LaTeX');
% ylabel('$$\bar{C_p}$$','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 0,'Units', 'Normalize','Position', ylabelpos,'Interpreter', 'LaTeX');
% 
% % end
% 
% 
% 
% 
% 









%% Airfoil plotting
% FontSizeAx = 24;
% FontSizeLb = 32;
% afFigurePosition = [0 2 30 20];
% axespos = [0.15 0.2 0.84 0.8];
% ylabelpos = [-0.12 0.5];
% xlabelpos = [0.5 -0.25];


% % Plotting airfoil with LE vortex sheets and the TE vortex sheet.
% fighand = figure;
% set(gcf, 'Units', 'centimeters','PaperPositionMode', 'auto','Position', afFigurePosition);
% set(gcf,'DefaultAxesFontSize',FontSizeAx,'DefaultAxesFontName','TimesNewRoman','DefaultAxesGridLineStyle','-.','DefaultAxesLineWidth',2,'DefaultAxesFontWeight','Normal')
% 
% hold on
% axis equal
% 
% plot(Xc(:,end),Zc(:,end),'xr','linewidth',2)
% plot(Xc(1,end),Zc(1,end),'xb','linewidth',2)
% plot(Xc(end,end),Zc(end,end),'xg','linewidth',2)
% plot(xp(:,end),zp(:,end),'.-k','linewidth',2)
% plot(xTE(:,end),zTE(:,end),'.-b','linewidth',2)
% % quiver(xc(:,end),zc(:,end),vn(:,1,end),vn(:,2,end),'b','linewidth',2)
% val = 1/10;
% 
% if grd == 1
% %         plot(xp_2(:,i_t),zp_2(:,i_t),'-k','linewidth',2)
% %         plot(xTE_2(:,i_t),zTE_2(:,i_t),'.-b','linewidth',2)
%     plot([min(xp(:,end))-10*c; min(xp(:,end))+10*c],[0 0],'-k','linewidth',4)
% %         if t > 0
% %             plot(xw_2(Ncyc*Nstep + 2 - i_t:end),zw_2(Ncyc*Nstep + 2 - i_t:end),'.-b','linewidth',2,'markersize',14)
% %         end
% end
% 
% if t > 0
%     if LES == 1
%         plot(xpt_LES(1,Ncyc*Nstep + 2 - i_t),zpt_LES(1,Ncyc*Nstep + 2 - i_t),'og','linewidth',2)
%         plot(xpb_LES(1,Ncyc*Nstep + 2 - i_t),zpb_LES(1,Ncyc*Nstep + 2 - i_t),'or','linewidth',2)
% 
%         for j = i_t:-1:2
%             if muLEt(Ncyc*Nstep + 2 - j) == 0
%             else
%                 plot(xpt_LES(:,Ncyc*Nstep + 2 - j),zpt_LES(:,Ncyc*Nstep + 2 - j),'.-g','linewidth',2,'markersize',14)
% %                 quiver(xpt_LES(1,Ncyc*Nstep + 1 - j),zpt_LES(1,Ncyc*Nstep + 1 - j),val*vtLEt(Ncyc*Nstep + 1 - j,1),val*vtLEt(Ncyc*Nstep + 1 - j,2),'b')
% %                 quiver(xpt_LES(1,Ncyc*Nstep + 1 - j),zpt_LES(1,Ncyc*Nstep + 1 - j),val*vnLEt(Ncyc*Nstep + 1 - j,1),val*vnLEt(Ncyc*Nstep + 1 - j,2),'k')
% 
%             end
%         end
% 
%         for j = i_t:-1:2
%             if muLEb(Ncyc*Nstep + 2 - j) == 0  
%             else
%                 plot(xpb_LES(:,Ncyc*Nstep + 2 - j),zpb_LES(:,Ncyc*Nstep + 2 - j),'.-r','linewidth',2,'markersize',14)
% %                 quiver(xpb_LES(1,Ncyc*Nstep + 1 - j),zpb_LES(1,Ncyc*Nstep + 1 - j),val*vtLEb(Ncyc*Nstep + 1 - j,1),val*vtLEb(Ncyc*Nstep + 1 - j,2),'b')
% %                 quiver(xpb_LES(1,Ncyc*Nstep + 1 - j),zpb_LES(1,Ncyc*Nstep + 1 - j),val*vnLEb(Ncyc*Nstep + 1 - j,1),val*vnLEb(Ncyc*Nstep + 1 - j,2),'k')
% 
%             end
%         end
%     end
% %         plot(xp(Stagpt,i_t),zp(Stagpt,i_t),'xk','linewidth',2)
% 
%     bluevec = [linspace(0,1,100)' linspace(0,1,100)' ones(100,1)];
%     redvec = [ones(100,1) linspace(1,0,100)' linspace(1,0,100)'];
%     vec = [bluevec; redvec(2:end,:)];
%     vec = vec(end:-1:1,:);
% 
%     WakeCirc = GammaW'/max(GammaW)/(1/2);
%     WakeCirc(WakeCirc > 1) = 1;  
%     WakeCirc(WakeCirc < -1) = -1;
%     WakeCirc = WakeCirc + 1;
%     WakeCirc = round(WakeCirc*100) + 1;
%     WakeCirc(WakeCirc > 199) = 199;
% 
%     for i_w = 1:Ncyc*Nstep+1
%         plot(xw(i_w),zw(i_w),'.','color',vec(WakeCirc(i_w),:),'linewidth',2,'markersize',14)
%     end
% end
% 
% axis equal
% %     axis([-c/2 + x_b(i_t) 2*c + x_b(i_t) -1/2*c + z_b(i_t) 1/2*c + z_b(i_t)])
% axis([-c/2 + x_b(end-1) 7*c + x_b(end-1) z_b(end-1) - 1*c z_b(end-1) + 1*c])
% 
% set(gca, 'FontName', 'TimesNewRoman', 'FontSize', FontSizeAx)
% set(gca, 'Units', 'normalized', 'Position', axespos);
% xlabel('$$x$$ (m)','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 0,'Units', 'Normalize','Position',xlabelpos,'Interpreter', 'LaTeX');
% ylabel('$$z$$ (m)','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 90,'Units', 'Normalize','Position',ylabelpos,'Interpreter', 'LaTeX');



%% Analytical solution.

% [~,~,~,~,a_vdv,k,epsilon,theta1,theta2,theta,r1,r2] = VanDeVooren(c,600,tmax/2);
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
% [xpt,zpt,xpb,zpb,a_vdv,k,epsilon,theta1,theta2,theta,r1,r2] = VanDeVooren(c,600,tmax/2);
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
% 
% 
% plot((xc(1:Npanels/2,end) - x_b(end-1))/c,-Cp(1:Npanels/2,end-round(Nstep/4)),'sk')
% plot((xc(Npanels/2 + 1:end,end) - x_b(end-1))/c,-Cp(Npanels/2 + 1:end,end-round(Nstep/4)),'^k')
% 
% % plot(xp_a(round(599/2):end)/c,-Cp_anltc_t,'-k')
% % plot(xp_a(round(599/2):end)/c,-Cp_anltc_b,'-k')
% 
% % plot(xp(2:round((Npanels + 1)/2)),-Cp(1:round((Npanels+1)/2)-1),'sk')
% % plot(xp(round((Npanels + 1)/2):end-1),-Cp(round((Npanels+1)/2)-1:end),'^k')
% xlabel('$$\frac{x}{c}$$','interpreter','latex','fontname','TimesNewRoman','fontsize',12)
% ylabel('$$-C_p$$','interpreter','latex','fontname','TimesNewRoman','fontsize',12)
% % axis([0 c -1 2])
% % axis([-1.6 -1.3 -1 2])
% legend('Bottom','Top','Location','NorthEast')
% 
% 


%% Plotting coefficient of pressure.
% FontSizeAx = 24;
% FontSizeLb = 32;
% afFigurePosition = [1 1 20 15];
% ylabelpos = [-0.1 0.5];
% xlabelpos = [0.5 -0.15];
% axespos = [0.15 0.2 0.8 0.75];
% Panelt = 189;
% Panelb = 21;
% 
% t = ([1:Ncyc*Nstep+1] - 1)*delT;  
% delt = (t((Ncyc-1)*Nstep+1:end) - t((Ncyc-1)*Nstep+1))*f;
% 
% Press = (Cp(Panelt,(Ncyc-1)*Nstep+1:end) + 1)*(1/2*rho*Qinf^2);
% Press_avg = mean(Press);
% Press_var = Press - Press_avg;
% 
% Press_max_t = (max(Cp(Npanels/2 + 1:end,(Ncyc-1)*Nstep+1:end),[],2))*(1/2*rho*Qinf^2);
% Press_max = Press_max_t - mean(Press_max_t);
% 
% figure
% set(gcf, 'Units', 'centimeters','PaperPositionMode', 'auto','Position', afFigurePosition);
% set(gcf,'DefaultAxesFontSize',FontSizeAx,'DefaultAxesFontName','TimesNewRoman','DefaultAxesGridLineStyle','-.','DefaultAxesLineWidth',2,'DefaultAxesFontWeight','Normal')
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
% figure
% set(gcf, 'Units', 'centimeters','PaperPositionMode', 'auto','Position', afFigurePosition);
% set(gcf,'DefaultAxesFontSize',FontSizeAx,'DefaultAxesFontName','TimesNewRoman','DefaultAxesGridLineStyle','-.','DefaultAxesLineWidth',2,'DefaultAxesFontWeight','Normal')
% 
% hold on
% %     plot(delt,Cp_s(Panel,(Ncyc-1)*Nstep+1:end)*(1/2*rho*Qinf^2),'--b','linewidth',3)
% %     plot(delt,Cp_us(Panel,(Ncyc-1)*Nstep+1:end)*(1/2*rho*Qinf^2),'--r','linewidth',3)
%     plot(xc(Npanels/2+1:end,1)/c,Press_max,'-k','linewidth',3)
% hold off
% 
% axis([0 1 1.2*min(Press_max) 1.2*max(Press_max)])
% set(gca, 'FontName', 'TimesNewRoman', 'FontSize', FontSizeAx)
% set(gca, 'Units', 'normalized', 'Position', axespos);
% xlabel('$$x/c$$','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 0,'Units', 'Normalize','Position', xlabelpos,'Interpreter', 'LaTeX');
% ylabel('$$P$$ (Pa)','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 90,'Units', 'Normalize','Position', ylabelpos,'Interpreter', 'LaTeX');

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

% %% Intantaneous Power
% 
% FontSizeAx = 24;
% FontSizeLb = 32;
% afFigurePosition = [0 2 30 20];
% axespos = [0.12 0.15 0.84 0.8];
% ylabelpos = [-0.075 0.5];
% xlabelpos = [0.5 -0.1];
% 
% t = ([1:Ncyc*Nstep+1] - 1)*delT;  
% delt = (t((Ncyc-2)*Nstep+1:end) - t((Ncyc-2)*Nstep+1))*f;
% 
% figure
% set(gcf, 'Units', 'centimeters','PaperPositionMode', 'auto','Position', afFigurePosition);
% set(gcf,'DefaultAxesFontSize',FontSizeAx,'DefaultAxesFontName','TimesNewRoman','DefaultAxesGridLineStyle','-.','DefaultAxesLineWidth',2,'DefaultAxesFontWeight','Normal')
% % 
% a_TE = [(zp(end,(Ncyc-2)*Nstep+1:end) - zp(Npanels/2+1,(Ncyc-2)*Nstep+1:end))/A_TE]';
% 
% Cpow_last = Cpow((Ncyc-2)*Nstep+1:end);
% hold on
%     plot(1/4*[1 1],[1.5*min(Cpow_last) 1.5*max(Cpow_last)],'-.k','linewidth',1)
%     plot(2/4*[1 1],[1.5*min(Cpow_last) 1.5*max(Cpow_last)],'-.k','linewidth',1)
%     plot(3/4*[1 1],[1.5*min(Cpow_last) 1.5*max(Cpow_last)],'-.k','linewidth',1)
%     plot(4/4*[1 1],[1.5*min(Cpow_last) 1.5*max(Cpow_last)],'-.k','linewidth',2)
%     plot(5/4*[1 1],[1.5*min(Cpow_last) 1.5*max(Cpow_last)],'-.k','linewidth',1)
%     plot(6/4*[1 1],[1.5*min(Cpow_last) 1.5*max(Cpow_last)],'-.k','linewidth',1)
%     plot(7/4*[1 1],[1.5*min(Cpow_last) 1.5*max(Cpow_last)],'-.k','linewidth',1)
% 
%     
%     plot([0 2],[0 0],'-k','linewidth',1)
%     plot(delt,3*a_TE,'--k','linewidth',3)
%     plot(delt,Cpow_last,'-b','linewidth',3)
% hold off
% 
% axis([0 2 1.2*min(Cpow_last) 1.2*max(Cpow_last)])
% set(gca, 'FontName', 'TimesNewRoman', 'FontSize', FontSizeAx)
% set(gca, 'Units', 'normalized', 'Position', axespos);
% xlabel('$$t/T$$','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 0,'Units', 'Normalize','Position', xlabelpos,'Interpreter', 'LaTeX');
% ylabel('$$C_p$$','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 90,'Units', 'Normalize','Position', ylabelpos,'Interpreter', 'LaTeX');

%% Intantaneous Lift
% 
% FontSizeAx = 24;
% FontSizeLb = 32;
% afFigurePosition = [0 2 30 20];
% axespos = [0.12 0.15 0.84 0.8];
% ylabelpos = [-0.075 0.5];
% xlabelpos = [0.5 -0.1];
% 
% t = ((1:Ncyc*Nstep+1) - 1)*delT;  
% delt = t*f;
% 
% figure
% set(gcf, 'Units', 'centimeters','PaperPositionMode', 'auto','Position', afFigurePosition);
% set(gcf,'DefaultAxesFontSize',FontSizeAx,'DefaultAxesFontName','TimesNewRoman','DefaultAxesGridLineStyle','-.','DefaultAxesLineWidth',2,'DefaultAxesFontWeight','Normal')
% % 
% % a_TE = [(zp(end,(Ncyc-2)*Nstep+1:end) - zp(Npanels/2+1,(Ncyc-2)*Nstep+1:end))/A_TE]';
% 
% hold on   
%     plot([0 Ncyc],[0 0],'-k','linewidth',1)
%     plot(delt,Ct_s,'-b','linewidth',3)
%     plot(delt,Ct_us,'-r','linewidth',3)
%     plot(delt,Ct,'-k','linewidth',3)
% hold off
% 
% axis([0 Ncyc min(Ct_us) max(Ct_us)])
% set(gca, 'FontName', 'TimesNewRoman', 'FontSize', FontSizeAx)
% set(gca, 'Units', 'normalized', 'Position', axespos);
% xlabel('$$t/T$$','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 0,'Units', 'Normalize','Position', xlabelpos,'Interpreter', 'LaTeX');
% ylabel('$$C_l$$','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 90,'Units', 'Normalize','Position', ylabelpos,'Interpreter', 'LaTeX');

%% Calculating the flowfield around the wing.
% nx = 50;
% nz = 50;
% 
% U = Uinf*ones(nz,nx);
% W = Winf*ones(nz,nx);
% 
% %         xf = linspace(-c/2 + x_b(i_t),3*c + x_b(i_t),nx)';
% %         zf = linspace(-3*c/4 + z_b(i_t),3*c/4 + z_b(i_t),nz)';
% 
% % Flow field for ground effect calculations
% xf = linspace(-c/2 + x_b(i_t),3*c + x_b(i_t),nx)';
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
% % c_range = max(max(abs(omega_y)));
% plot(xp(:,end),zp(:,end),'.-k')
% plot(xTE(:,end),zTE(:,end),'.-g')
% if grd == 1
%     plot([min(xp(:,i_t))-10*c; min(xp(:,i_t))+10*c],[0 0],'-k','linewidth',4)
% end
% % plot(xw,zw,'.-g')
% % plot(Xc,Zc,'xb')
% % streamline(Xf,Zf,U,W,Xf(:,1),Zf(:,1))
% 
% % bluevec = [linspace(0,1,100)' linspace(0,1,100)' ones(100,1)];
% % redvec = [ones(100,1) linspace(1,0,100)' linspace(1,0,100)'];
% % vec = [bluevec; redvec(2:end,:)];
% % colormap(vec)
%     
% % pcolor(Xstar,Zstar,-omega_y)
% quiver(Xf,Zf,U,W,'k')
% shading interp
% % caxis(0.7*[-c_range c_range])
% axis equal
% axis([-c/2 + x_b(end-1) 3*c + x_b(end-1) -c/10 c])
% 
% axis([-c/2 + x_b(end-1) 3*c + x_b(end-1) -c + z_b(end-1) c + z_b(end-1)])
% % title('2D Solution','interpreter','latex','fontsize',14,'fontname','TimesNewRoman')
% xlabel('$$x$$','interpreter','latex','fontsize',12,'fontname','TimesNewRoman')
% ylabel('$$z$$','interpreter','latex','fontsize',12,'fontname','TimesNewRoman')
% % print('-dpng','-r300',['2D_Wake_VF']);




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


