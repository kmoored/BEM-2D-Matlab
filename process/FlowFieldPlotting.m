clear
close all
clc

%%%% Source paths on local machine
addpath('../src/geom')
addpath('../src/kine')
addpath('../src/util')

%% Parameters
% Parameters for identifying data file.
% c = 0.08;
% Qinf = 0.06;
% St_vec = linspace(0.01,0.4,12)';
% A_c_vec = [0.015 0.02 0.025]'/c;
% d_c_vec = [0.02 0.023333 0.026667 0.03 0.04 0.05 0.06 0.07 0.08 1000]'/c;

c = 0.08;
Qinf = 0.06;
f = 1;
A_c = 0.25;
d_c = 0.5033;

    
%% Loading Data
savefilename = ['_PitchGrd_f',num2str(f),...
            '_A_c',num2str(A_c),...
            '_d_c',num2str(d_c)];
        
load(['FlowfieldData/Figure9/Processed_',savefilename,'.mat']);   
load(['FlowfieldData/Figure9/Flowfield',savefilename,'.mat']);

%% Calculating
Phi = 1/2*pi;
i_t = 1 + ceil(Nstep/2/pi*Phi)

% Free-stream velocity
alpha = 0;
Uinf = Qinf*cos(alpha);
Winf = Qinf*sin(alpha);

nx = 91;
nz = 91;
U = Uinf*ones(nz,nx);
W = Winf*ones(nz,nx);
        
% Flow field for ground effect calculations
xf = linspace(-c/2 + x_b(1,i_t),6*c + x_b(1,i_t),nx)';
zf = linspace(z_b(1,i_t) - 1.5*c,z_b(1,i_t) + 1.5*c,nz)';

[Xf,Zf] = meshgrid(xf,zf);

u_p = Ut(:,:,i_t) - U;
w_p = Wt(:,:,i_t) - W;

[omega_y,Lambda2,Xstar,Zstar] = Lambda2Crit2D(Xf,Zf,u_p,w_p);
Xstar = (xf(3:end) + xf(1:end-2))/2;
Zstar = (zf(3:end) + zf(1:end-2))/2;

% Interpolating data points for imagesc plotting
num = 5;
omega_yI = interp2(omega_y,num);
[row,col] = size(omega_yI);
xstarI = linspace(Xstar(1),Xstar(end),col)';
zstarI = linspace(Zstar(1),Zstar(end),row)';


%% Plotting
% Plotting Vorticity
% Plotting airfoil with LE vortex sheets and the TE vortex
% sheet.

% Figure 1: Vorticity plots
flowfig = figure;
FontSizeAx = 24;

set(gcf, 'Units', 'centimeters');
afFigurePosition = [15 7 23 13]; % [pos_x pos_y width_x width_y]
set(gcf, 'Position', afFigurePosition); % [left bottom width height]
set(gcf, 'PaperPositionMode', 'auto')
set(gcf,'DefaultAxesFontSize',FontSizeAx,'DefaultAxesFontName','TimesNewRoman','DefaultAxesGridLineStyle','-.','DefaultAxesLineWidth',2,'DefaultAxesFontWeight','Normal')
set(gcf,'DefaultAxesTickDir', 'out')

bluevec = [linspace(0,1,100)' linspace(0,1,100)' ones(100,1)];
redvec = [ones(100,1) linspace(1,0,100)' linspace(1,0,100)'];
vec = [bluevec; redvec(2:end,:)];
colormap(vec)

%         colormap hot
hold on
axis image
%         pcolor(Xstar,Zstar,-Lambda2.*omega_y)
    imagesc(xstarI/c,zstarI/c,-omega_yI)
    caxis(10*[-1 1])
    colorbar

    quiver(Xf/c,Zf/c,u_p,w_p,'k')
    plot(xp(:,i_t)/c,zp(:,i_t)/c,'-k','linewidth',2)

    if grd == 1
        plot([min(xp(:,i_t))/c-10; min(xp(:,i_t))/c+10],[0 0],'-k','linewidth',4)
    end

axis equal

set(gca, 'Units', 'normalized', 'Position', [0.03 0.18 0.9 0.7]);
set(gca,'XTick',[0 1 2 3 4],'YTick',[0 0.5 1 1.5])  
%         set(gca,'XTick',[0 1 2 3 4],'YTick',[12499 12500 12501],'YTickLabel',{'-1','0','1'}) 

annotation(flowfig,'textbox',[0.875 0.845 0.11 0.16],'Interpreter','LaTeX','String',{'$$\omega,\;s^{-1}$$'},'HorizontalAlignment','center','FontSize',30,'FontName','TimesNewRoman','FitBoxToText','off','LineStyle','none');
xlabel('$$x/c$$','interpreter','latex','fontsize',30,'fontname','TimesNewRoman')
ylabel('$$z/c$$','interpreter','latex','fontsize',30,'fontname','TimesNewRoman')
%         axis([-1/4 + x_b(1,i_t)/c 3.5 + x_b(1,i_t)/c z_b(1,i_t)/c - 1 z_b(1,i_t)/c + 1]) 
axis([-1/4 + x_b(1,i_t)/c 3.5 + x_b(1,i_t)/c -0.05 z_b(1,i_t)/c + 1.45])
print('-depsc','-r600',['FlowFieldData/Figure9/Vort',savefilename,'.eps']);
save(['FlowFieldData/Figure9/FlowfieldData_Vort_1_2pi',savefilename,'.mat'],'-v7.3','xstarI','zstarI','omega_yI','xp','zp')    



%% Calculating time-averaged values
u_p_mean = mean(Ut(:,:,351:end),3) - U;
w_p_mean = mean(Wt(:,:,351:end),3) - W;
Vmean = sqrt((u_p_mean + U).^2 + (w_p_mean + W).^2);

% Interpolating data points for imagesc plotting
num = 3;
% xfI = linspace(xf(1),xf(end),num*71)';

VmeanI = interp2(Vmean,num);
[row,col] = size(VmeanI);
xfI = linspace(xf(1),xf(end),col)';
zfI = linspace(zf(1),zf(end),row)';

%% Plotting time-averaged
% Figure 2: Time-average velocity plots
flowfig = figure;
FontSizeAx = 24;

set(gcf, 'Units', 'centimeters');
afFigurePosition = [15 7 23 13]; % [pos_x pos_y width_x width_y]
set(gcf, 'Position', afFigurePosition); % [left bottom width height]
set(gcf, 'PaperPositionMode', 'auto')
set(gcf,'DefaultAxesFontSize',FontSizeAx,'DefaultAxesFontName','TimesNewRoman','DefaultAxesGridLineStyle','-.','DefaultAxesLineWidth',2,'DefaultAxesFontWeight','Normal')
set(gcf,'DefaultAxesTickDir', 'out')

bluevec = [linspace(0,1,100)' linspace(0,1,100)' ones(100,1)];
redvec = [ones(100,1) linspace(1,0,100)' linspace(1,0,100)'];
vec = [bluevec; redvec(2:end,:)];
colormap(vec)

        colormap hot

hold on
axis image
%         pcolor(Xstar,Zstar,-Lambda2.*omega_y)
    imagesc(xfI/c,zfI/c,VmeanI)
    caxis(0.1*[0 1])
    colorbar

    quiver(Xf/c,Zf/c,u_p_mean,w_p_mean,'k')
    plot(xp(:,i_t)/c,zp(:,i_t)/c,'-k','linewidth',2)

    if grd == 1
        plot([min(xp(:,i_t))/c-10; min(xp(:,i_t))/c+10],[0 0],'-k','linewidth',4)
    end

axis equal

set(gca, 'Units', 'normalized', 'Position', [0.03 0.18 0.9 0.7]);
set(gca,'XTick',[0 1 2 3 4],'YTick',[0 0.5 1 1.5])  
%         set(gca,'XTick',[0 1 2 3 4],'YTick',[12499 12500 12501],'YTickLabel',{'-1','0','1'}) 

annotation(flowfig,'textbox',[0.875 0.845 0.11 0.16],'Interpreter','LaTeX','String',{'$$\omega,\;s^{-1}$$'},'HorizontalAlignment','center','FontSize',30,'FontName','TimesNewRoman','FitBoxToText','off','LineStyle','none');
xlabel('$$x/c$$','interpreter','latex','fontsize',30,'fontname','TimesNewRoman')
ylabel('$$z/c$$','interpreter','latex','fontsize',30,'fontname','TimesNewRoman')
        axis([-1/4 + x_b(1,i_t)/c 3.5 + x_b(1,i_t)/c z_b(1,i_t)/c - 1 z_b(1,i_t)/c + 1]) 
% axis([-1/4 + x_b(1,i_t)/c 3.5 + x_b(1,i_t)/c -0.05 z_b(1,i_t)/c + 1.45])
print('-depsc','-r600',['FlowFieldData/Figure9/TA',savefilename,'.eps']);


%% Saving Data
save(['FlowFieldData/Figure9/FlowfieldData_TA',savefilename,'.mat'],'-v7.3','xstarI','zstarI','omega_yI','xfI','zfI','VmeanI','xp','zp')    
