clear
clc

% This data is for St = 0.5, A_c = 0.25, Re = 7906.
Data_1deg = load('Cpow_1deg.mat');
Data_15deg = load('Cpow_15deg.mat');
load('InstantPowerData.mat')

Cpow_1deg = Data_1deg.Cpow;
Cpow_15deg = Data_15deg.Cpow;

FontSizeAx = 24;
FontSizeLb = 32;
afFigurePosition = [0 2 30 20];
axespos = [0.12 0.15 0.84 0.8];
ylabelpos = [-0.075 0.5];
xlabelpos = [0.5 -0.1];

t = ([1:Ncyc*Nstep+1] - 1)*delT;  
delt = (t((Ncyc-2)*Nstep+1:end) - t((Ncyc-2)*Nstep+1))*f;

figure
set(gcf, 'Units', 'centimeters','PaperPositionMode', 'auto','Position', afFigurePosition);
set(gcf,'DefaultAxesFontSize',FontSizeAx,'DefaultAxesFontName','TimesNewRoman','DefaultAxesGridLineStyle','-.','DefaultAxesLineWidth',2,'DefaultAxesFontWeight','Normal')
% 
a_TE = [(zp(end,(Ncyc-2)*Nstep+1:end) - zp(Npanels/2+1,(Ncyc-2)*Nstep+1:end))/A_TE]';

Cpow_last = Cpow((Ncyc-2)*Nstep+1:end);
Cpow_last_1deg = Cpow_1deg((Ncyc-2)*Nstep+1:end);
Cpow_last_15deg = Cpow_15deg((Ncyc-2)*Nstep+1:end);
hold on
    plot(1/4*[1 1],[1.5*min(Cpow_last) 1.5*max(Cpow_last)],'-.k','linewidth',1)
    plot(2/4*[1 1],[1.5*min(Cpow_last) 1.5*max(Cpow_last)],'-.k','linewidth',1)
    plot(3/4*[1 1],[1.5*min(Cpow_last) 1.5*max(Cpow_last)],'-.k','linewidth',1)
    plot(4/4*[1 1],[1.5*min(Cpow_last) 1.5*max(Cpow_last)],'-.k','linewidth',2)
    plot(5/4*[1 1],[1.5*min(Cpow_last) 1.5*max(Cpow_last)],'-.k','linewidth',1)
    plot(6/4*[1 1],[1.5*min(Cpow_last) 1.5*max(Cpow_last)],'-.k','linewidth',1)
    plot(7/4*[1 1],[1.5*min(Cpow_last) 1.5*max(Cpow_last)],'-.k','linewidth',1)

    
    plot([0 2],[0 0],'-k','linewidth',1)
    plot(delt,3*a_TE,'--k','linewidth',3)
    plot(delt,Cpow_last_15deg,'-r','linewidth',3)
    plot(delt,Cpow_last_1deg,'-b','linewidth',3)
    plot(delt,Cpow_last,'-k','linewidth',3)
    
    
hold off

axis([0 2 1.2*min(Cpow_last) 1.2*max(Cpow_last)])
set(gca, 'FontName', 'TimesNewRoman', 'FontSize', FontSizeAx)
set(gca, 'Units', 'normalized', 'Position', axespos);
xlabel('$$t/T$$','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 0,'Units', 'Normalize','Position', xlabelpos,'Interpreter', 'LaTeX');
ylabel('$$C_p$$','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 90,'Units', 'Normalize','Position', ylabelpos,'Interpreter', 'LaTeX');
