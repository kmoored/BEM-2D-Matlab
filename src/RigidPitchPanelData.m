clear
clc

St_vec = linspace(0.1, 0.45, 8);

% c = 0.12 m, Uinf = 0.06 m/s, NACA 0005, pure pitch, A/c = 0.25 (peak-to-peak)
% Case 1: d/c = infty.


%% Experimental Data

% Case 1, h/c =  inf
Ct_exp = [-0.0363 0.0101 0.0630 0.1553 0.2754 0.4247 0.5841 0.7811];
Cpow_exp = [0.2080 0.2984 0.6896 1.1665 2.1640 3.4181 5.6112 8.2545];
eta_exp = Ct_exp./Cpow_exp;


I_Data = load('Data2D_PitchComp_N210_Nstep80_Ep1e-3_Inviscid.mat');
V_Data = load('Data2D_PitchComp_N210_Nstep80_Ep1e-3_Viscous.mat');

%% Plotting________________________________________________________________
FontSizeAx = 24;
FontSizeLb = 32;
afFigurePosition = [1 1 25 15];
ylabelpos = [-0.08 0.5];
xlabelpos = [0.5 -0.13];
axespos = [0.125 0.175 0.8 0.8];


% Plotting the thrust
figure
set(gcf, 'Units', 'centimeters','PaperPositionMode', 'auto','Position', afFigurePosition);
set(gcf,'DefaultAxesFontSize',FontSizeAx,'DefaultAxesFontName','TimesNewRoman','DefaultAxesGridLineStyle','-.','DefaultAxesLineWidth',2,'DefaultAxesFontWeight','Normal')

hold on

    plot(St_vec,Ct_exp,'ok','linewidth',2,'MarkerSize',10)
    plot(St_vec,Ct_avg,'-k','linewidth',2,'MarkerSize',10)
    plot([0,0.5], [0,0], '-k', 'LineWidth', 1, 'MarkerSize',6)
      
hold off

axis([0 0.5 1.1*min(Ct_exp) 1.1*max(Ct_exp)])
set(gca, 'FontName', 'TimesNewRoman', 'FontSize', FontSizeAx)
set(gca, 'Units', 'normalized', 'Position', axespos);
legend('Exp','Panel','Location','NorthWest')
xlabel('$$St$$','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 0,'Units', 'Normalize','Position',xlabelpos,'Interpreter', 'LaTeX');
ylabel('$$C_t$$','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 90,'Units', 'Normalize','Position',ylabelpos,'Interpreter', 'LaTeX');


% Plotting the power input
figure
set(gcf, 'Units', 'centimeters','PaperPositionMode', 'auto','Position', afFigurePosition);
set(gcf,'DefaultAxesFontSize',FontSizeAx,'DefaultAxesFontName','TimesNewRoman','DefaultAxesGridLineStyle','-.','DefaultAxesLineWidth',2,'DefaultAxesFontWeight','Normal')

hold on
    plot(St_vec, Cpow_exp, 'ok', 'LineWidth', 2, 'MarkerSize',10)
    plot(St_vec, Cpow_avg, '-k', 'LineWidth', 2, 'MarkerSize',10)
    plot([0,0.5], [0,0], '-k', 'LineWidth', 1, 'MarkerSize',6)
hold off

axis([0 0.5 -1.1*min(Cpow_exp) 1.1*max(Cpow_exp)])
set(gca, 'FontName', 'TimesNewRoman', 'FontSize', FontSizeAx)
set(gca, 'Units', 'normalized', 'Position', axespos);
legend('Exp','Panel','Location','NorthWest')
xlabel('$$St$$','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 0,'Units', 'Normalize','Position',xlabelpos,'Interpreter', 'LaTeX');
ylabel('$$C_{pow}$$','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 90,'Units', 'Normalize','Position',ylabelpos,'Interpreter', 'LaTeX');





		
