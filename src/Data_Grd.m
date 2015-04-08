clear
clc

St_vec = linspace(0.1, 0.45, 8);

% c = 0.12 m, Uinf = 0.06 m/s, NACA 0005, pure pitch, A/c = 0.25 (peak-to-peak)
% Case 1: d/c = infty.
% Case 2: d/c = 1/2.
% Case 3: d/c = 1/4.

%% Experimental Data

% Case 1, h/c =  inf
Ct_inf = [-0.0363 0.0101 0.0630 0.1553 0.2754 0.4247 0.5841 0.7811];
Cpow_inf = [0.2080 0.2984 0.6896 1.1665 2.1640 3.4181 5.6112 8.2545];
eta_inf = Ct_inf./Cpow_inf;

% Case 3, h/c = 0.25
CtExp_1 = [-0.0128 0.0149 0.0591 0.1232 0.2122 0.3338 0.5015 0.6719]';
CtExp_2 = [-0.0232 0.0107 0.0817 0.178 0.3095 0.4975 0.6964 0.9347]';
CtExp_3 = [-0.0128 0.0479 0.1425 0.2595 0.4075 0.5952 0.8675 1.1741]';


% load('Data2D_Grd_Inf_N210_Nstep80_Ep1e-3_Invscd.mat')
load('Data2D_Grd_Inf_N210_Nstep80_Ep1e-3_Visc.mat')

%% Plotting________________________________________________________________
FontSizeAx = 24;
FontSizeLb = 32;
afFigurePosition = [1 1 25 15];
ylabelpos = [-0.08 0.5];
xlabelpos = [0.5 -0.13];
axespos = [0.125 0.175 0.8 0.8];


% figure
% 
% set(gcf, 'Units', 'centimeters');
% set(gcf, 'Position', afFigurePosition); % [left bottom width height]
% set(gcf, 'PaperPositionMode', 'auto')
% set(gcf,'DefaultAxesFontSize',FontSizeAx,'DefaultAxesFontName','TimesNewRoman','DefaultAxesGridLineStyle','-.','DefaultAxesLineWidth',2,'DefaultAxesFontWeight','Normal')
% 
% hold on
%     plot(St_vec, np, 's-k', 'LineWidth', 2, 'MarkerSize',6)
% hold off
% 
% axis([min(St_vec) max(St_vec) 0 0.31])
% set(gca, 'FontName', 'TimesNewRoman', 'FontSize', FontSizeAx)
% set(gca, 'Units', 'normalized', 'Position', [0.19 0.2 0.78 0.7]);
% xlabel('$$St$$','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 0,'Units', 'Normalize','Position', [0.5 -0.15],'Interpreter', 'LaTeX');
% ylabel('$$\eta$$','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 0,'Units', 'Normalize','Position', [-0.15 0.6],'Interpreter', 'LaTeX');



figure
set(gcf, 'Units', 'centimeters','PaperPositionMode', 'auto','Position', afFigurePosition);
set(gcf,'DefaultAxesFontSize',FontSizeAx,'DefaultAxesFontName','TimesNewRoman','DefaultAxesGridLineStyle','-.','DefaultAxesLineWidth',2,'DefaultAxesFontWeight','Normal')

hold on

    plot(St_vec,Ct_inf,'ok','linewidth',2,'MarkerSize',10)
    plot(St_vec,Ct_avg,'-k','linewidth',2,'MarkerSize',10)
%     plot(St_vec,CtExp_2,'ob','linewidth',2,'MarkerSize',10)
%     plot(St_vec,CtExp_3,'og','linewidth',2,'MarkerSize',10)
    
%     plot(St_vec, Ct_avg, '-k', 'LineWidth', 2, 'MarkerSize',6)
    plot([0,0.5], [0,0], '-k', 'LineWidth', 1, 'MarkerSize',6)
      
hold off

axis([0 0.5 1.1*min(Ct_inf) 1.1*max(Ct_inf)])

set(gca, 'FontName', 'TimesNewRoman', 'FontSize', FontSizeAx)
set(gca, 'Units', 'normalized', 'Position', axespos);
legend('Exp','Panel','Location','NorthWest')
xlabel('$$St$$','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 0,'Units', 'Normalize','Position',xlabelpos,'Interpreter', 'LaTeX');
ylabel('$$C_t$$','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 90,'Units', 'Normalize','Position',ylabelpos,'Interpreter', 'LaTeX');




figure
set(gcf, 'Units', 'centimeters','PaperPositionMode', 'auto','Position', afFigurePosition);
set(gcf,'DefaultAxesFontSize',FontSizeAx,'DefaultAxesFontName','TimesNewRoman','DefaultAxesGridLineStyle','-.','DefaultAxesLineWidth',2,'DefaultAxesFontWeight','Normal')

hold on
    plot(St_vec, Cpow_inf, 'ok', 'LineWidth', 2, 'MarkerSize',10)
    plot(St_vec, Cpow_avg, '-k', 'LineWidth', 2, 'MarkerSize',10)
    plot([0,0.5], [0,0], '-k', 'LineWidth', 1, 'MarkerSize',6)
hold off

axis([0 0.5 -1.1*min(Cpow_inf) 1.1*max(Cpow_inf)])
set(gca, 'FontName', 'TimesNewRoman', 'FontSize', FontSizeAx)
set(gca, 'Units', 'normalized', 'Position', axespos);
legend('Exp','Panel','Location','NorthWest')
xlabel('$$St$$','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 0,'Units', 'Normalize','Position',xlabelpos,'Interpreter', 'LaTeX');
ylabel('$$C_{pow}$$','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 90,'Units', 'Normalize','Position',ylabelpos,'Interpreter', 'LaTeX');





		
