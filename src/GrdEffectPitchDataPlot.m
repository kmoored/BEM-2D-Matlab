clear
clc

% Plotting Ground Effect Data
% load('Data2D_GrdEffect_N150_Nstep150_Viscous.mat')
% load('Data2D_GrdEffect_N150_Nstep200_Viscous.mat')
load('Data2D_Grd_N150_Nstep150_Viscous_2.mat')

% Plotting________________________________________________________________
FontSizeAx = 24;
FontSizeLb = 32;
% axis_vec = [0 0.5 0 0.1];

%% Thrust___________________________________________________________________
figure

set(gcf, 'Units', 'centimeters');
afFigurePosition = [1 1 25 20]; % [pos_x pos_y width_x width_y]
set(gcf, 'Position', afFigurePosition); % [left bottom width height]
set(gcf, 'PaperPositionMode', 'auto')
set(gcf,'DefaultAxesFontSize',FontSizeAx,'DefaultAxesFontName','TimesNewRoman','DefaultAxesGridLineStyle','-.','DefaultAxesLineWidth',2,'DefaultAxesFontWeight','Normal')

hold on
    plot(St_vec,Ct_avg(:,1), '-k', 'LineWidth', 2, 'MarkerSize',6)
    plot(St_vec,Ct_avg(:,2), '-b', 'LineWidth', 2, 'MarkerSize',6)
%     plot(St_vec,Ct_avg(:,3), '-r', 'LineWidth', 2, 'MarkerSize',6)
%     plot(St_vec,Ct_avg(:,4), '-m', 'LineWidth', 2, 'MarkerSize',6)
%     plot(St_vec,Ct_avg(:,5), '-g', 'LineWidth', 2, 'MarkerSize',6)
%     plot(St_vec,Ct_avg(:,6), 'o-y', 'LineWidth', 2, 'MarkerSize',6)
hold off

axis([0 0.5 -0.1 1.3])
set(gca, 'FontName', 'TimesNewRoman', 'FontSize', FontSizeAx)
set(gca, 'Units', 'normalized', 'Position', [0.19 0.2 0.78 0.7]);
% legend('d/c = 100','d/c = 1','d/c = 5/6','d/c = 2/3','d/c = 1/2','d/c = 1/3','Location','NorthWest')
xlabel('$$St$$','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 0,'Units', 'Normalize','Position', [0.5 -0.15],'Interpreter', 'LaTeX');
ylabel('$$C_t$$','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 0,'Units', 'Normalize','Position', [-0.15 0.6],'Interpreter', 'LaTeX');

%% Lift___________________________________________________________________
figure

set(gcf, 'Units', 'centimeters');
afFigurePosition = [1 1 25 20]; % [pos_x pos_y width_x width_y]
set(gcf, 'Position', afFigurePosition); % [left bottom width height]
set(gcf, 'PaperPositionMode', 'auto')
set(gcf,'DefaultAxesFontSize',FontSizeAx,'DefaultAxesFontName','TimesNewRoman','DefaultAxesGridLineStyle','-.','DefaultAxesLineWidth',2,'DefaultAxesFontWeight','Normal')

hold on
    plot(St_vec,Cl_avg(:,1), '-k', 'LineWidth', 2, 'MarkerSize',6)
    plot(St_vec,Cl_avg(:,2), '-b', 'LineWidth', 2, 'MarkerSize',6)
%     plot(St_vec,Cl_avg(:,3), '-r', 'LineWidth', 2, 'MarkerSize',6)
%     plot(St_vec,Cl_avg(:,4), '-m', 'LineWidth', 2, 'MarkerSize',6)
%     plot(St_vec,Cl_avg(:,5), '-g', 'LineWidth', 2, 'MarkerSize',6)
hold off

axis([0 0.5 -0.5 0.05])
set(gca, 'FontName', 'TimesNewRoman', 'FontSize', FontSizeAx)
set(gca, 'Units', 'normalized', 'Position', [0.19 0.2 0.78 0.7]);
% legend('d/c = 100','d/c = 1','d/c = 5/6','d/c = 2/3','d/c = 1/2','d/c = 1/3','Location','SouthWest')
xlabel('$$St$$','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 0,'Units', 'Normalize','Position', [0.5 -0.15],'Interpreter', 'LaTeX');
ylabel('$$C_l$$','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 0,'Units', 'Normalize','Position', [-0.15 0.6],'Interpreter', 'LaTeX');

%% Power___________________________________________________________________
figure

set(gcf, 'Units', 'centimeters');
afFigurePosition = [1 1 25 20]; % [pos_x pos_y width_x width_y]
set(gcf, 'Position', afFigurePosition); % [left bottom width height]
set(gcf, 'PaperPositionMode', 'auto')
set(gcf,'DefaultAxesFontSize',FontSizeAx,'DefaultAxesFontName','TimesNewRoman','DefaultAxesGridLineStyle','-.','DefaultAxesLineWidth',2,'DefaultAxesFontWeight','Normal')

hold on
    plot(St_vec,Cpow_avg(:,1), '-k', 'LineWidth', 2, 'MarkerSize',6)
    plot(St_vec,Cpow_avg(:,2), '-b', 'LineWidth', 2, 'MarkerSize',6)
%     plot(St_vec,Cpow_avg(:,3), '-r', 'LineWidth', 2, 'MarkerSize',6)
%     plot(St_vec,Cpow_avg(:,4), '-m', 'LineWidth', 2, 'MarkerSize',6)
%     plot(St_vec,Cpow_avg(:,5), '-g', 'LineWidth', 2, 'MarkerSize',6)
hold off

% axis(axis_vec)
set(gca, 'FontName', 'TimesNewRoman', 'FontSize', FontSizeAx)
set(gca, 'Units', 'normalized', 'Position', [0.19 0.2 0.78 0.7]);
% legend('d/c = 100','d/c = 1','d/c = 5/6','d/c = 2/3','d/c = 1/2','d/c = 1/3','Location','NorthWest')
xlabel('$$St$$','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 0,'Units', 'Normalize','Position', [0.5 -0.15],'Interpreter', 'LaTeX');
ylabel('$$C_{pow}$$','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 0,'Units', 'Normalize','Position', [-0.15 0.6],'Interpreter', 'LaTeX');


%% Efficiency___________________________________________________________________
figure

set(gcf, 'Units', 'centimeters');
afFigurePosition = [1 1 25 20]; % [pos_x pos_y width_x width_y]
set(gcf, 'Position', afFigurePosition); % [left bottom width height]
set(gcf, 'PaperPositionMode', 'auto')
set(gcf,'DefaultAxesFontSize',FontSizeAx,'DefaultAxesFontName','TimesNewRoman','DefaultAxesGridLineStyle','-.','DefaultAxesLineWidth',2,'DefaultAxesFontWeight','Normal')

hold on
    plot(St_vec,np(:,1), '-k', 'LineWidth', 2, 'MarkerSize',6)
    plot(St_vec,np(:,2), '-b', 'LineWidth', 2, 'MarkerSize',6)
%     plot(St_vec,np(:,3), '-r', 'LineWidth', 2, 'MarkerSize',6)
%     plot(St_vec,np(:,4), '-m', 'LineWidth', 2, 'MarkerSize',6)
%     plot(St_vec,np(:,5), '-g', 'LineWidth', 2, 'MarkerSize',6)
hold off

axis([0 0.5 0 0.4])
set(gca, 'FontName', 'TimesNewRoman', 'FontSize', FontSizeAx)
set(gca, 'Units', 'normalized', 'Position', [0.19 0.2 0.78 0.7]);
% legend('d/c = 100','d/c = 1','d/c = 5/6','d/c = 2/3','d/c = 1/2','d/c = 1/3','Location','NorthWest')
xlabel('$$St$$','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 0,'Units', 'Normalize','Position', [0.5 -0.15],'Interpreter', 'LaTeX');
ylabel('$$\eta_p$$','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 0,'Units', 'Normalize','Position', [-0.15 0.6],'Interpreter', 'LaTeX');
