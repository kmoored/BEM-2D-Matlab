clear
clc

%% Nstep vary______________________________________________________________


load('Data2D_Convergence_N_vary150_Nstep_vary_Viscous.mat')
max_val = 15;

% load('Data2D_Convergence_N_vary_Nstep100_Viscous.mat')
% max_val = 15;
% load('Data2D_Convergence_N_vary_Nstep100_Viscous2.mat')
% max_val = 16;


for j = 1:max_val
    ind = find(Nstep_vec == 2*Nstep_vec(j));
    Ct_PctDiff(j) = abs(Ct_avg(j) - Ct_avg(ind))./Ct_avg(ind);
    Cpow_PctDiff(j) = abs(Cpow_avg(j) - Cpow_avg(ind))./Cpow_avg(ind);
    np_PctDiff(j) = abs(np(j) - np(ind))./np(ind);
end

% Plotting________________________________________________________________
FontSizeAx = 24;
FontSizeLb = 32;

crit = 0.01;
axis_vec = [Nstep_vec(1) Nstep_vec(max_val) 0 0.1];

% Thrust___________________________________________________________________
figure

set(gcf, 'Units', 'centimeters');
afFigurePosition = [1 1 25 10]; % [pos_x pos_y width_x width_y]
set(gcf, 'Position', afFigurePosition); % [left bottom width height]
set(gcf, 'PaperPositionMode', 'auto')
set(gcf,'DefaultAxesFontSize',FontSizeAx,'DefaultAxesFontName','TimesNewRoman','DefaultAxesGridLineStyle','-.','DefaultAxesLineWidth',2,'DefaultAxesFontWeight','Normal')

hold on
    plot(Nstep_vec(1:max_val), Ct_PctDiff, 's-k', 'LineWidth', 2, 'MarkerSize',6)
    plot([Nstep_vec(1) Nstep_vec(max_val)],crit*[1 1], '-.k', 'LineWidth', 1, 'MarkerSize',6)
    plot([Nstep_vec(1) Nstep_vec(max_val)],2*crit*[1 1], '-.k', 'LineWidth', 1, 'MarkerSize',6)
hold off

axis(axis_vec)
set(gca, 'FontName', 'TimesNewRoman', 'FontSize', FontSizeAx)
set(gca, 'Units', 'normalized', 'Position', [0.19 0.2 0.78 0.7]);
xlabel('$$N_{step}$$','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 0,'Units', 'Normalize','Position', [0.5 -0.15],'Interpreter', 'LaTeX');
ylabel('$$C_t$$','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 0,'Units', 'Normalize','Position', [-0.15 0.6],'Interpreter', 'LaTeX');


% Power___________________________________________________________________
figure

set(gcf, 'Units', 'centimeters');
afFigurePosition = [1 1 25 10]; % [pos_x pos_y width_x width_y]
set(gcf, 'Position', afFigurePosition); % [left bottom width height]
set(gcf, 'PaperPositionMode', 'auto')
set(gcf,'DefaultAxesFontSize',FontSizeAx,'DefaultAxesFontName','TimesNewRoman','DefaultAxesGridLineStyle','-.','DefaultAxesLineWidth',2,'DefaultAxesFontWeight','Normal')

hold on
    plot(Nstep_vec(1:max_val), Cpow_PctDiff, 's-k', 'LineWidth', 2, 'MarkerSize',6)
    plot([Nstep_vec(1) Nstep_vec(max_val)],crit*[1 1], '-.k', 'LineWidth', 1, 'MarkerSize',6)
    plot([Nstep_vec(1) Nstep_vec(max_val)],2*crit*[1 1], '-.k', 'LineWidth', 1, 'MarkerSize',6)
hold off

axis(axis_vec)
set(gca, 'FontName', 'TimesNewRoman', 'FontSize', FontSizeAx)
set(gca, 'Units', 'normalized', 'Position', [0.19 0.2 0.78 0.7]);
xlabel('$$N_{step}$$','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 0,'Units', 'Normalize','Position', [0.5 -0.15],'Interpreter', 'LaTeX');
ylabel('$$C_{pow}$$','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 0,'Units', 'Normalize','Position', [-0.15 0.6],'Interpreter', 'LaTeX');

% Efficiency___________________________________________________________________
figure

set(gcf, 'Units', 'centimeters');
afFigurePosition = [1 1 25 10]; % [pos_x pos_y width_x width_y]
set(gcf, 'Position', afFigurePosition); % [left bottom width height]
set(gcf, 'PaperPositionMode', 'auto')
set(gcf,'DefaultAxesFontSize',FontSizeAx,'DefaultAxesFontName','TimesNewRoman','DefaultAxesGridLineStyle','-.','DefaultAxesLineWidth',2,'DefaultAxesFontWeight','Normal')

hold on
    plot(Nstep_vec(1:max_val), np_PctDiff, 's-k', 'LineWidth', 2, 'MarkerSize',6)
    plot([Nstep_vec(1) Nstep_vec(max_val)],crit*[1 1], '-.k', 'LineWidth', 1, 'MarkerSize',6)
    plot([Nstep_vec(1) Nstep_vec(max_val)],2*crit*[1 1], '-.k', 'LineWidth', 1, 'MarkerSize',6)
hold off

axis(axis_vec)
set(gca, 'FontName', 'TimesNewRoman', 'FontSize', FontSizeAx)
set(gca, 'Units', 'normalized', 'Position', [0.19 0.2 0.78 0.7]);
xlabel('$$N_{step}$$','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 0,'Units', 'Normalize','Position', [0.5 -0.15],'Interpreter', 'LaTeX');
ylabel('$$\eta_p$$','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 0,'Units', 'Normalize','Position', [-0.15 0.6],'Interpreter', 'LaTeX');

%% N vary______________________________________________________________


load('Data2D_Convergence_N_vary_Nstep100_Viscous.mat')
max_val = 15;
% load('Data2D_Convergence_N_vary_Nstep100_Viscous2.mat')
% max_val = 16;


for j = 1:max_val
    ind = find(N_vec == 2*N_vec(j));
    Ct_PctDiff(j) = abs(Ct_avg(j) - Ct_avg(ind))./Ct_avg(ind);
    Cpow_PctDiff(j) = abs(Cpow_avg(j) - Cpow_avg(ind))./Cpow_avg(ind);
    np_PctDiff(j) = abs(np(j) - np(ind))./np(ind);
end

% Plotting________________________________________________________________
FontSizeAx = 24;
FontSizeLb = 32;

crit = 0.01;
axis_vec = [N_vec(1) N_vec(max_val) 0 0.1];

% Thrust___________________________________________________________________
figure

set(gcf, 'Units', 'centimeters');
afFigurePosition = [1 1 25 10]; % [pos_x pos_y width_x width_y]
set(gcf, 'Position', afFigurePosition); % [left bottom width height]
set(gcf, 'PaperPositionMode', 'auto')
set(gcf,'DefaultAxesFontSize',FontSizeAx,'DefaultAxesFontName','TimesNewRoman','DefaultAxesGridLineStyle','-.','DefaultAxesLineWidth',2,'DefaultAxesFontWeight','Normal')

hold on
    plot(N_vec(1:max_val), Ct_PctDiff, 's-k', 'LineWidth', 2, 'MarkerSize',6)
    plot([N_vec(1) N_vec(max_val)],crit*[1 1], '-.k', 'LineWidth', 1, 'MarkerSize',6)
    plot([N_vec(1) N_vec(max_val)],2*crit*[1 1], '-.k', 'LineWidth', 1, 'MarkerSize',6)
hold off

axis(axis_vec)
set(gca, 'FontName', 'TimesNewRoman', 'FontSize', FontSizeAx)
set(gca, 'Units', 'normalized', 'Position', [0.19 0.2 0.78 0.7]);
xlabel('$$N$$','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 0,'Units', 'Normalize','Position', [0.5 -0.15],'Interpreter', 'LaTeX');
ylabel('$$C_t$$','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 0,'Units', 'Normalize','Position', [-0.15 0.6],'Interpreter', 'LaTeX');


% Power___________________________________________________________________
figure

set(gcf, 'Units', 'centimeters');
afFigurePosition = [1 1 25 10]; % [pos_x pos_y width_x width_y]
set(gcf, 'Position', afFigurePosition); % [left bottom width height]
set(gcf, 'PaperPositionMode', 'auto')
set(gcf,'DefaultAxesFontSize',FontSizeAx,'DefaultAxesFontName','TimesNewRoman','DefaultAxesGridLineStyle','-.','DefaultAxesLineWidth',2,'DefaultAxesFontWeight','Normal')

hold on
    plot(N_vec(1:max_val), Cpow_PctDiff, 's-k', 'LineWidth', 2, 'MarkerSize',6)
    plot([N_vec(1) N_vec(max_val)],crit*[1 1], '-.k', 'LineWidth', 1, 'MarkerSize',6)
    plot([N_vec(1) N_vec(max_val)],2*crit*[1 1], '-.k', 'LineWidth', 1, 'MarkerSize',6)
hold off

axis(axis_vec)
set(gca, 'FontName', 'TimesNewRoman', 'FontSize', FontSizeAx)
set(gca, 'Units', 'normalized', 'Position', [0.19 0.2 0.78 0.7]);
xlabel('$$N$$','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 0,'Units', 'Normalize','Position', [0.5 -0.15],'Interpreter', 'LaTeX');
ylabel('$$C_{pow}$$','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 0,'Units', 'Normalize','Position', [-0.15 0.6],'Interpreter', 'LaTeX');

% Efficiency___________________________________________________________________
figure

set(gcf, 'Units', 'centimeters');
afFigurePosition = [1 1 25 10]; % [pos_x pos_y width_x width_y]
set(gcf, 'Position', afFigurePosition); % [left bottom width height]
set(gcf, 'PaperPositionMode', 'auto')
set(gcf,'DefaultAxesFontSize',FontSizeAx,'DefaultAxesFontName','TimesNewRoman','DefaultAxesGridLineStyle','-.','DefaultAxesLineWidth',2,'DefaultAxesFontWeight','Normal')

hold on
    plot(N_vec(1:max_val), np_PctDiff, 's-k', 'LineWidth', 2, 'MarkerSize',6)
    plot([N_vec(1) N_vec(max_val)],crit*[1 1], '-.k', 'LineWidth', 1, 'MarkerSize',6)
    plot([N_vec(1) N_vec(max_val)],2*crit*[1 1], '-.k', 'LineWidth', 1, 'MarkerSize',6)
hold off

axis(axis_vec)
set(gca, 'FontName', 'TimesNewRoman', 'FontSize', FontSizeAx)
set(gca, 'Units', 'normalized', 'Position', [0.19 0.2 0.78 0.7]);
xlabel('$$N$$','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 0,'Units', 'Normalize','Position', [0.5 -0.15],'Interpreter', 'LaTeX');
ylabel('$$\eta_p$$','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 0,'Units', 'Normalize','Position', [-0.15 0.6],'Interpreter', 'LaTeX');

%% val vary______________________________________________________________


load('Data2D_convergence_N150_Nstep150_val_vary_Viscous.mat')



% Plotting________________________________________________________________
FontSizeAx = 24;
FontSizeLb = 32;

crit = 0.01;
N_vec = 3:8;
% axis_vec = [N_vec(1) N_vec(max_val) 0 0.1];

% Thrust___________________________________________________________________
figure

set(gcf, 'Units', 'centimeters');
afFigurePosition = [1 1 25 10]; % [pos_x pos_y width_x width_y]
set(gcf, 'Position', afFigurePosition); % [left bottom width height]
set(gcf, 'PaperPositionMode', 'auto')
set(gcf,'DefaultAxesFontSize',FontSizeAx,'DefaultAxesFontName','TimesNewRoman','DefaultAxesGridLineStyle','-.','DefaultAxesLineWidth',2,'DefaultAxesFontWeight','Normal')

hold on
    plot(N_vec, Ct_avg, 's-k', 'LineWidth', 2, 'MarkerSize',6)
%     plot([N_vec(1) N_vec(max_val)],crit*[1 1], '-.k', 'LineWidth', 1, 'MarkerSize',6)
%     plot([N_vec(1) N_vec(max_val)],2*crit*[1 1], '-.k', 'LineWidth', 1, 'MarkerSize',6)
hold off

% axis(axis_vec)
set(gca, 'FontName', 'TimesNewRoman', 'FontSize', FontSizeAx)
set(gca, 'Units', 'normalized', 'Position', [0.19 0.2 0.78 0.7]);
xlabel('$$N$$','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 0,'Units', 'Normalize','Position', [0.5 -0.15],'Interpreter', 'LaTeX');
ylabel('$$C_t$$','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 0,'Units', 'Normalize','Position', [-0.15 0.6],'Interpreter', 'LaTeX');


% Power___________________________________________________________________
figure

set(gcf, 'Units', 'centimeters');
afFigurePosition = [1 1 25 10]; % [pos_x pos_y width_x width_y]
set(gcf, 'Position', afFigurePosition); % [left bottom width height]
set(gcf, 'PaperPositionMode', 'auto')
set(gcf,'DefaultAxesFontSize',FontSizeAx,'DefaultAxesFontName','TimesNewRoman','DefaultAxesGridLineStyle','-.','DefaultAxesLineWidth',2,'DefaultAxesFontWeight','Normal')

hold on
    plot(N_vec, Cpow_avg, 's-k', 'LineWidth', 2, 'MarkerSize',6)
%     plot([N_vec(1) N_vec(max_val)],crit*[1 1], '-.k', 'LineWidth', 1, 'MarkerSize',6)
%     plot([N_vec(1) N_vec(max_val)],2*crit*[1 1], '-.k', 'LineWidth', 1, 'MarkerSize',6)
hold off

% axis(axis_vec)
set(gca, 'FontName', 'TimesNewRoman', 'FontSize', FontSizeAx)
set(gca, 'Units', 'normalized', 'Position', [0.19 0.2 0.78 0.7]);
xlabel('$$N$$','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 0,'Units', 'Normalize','Position', [0.5 -0.15],'Interpreter', 'LaTeX');
ylabel('$$C_{pow}$$','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 0,'Units', 'Normalize','Position', [-0.15 0.6],'Interpreter', 'LaTeX');

% Efficiency___________________________________________________________________
figure

set(gcf, 'Units', 'centimeters');
afFigurePosition = [1 1 25 10]; % [pos_x pos_y width_x width_y]
set(gcf, 'Position', afFigurePosition); % [left bottom width height]
set(gcf, 'PaperPositionMode', 'auto')
set(gcf,'DefaultAxesFontSize',FontSizeAx,'DefaultAxesFontName','TimesNewRoman','DefaultAxesGridLineStyle','-.','DefaultAxesLineWidth',2,'DefaultAxesFontWeight','Normal')

hold on
    plot(N_vec, np, 's-k', 'LineWidth', 2, 'MarkerSize',6)
%     plot([N_vec(1) N_vec(max_val)],crit*[1 1], '-.k', 'LineWidth', 1, 'MarkerSize',6)
%     plot([N_vec(1) N_vec(max_val)],2*crit*[1 1], '-.k', 'LineWidth', 1, 'MarkerSize',6)
hold off

% axis(axis_vec)
set(gca, 'FontName', 'TimesNewRoman', 'FontSize', FontSizeAx)
set(gca, 'Units', 'normalized', 'Position', [0.19 0.2 0.78 0.7]);
xlabel('$$N$$','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 0,'Units', 'Normalize','Position', [0.5 -0.15],'Interpreter', 'LaTeX');
ylabel('$$\eta_p$$','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 0,'Units', 'Normalize','Position', [-0.15 0.6],'Interpreter', 'LaTeX');
