clc; clear all; close all

%% Plot panel thrust data

clc; clear all; close all

St = linspace(0.1, 0.5, 9);
k = (6.28*4)*St;

Ct = [-0.0293 0.0237 0.0902 0.1471 0.2425 0.3312 0.4401 0.5638 0.7177];
Cp = [0.1054 0.2189 0.4370 0.6610 1.1792 1.9818 3.2044 4.3598 5.8552];
eff = Ct./Cp;


% St_m = [0.1 0.25 0.4 0.5];
% Ct_m = [-0.045329547619981 0.12523752373967 0.491505957680192 0.851863087565183];
% Cp_m = [0.069400646551917 0.692327696554625 2.71220975930077 5.2043200319988];
% eff_m = Ct_m./Cp_m;

St_m =	[0.1 0.25 0.35 0.5]';
Ct_m =	[-0.0378 0.1427 0.3749 0.8852]';
Cp_m = [0.0708 0.7024 1.8499 5.2401]';
eff_m =	[-10 20.32 20.26 16.89]'/100;

St_p = [0.025 0.05 0.0750 0.1 0.125 0.15 0.175 0.2 0.225 0.25 0.275 0.3 0.325 0.35 0.375 0.4 0.425 0.45 0.475 0.5];
Ct_p = [-0.05969785260716 -0.052235540346368 -0.038706405275524 -0.020086042566415 0.003371371871116 0.031562820336695 0.064739140322591 0.103050240638928 0.14611934191749 0.193990470711937 0.246904956455598 0.304947359822924 0.368047614676355 0.435060728119264 0.506838287900665 0.584040425298027 0.666634331565124 0.754494644219231 0.848024340136759 0.945023797030371];
Cp_p = [0.006445120965191 0.02580463408519 0.058225589889916 0.105040074499663 0.168491081211162 0.251455852148699 0.356911953580825 0.488543226890541 0.649632469888696 0.843473828742505 1.0738834428935 1.34426549765953 1.65835980039114 2.01900137058185 2.43018571595691 2.89493039521927 3.4169539245982 3.99915876464462 4.64626730264549 5.35946358692507];
eff_p = Ct_p./Cp_p;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         Thrust vs k          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
set(gcf, 'Units', 'centimeters');
afFigurePosition = [1 1 20 15]; % [pos_x pos_y width_x width_y]
set(gcf, 'Position', afFigurePosition); % [left bottom width height]
set(gcf, 'PaperPositionMode', 'auto')

hold on
plot(St, Ct, 's',   'Color', [0 0 0], 'LineWidth', 1.8, 'MarkerSize', 9)
plot(St_m, Ct_m, 'o-r', 'LineWidth', 1.8, 'MarkerSize', 9)
plot(St_p, Ct_p, '-b', 'LineWidth', 1.8, 'MarkerSize', 9)
hold off

axis([0.05 0.55 -0.1 1])
set(gca, 'FontName', 'Times', 'FontSize', 24)
set(gca, 'Units', 'normalized', 'Position', [0.15 0.2 0.8 0.7]);
       
xlabel('$$St$$', ...
    'FontName', 'Times', ...
    'FontUnit', 'points', ...
    'FontSize', 28, ...
    'FontWeight', 'normal', ...
    'Rotation', 0, ...
    'Units', 'Normalize', ...
    'Position', [0.5 -0.15], ...
    'Interpreter', 'LaTeX');

ylabel('$$C_T$$', ...
    'FontName', 'Times', ...
    'FontUnit', 'points', ...
    'FontSize', 28, ...
    'FontWeight', 'normal', ...
    'Rotation', 0, ...
    'Units', 'Normalize', ...
    'Position', [-0.12 0.45], ...
    'Interpreter', 'LaTeX');

printname = ['thrust_vs_k_single.eps'];
print('-depsc', '-r600',printname);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%          Power vs k          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
set(gcf, 'Units', 'centimeters');
afFigurePosition = [1 1 20 15]; % [pos_x pos_y width_x width_y]
set(gcf, 'Position', afFigurePosition); % [left bottom width height]
set(gcf, 'PaperPositionMode', 'auto')

hold on
plot(St, Cp, 's',   'Color', [0 0 0], 'LineWidth', 1.8, 'MarkerSize', 9)
plot(St_m, Cp_m, 'o-r', 'LineWidth', 1.8, 'MarkerSize', 9)
plot(St_p, Cp_p, '-b', 'LineWidth', 1.8, 'MarkerSize', 9)
hold off

axis([0 0.55 0 7])
set(gca, 'FontName', 'Times', 'FontSize', 24)
set(gca, 'Units', 'normalized', 'Position', [0.15 0.2 0.8 0.7]);


xlabel('$$St$$', ...
    'FontName', 'Times', ...
    'FontUnit', 'points', ...
    'FontSize', 28, ...
    'FontWeight', 'normal', ...
    'Rotation', 0, ...
    'Units', 'Normalize', ...
    'Position', [0.5 -0.15], ...
    'Interpreter', 'LaTeX');

ylabel('$$C_P$$', ...
    'FontName', 'Times', ...
    'FontUnit', 'points', ...
    'FontSize', 28, ...
    'FontWeight', 'normal', ...
    'Rotation', 0, ...
    'Units', 'Normalize', ...
    'Position', [-0.12 0.45], ...
    'Interpreter', 'LaTeX');

printname = ['power_vs_k_single.eps'];
print('-depsc', '-r600',printname);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       Efficiency vs k        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
set(gcf, 'Units', 'centimeters');
afFigurePosition = [1 1 20 15]; % [pos_x pos_y width_x width_y]
set(gcf, 'Position', afFigurePosition); % [left bottom width height]
set(gcf, 'PaperPositionMode', 'auto')

hold on
plot(St, eff, 's',   'Color', [0 0 0], 'LineWidth', 1.8, 'MarkerSize', 9)
plot(St_m, eff_m, 'o-r', 'LineWidth', 1.8, 'MarkerSize', 9)
plot(St_p, eff_p, '-b', 'LineWidth', 1.8, 'MarkerSize', 9)
hold off

axis([0.05 0.55 0 0.3])
set(gca, 'FontName', 'Times', 'FontSize', 24)
set(gca, 'Units', 'normalized', 'Position', [0.15 0.2 0.8 0.7]);


xlabel('$$St$$', ...
    'FontName', 'Times', ...
    'FontUnit', 'points', ...
    'FontSize', 28, ...
    'FontWeight', 'normal', ...
    'Rotation', 0, ...
    'Units', 'Normalize', ...
    'Position', [0.5 -0.15], ...
    'Interpreter', 'LaTeX');

ylabel('$$\eta$$', ...
    'FontName', 'Times', ...
    'FontUnit', 'points', ...
    'FontSize', 28, ...
    'FontWeight', 'normal', ...
    'Rotation', 0, ...
    'Units', 'Normalize', ...
    'Position', [-0.12 0.45], ...
    'Interpreter', 'LaTeX');

printname = ['efficiency_vs_k_single.eps'];
print('-depsc', '-r600',printname);

