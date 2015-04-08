clear
close all
clc


%% Data
load('FreeV/FreeV_Pitch_Data.mat')
Beta = 1;
Cd_bod = 0.05;
Ct_prop = 0.6;


%% Plotting
FontSizeAx = 24;
FontSizeLb = 32;
afFigurePosition = [1 1 22 15];
ylabelpos = [-0.075 0.5];
xlabelpos = [0.5 -0.14];
axespos = [0.12 0.2 0.85 0.75];

num_vec = 1:4;
colorvec = -((num_vec-1)/3 - 1);

delT = 1/f/Nstep;
t = ((1:Ncyc*Nstep+1) - 1)*delT;


figure

set(gcf, 'Units', 'centimeters');
set(gcf, 'Position', afFigurePosition); % [left bottom width height]
set(gcf, 'PaperPositionMode', 'auto')
set(gcf,'DefaultAxesFontSize',FontSizeAx,'DefaultAxesFontName','TimesNewRoman','DefaultAxesGridLineStyle','-.','DefaultAxesLineWidth',2,'DefaultAxesFontWeight','Normal')

hold on
for i = 1:length(A_bp_vec)
    for j = 1:length(M_star_vec)
        if i == 3 && j == 4
        else
            U_t = reshape(Q0(i,j,:),length(t),1);
            plot(t,U_t, '-','color',[(1 - colorvec(i)) 0 colorvec(i)], 'LineWidth', 2, 'MarkerSize',6)
        end
    end
end
hold off

axis([0 15 0 2.75])
set(gca, 'FontName', 'TimesNewRoman', 'FontSize', FontSizeAx)
set(gca, 'Units', 'normalized', 'Position', axespos);
xlabel('$$t$$, s','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 0,'Units', 'Normalize','Position', xlabelpos,'Interpreter', 'LaTeX');
ylabel('$$U(t)$$, m/s','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 90,'Units', 'Normalize','Position', ylabelpos,'Interpreter', 'LaTeX');
% print('-depsc', '-r600','Ut_v_t.eps');


figure

set(gcf, 'Units', 'centimeters');
set(gcf, 'Position', afFigurePosition); % [left bottom width height]
set(gcf, 'PaperPositionMode', 'auto')
set(gcf,'DefaultAxesFontSize',FontSizeAx,'DefaultAxesFontName','TimesNewRoman','DefaultAxesGridLineStyle','-.','DefaultAxesLineWidth',2,'DefaultAxesFontWeight','Normal')

hold on
for i = length(A_bp_vec):-1:1
    for j = 1:length(M_star_vec)
        if i == 3 && j == 4
        else
            U_t = reshape(Q0(i,j,:),length(t),1);
            plot(t/(M_star_vec(j)/(Ct_prop*f)*sqrt(2*(1/A_bp_vec(i))*Ct_prop/Cd_bod)),U_t/(f*A_c*c*sqrt(2*(1/A_bp_vec(i))*Ct_prop/Cd_bod)), '-','color',[(1 - colorvec(i)) 0 colorvec(i)], 'LineWidth', 2, 'MarkerSize',6)
        end
    end
end
hold off

axis([0 10 0 3])
set(gca, 'FontName', 'TimesNewRoman', 'FontSize', FontSizeAx)
set(gca, 'Units', 'normalized', 'Position', axespos);
xlabel('$$\hat{t}$$','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 0,'Units', 'Normalize','Position', xlabelpos,'Interpreter', 'LaTeX');
ylabel('$$\hat{U}$$','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 90,'Units', 'Normalize','Position', ylabelpos,'Interpreter', 'LaTeX');
% print('-depsc', '-r600','Uhat_v_that.eps');


figure

set(gcf, 'Units', 'centimeters');
set(gcf, 'Position', afFigurePosition); % [left bottom width height]
set(gcf, 'PaperPositionMode', 'auto')
set(gcf,'DefaultAxesFontSize',FontSizeAx,'DefaultAxesFontName','TimesNewRoman','DefaultAxesGridLineStyle','-.','DefaultAxesLineWidth',2,'DefaultAxesFontWeight','Normal')

hold on
% for i = length(A_bp_vec):-1:1
    for j = 1:length(M_star_vec)
        plot(A_bp_vec,Ec(:,j),'s-','color',[(1 - colorvec(j)) 0 colorvec(j)], 'LineWidth', 2, 'MarkerSize',6) 
    end
% end
hold off

% axis([0 10 0 3])
set(gca, 'FontName', 'TimesNewRoman', 'FontSize', FontSizeAx)
set(gca, 'Units', 'normalized', 'Position', axespos);
xlabel('$$S_{bp}$$','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 0,'Units', 'Normalize','Position', xlabelpos,'Interpreter', 'LaTeX');
ylabel('$$\zeta$$','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 90,'Units', 'Normalize','Position', ylabelpos,'Interpreter', 'LaTeX');
% print('-depsc', '-r600','Uhat_v_that.eps');
