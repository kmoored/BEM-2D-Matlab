clear
close all
clc



%% Data -- Matrices are (12 frequencies X 3 amplitudes X 10 wall proximities) 

load('nearWallData.mat')
load('FixedV_Pitch_Grd_Data_v4v5_combined.mat')
% Ct_avg_low = Ct_avg;
% Cpow_avg_low = Cpow_avg;
% Cl_avg_low = Cl_avg;
% np_low = np;
% 
% load('FixedV_Pitch_Grd_Data_v4_combined.mat')
% Ct_avg(:,:,1:4) = Ct_avg_low;
% Cpow_avg(:,:,1:4) = Cpow_avg_low;
% Cl_avg(:,:,1:4) = Cl_avg_low;
% np(:,:,1:4) = np_low;

% %%%%%
% Ct_avg([1 2],3,3) = NaN*Ct_avg([1 2],3,3);
% Cpow_avg([1 2],3,3) = NaN*Cpow_avg([1 2],3,3);
% Cl_avg([1 2],3,3) = NaN*Cl_avg([1 2],3,3);
% np([1 2],3,3) = NaN*np([1 2],3,3);
% 
% %%%%%
% Ct_avg(1:3,3,2) = NaN*Ct_avg(1:3,3,2) ;
% Cpow_avg(1:3,3,2) = NaN*Cpow_avg(1:3,3,2);
% Cl_avg(1:3,3,2) = NaN*Cl_avg(1:3,3,2);
% np(1:3,3,2) = NaN*np(1:3,3,2);
% 
% Ct_avg(1:2,2,2) = NaN*Ct_avg(1:2,2,2) ;
% Cpow_avg(1:2,2,2) = NaN*Cpow_avg(1:2,2,2);
% Cl_avg(1:2,2,2) = NaN*Cl_avg(1:2,2,2);
% np(1:2,2,2) = NaN*np(1:2,2,2);
% 
% %%%%%
% % Ct_avg(:,[2 3],1) = NaN*Ct_avg(:,[2 3],1);
% % Cpow_avg(:,[2 3],1) = NaN*Cpow_avg(:,[2 3],1);
% % Cl_avg(:,[2 3],1) = NaN*Cl_avg(:,[2 3],1);
% % np(:,[2 3],1) = NaN*np(:,[2 3],1);
% 
% Ct_avg(:,[3],1) = NaN*Ct_avg(:,[3],1);
% Cpow_avg(:,[3],1) = NaN*Cpow_avg(:,[3],1);
% Cl_avg(:,[3],1) = NaN*Cl_avg(:,[3],1);
% np(:,[3],1) = NaN*np(:,[3],1);
% 
% Ct_avg([1:4],[2],1) = NaN*Ct_avg([1:4],[2],1);
% Cpow_avg([1:4],[2],1) = NaN*Cpow_avg([1:4],[2],1);
% Cl_avg([1:4],[2],1) = NaN*Cl_avg([1:4],[2],1);
% np([1:4],[2],1) = NaN*np([1:4],[2],1);
% 
% Ct_avg(end,[2],1) = NaN*Ct_avg(end,[2],1);
% Cpow_avg(end,[2],1) = NaN*Cpow_avg(end,[2],1);
% Cl_avg(end,[2],1) = NaN*Cl_avg(end,[2],1);
% np(end,[2],1) = NaN*np(end,[2],1);

% save('FixedV_Pitch_Grd_Data_v4v5_combined.mat','-v7.3')

%% Plotting
FontSizeAx = 24;
FontSizeLb = 32;
afFigurePosition = [1 1 20 15];
ylabelpos = [-0.075 0.5];
xlabelpos = [0.5 -0.14];
axespos = [0.12 0.2 0.85 0.75];

num_vec = 10:-1:1;
colorvec = -((num_vec-1)/9 - 1);

% Figure 1: Classic Thrust Coefficient vs St
figure
set(gcf, 'Units', 'centimeters');
set(gcf, 'Position', afFigurePosition); % [left bottom width height]
set(gcf, 'PaperPositionMode', 'auto')
set(gcf,'DefaultAxesFontSize',FontSizeAx,'DefaultAxesFontName','TimesNewRoman','DefaultAxesGridLineStyle','-.','DefaultAxesLineWidth',2,'DefaultAxesFontWeight','Normal')

hold on
num_vec = 1:8;
cv = -((num_vec-1)/7 - 1);
    
for i = 1:length(A_c_vec)
    A_c = A_c_vec(i);
    if i == 1
        vec = 'v-';
    elseif i == 2
        vec = 's-';
    elseif i == 3
        vec = '-';
    end

    for j = 1:length(d_c_vec)
        d_c = d_c_vec(j);
       
%         if i && j == 1
%         else
            plot(St_vec, Ct_avg(:,i,j), vec,'color',[(1 - colorvec(j)) 0 colorvec(j)], 'LineWidth', 2, 'MarkerSize',6)
%         end
    end
    
    St = f*(0.3125*c)/0.058;
    l = 1.3;
    m = 9;
    plot(St, Ct_d20_A25, 'd', 'Color', [cv(1) 0 (1 - cv(1))], 'LineWidth', l, 'MarkerSize', m)
    plot(St, Ct_d30_A25, '*', 'Color', [cv(2) 0 (1 - cv(2))], 'LineWidth', l, 'MarkerSize', m)
    plot(St, Ct_d40_A25, '^', 'Color', [cv(3) 0 (1 - cv(3))], 'LineWidth', l, 'MarkerSize', m)
    plot(St, Ct_d50_A25, 'v', 'Color', [cv(4) 0 (1 - cv(4))], 'LineWidth', l, 'MarkerSize', m)
    plot(St, Ct_d60_A25, '<', 'Color', [cv(5) 0 (1 - cv(5))], 'LineWidth', l, 'MarkerSize', m)
    plot(St, Ct_d70_A25, 'x', 'Color', [cv(6) 0 (1 - cv(6))], 'LineWidth', l, 'MarkerSize', m)
    plot(St, Ct_d80_A25, 's', 'Color', [cv(7) 0 (1 - cv(7))], 'LineWidth', l, 'MarkerSize', m)
    plot(St, Ct_free_A25, 'd', 'Color', [cv(8) 0 (1 - cv(8))], 'LineWidth', l, 'MarkerSize', m)

end


plot([0 20], [0 0], '-k', 'LineWidth', 1, 'MarkerSize',6)

hold off

% axis([0 0.09 -0.01 0.125])
axis([0 0.55 -0.1 1.5])
set(gca, 'FontName', 'TimesNewRoman', 'FontSize', FontSizeAx)
set(gca, 'Units', 'normalized', 'Position', axespos);
xlabel('$$St$$','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 0,'Units', 'Normalize','Position', xlabelpos,'Interpreter', 'LaTeX');
ylabel('$$C_t$$','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 90,'Units', 'Normalize','Position', ylabelpos,'Interpreter', 'LaTeX');
% print('-depsc', '-r600','T_v_f.eps');



% Figure 2: Lift Coefficient vs St
figure
set(gcf, 'Units', 'centimeters');
set(gcf, 'Position', afFigurePosition); % [left bottom width height]
set(gcf, 'PaperPositionMode', 'auto')
set(gcf,'DefaultAxesFontSize',FontSizeAx,'DefaultAxesFontName','TimesNewRoman','DefaultAxesGridLineStyle','-.','DefaultAxesLineWidth',2,'DefaultAxesFontWeight','Normal')

hold on

    
for i = 1:length(A_c_vec)
    A_c = A_c_vec(i);
    if i == 1
        vec = 'v-';
    elseif i == 2
        vec = 's-';
    elseif i == 3
        vec = '-';
    end

    for j = 1:length(d_c_vec)
        d_c = d_c_vec(j);
%         St_vec = f_vec*(A_c*c)/Qinf;
        plot(St_vec, Cl_avg(:,i,j), vec,'color',[(1 - colorvec(j)) 0 colorvec(j)], 'LineWidth', 2, 'MarkerSize',6)
    end

end
   
% plot([0 20], [0 0], '-k', 'LineWidth', 1, 'MarkerSize',6)

hold off

axis([0 0.4 -1 1])
set(gca, 'FontName', 'TimesNewRoman', 'FontSize', FontSizeAx)
set(gca, 'Units', 'normalized', 'Position', axespos);
xlabel('$$St$$','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 0,'Units', 'Normalize','Position', xlabelpos,'Interpreter', 'LaTeX');
ylabel('$$C_l$$','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 90,'Units', 'Normalize','Position', ylabelpos,'Interpreter', 'LaTeX');
% print('-depsc', '-r600','T_v_f.eps');


% Figure 3: Power Coefficient vs St
figure
set(gcf, 'Units', 'centimeters');
set(gcf, 'Position', afFigurePosition); % [left bottom width height]
set(gcf, 'PaperPositionMode', 'auto')
set(gcf,'DefaultAxesFontSize',FontSizeAx,'DefaultAxesFontName','TimesNewRoman','DefaultAxesGridLineStyle','-.','DefaultAxesLineWidth',2,'DefaultAxesFontWeight','Normal')

hold on

    
for i = 1:length(A_c_vec)
    A_c = A_c_vec(i);
    if i == 1
        vec = 'v-';
    elseif i == 2
        vec = 's-';
    elseif i == 3
        vec = '-';
    end

    for j = 1:length(d_c_vec)
        d_c = d_c_vec(j);
%         St_vec = f_vec*(A_c*c)/Qinf;
        plot(St_vec, Cpow_avg(:,i,j), vec,'color',[(1 - colorvec(j)) 0 colorvec(j)], 'LineWidth', 2, 'MarkerSize',6)
    end

end

 St = f*(0.3125*c)/0.058;
    l = 1.3;
    m = 9;
    plot(St, Cp_d20_A25, 'd', 'Color', [cv(1) 0 (1 - cv(1))], 'LineWidth', l, 'MarkerSize', m)
    plot(St, Cp_d30_A25, '*', 'Color', [cv(2) 0 (1 - cv(2))], 'LineWidth', l, 'MarkerSize', m)
    plot(St, Cp_d40_A25, '^', 'Color', [cv(3) 0 (1 - cv(3))], 'LineWidth', l, 'MarkerSize', m)
    plot(St, Cp_d50_A25, 'v', 'Color', [cv(4) 0 (1 - cv(4))], 'LineWidth', l, 'MarkerSize', m)
    plot(St, Cp_d60_A25, '<', 'Color', [cv(5) 0 (1 - cv(5))], 'LineWidth', l, 'MarkerSize', m)
    plot(St, Cp_d70_A25, 'x', 'Color', [cv(6) 0 (1 - cv(6))], 'LineWidth', l, 'MarkerSize', m)
    plot(St, Cp_d80_A25, 's', 'Color', [cv(7) 0 (1 - cv(7))], 'LineWidth', l, 'MarkerSize', m)
    plot(St, Cp_free_A25, 'd', 'Color', [cv(8) 0 (1 - cv(8))], 'LineWidth', l, 'MarkerSize', m)

% plot([0 20], [0 0], '-k', 'LineWidth', 1, 'MarkerSize',6)

hold off

axis([0 0.55 -0.01 8])
set(gca, 'FontName', 'TimesNewRoman', 'FontSize', FontSizeAx)
set(gca, 'Units', 'normalized', 'Position', axespos);
xlabel('$$St$$','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 0,'Units', 'Normalize','Position', xlabelpos,'Interpreter', 'LaTeX');
ylabel('$$C_p$$','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 90,'Units', 'Normalize','Position', ylabelpos,'Interpreter', 'LaTeX');
% print('-depsc', '-r600','T_v_f.eps');



% Figure 4: Efficiency vs St
figure
set(gcf, 'Units', 'centimeters');
set(gcf, 'Position', afFigurePosition); % [left bottom width height]
set(gcf, 'PaperPositionMode', 'auto')
set(gcf,'DefaultAxesFontSize',FontSizeAx,'DefaultAxesFontName','TimesNewRoman','DefaultAxesGridLineStyle','-.','DefaultAxesLineWidth',2,'DefaultAxesFontWeight','Normal')

hold on

    
for i = 1:length(A_c_vec)
    A_c = A_c_vec(i);
    if i == 1
        vec = 'v-';
    elseif i == 2
        vec = 's-';
    elseif i == 3
        vec = '-';
    end

    for j = 1:length(d_c_vec)
        d_c = d_c_vec(j);
%         St_vec = f_vec*(A_c*c)/Qinf;
        plot(St_vec, np(:,i,j), vec,'color',[(1 - colorvec(j)) 0 colorvec(j)], 'LineWidth', 2, 'MarkerSize',6)
    end

end

 St = f*(0.3125*c)/0.058;
    l = 1.3;
    m = 9;
    plot(St, eff_d20_A25, 'd', 'Color', [cv(1) 0 (1 - cv(1))], 'LineWidth', l, 'MarkerSize', m)
    plot(St, eff_d30_A25, '*', 'Color', [cv(2) 0 (1 - cv(2))], 'LineWidth', l, 'MarkerSize', m)
    plot(St, eff_d40_A25, '^', 'Color', [cv(3) 0 (1 - cv(3))], 'LineWidth', l, 'MarkerSize', m)
    plot(St, eff_d50_A25, 'v', 'Color', [cv(4) 0 (1 - cv(4))], 'LineWidth', l, 'MarkerSize', m)
    plot(St, eff_d60_A25, '<', 'Color', [cv(5) 0 (1 - cv(5))], 'LineWidth', l, 'MarkerSize', m)
    plot(St, eff_d70_A25, 'x', 'Color', [cv(6) 0 (1 - cv(6))], 'LineWidth', l, 'MarkerSize', m)
    plot(St, eff_d80_A25, 's', 'Color', [cv(7) 0 (1 - cv(7))], 'LineWidth', l, 'MarkerSize', m)
    plot(St, eff_free_A25, 'd', 'Color', [cv(8) 0 (1 - cv(8))], 'LineWidth', l, 'MarkerSize', m)


% plot([0 20], [0 0], '-k', 'LineWidth', 1, 'MarkerSize',6)

hold off

axis([0.05 0.55 0 0.35])
set(gca, 'FontName', 'TimesNewRoman', 'FontSize', FontSizeAx)
set(gca, 'Units', 'normalized', 'Position', axespos);
xlabel('$$St$$','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 0,'Units', 'Normalize','Position', xlabelpos,'Interpreter', 'LaTeX');
ylabel('$$\eta$$','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 90,'Units', 'Normalize','Position', ylabelpos,'Interpreter', 'LaTeX');
% print('-depsc', '-r600','T_v_f.eps');
