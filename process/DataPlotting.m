clear
close all
clc

%%%% Source paths on local machine
addpath('../src/geom')
addpath('../src/kine')
addpath('../src/util')

%% Data

load('FixedV/FixedV_Pitch_Data.mat')
% load('FixedV/Fixedv_Pitch_Data_Re_1e4_2.mat')
% load('FixedV/fixedv_pitch_data_re_1e4.mat')
rho = 998;                                    % kg/m^3
nu = 1.004e-6;                                % m^2/s

Beta = 1;
% Beta = A_c_vec.^(1/2);
Sp = 2*c^2;


% Pow_avg = Cpow_avg.*(1/2*rho*c*U^3);
% A = 2*A;

%% Plotting
FontSizeAx = 24;
FontSizeLb = 32;
afFigurePosition = [1 1 20 15];
ylabelpos = [-0.075 0.5];
xlabelpos = [0.5 -0.14];
axespos = [0.12 0.2 0.85 0.75];

num_vec = 1:5;
colorvec = -((num_vec-1)/4 - 1);

T_avg = Tnet_avg + D_avg;
np_ND = T_avg*U./Pow_avg;

figure

set(gcf, 'Units', 'centimeters');
set(gcf, 'Position', afFigurePosition); % [left bottom width height]
set(gcf, 'PaperPositionMode', 'auto')
set(gcf,'DefaultAxesFontSize',FontSizeAx,'DefaultAxesFontName','TimesNewRoman','DefaultAxesGridLineStyle','-.','DefaultAxesLineWidth',2,'DefaultAxesFontWeight','Normal')

hold on
for j = 1:length(Qinf_vec)
    Qinf = Qinf_vec(j);
    for i = 1:length(A_c_vec)
        A_c = A_c_vec(i);
        if j == 1
            vec = 's-';
        elseif j == 2
            vec = 's--';
        elseif j == 3
            vec = 's:';
        elseif j == 4
            vec = 's-.';
        else 
            vec = 'o-';
        plot((Qinf/2/pi/c)*k_vec, Tnet_avg(:,i), 's-','color',[(1 - colorvec(i)) 0 colorvec(i)], 'LineWidth', 2, 'MarkerSize',6)
    end
end
plot([0.0007 0.08], [0 0], '-k', 'LineWidth', 2, 'MarkerSize',6)
hold off

% axis([0 0.09 -0.01 0.125])
set(gca, 'FontName', 'TimesNewRoman', 'FontSize', FontSizeAx)
set(gca, 'Units', 'normalized', 'Position', axespos);
xlabel('$$f, Hz$$','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 0,'Units', 'Normalize','Position', xlabelpos,'Interpreter', 'LaTeX');
ylabel('$$T, N$$','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 90,'Units', 'Normalize','Position', ylabelpos,'Interpreter', 'LaTeX');
% print('-depsc', '-r600','T_v_f.eps');



figure

set(gcf, 'Units', 'centimeters');
set(gcf, 'Position', afFigurePosition); % [left bottom width height]
set(gcf, 'PaperPositionMode', 'auto')
set(gcf,'DefaultAxesFontSize',FontSizeAx,'DefaultAxesFontName','TimesNewRoman','DefaultAxesGridLineStyle','-.','DefaultAxesLineWidth',2,'DefaultAxesFontWeight','Normal')

hold on
for i = 1:length(A_c_vec)
    plot(St_vec, T_avg(:,i), 's-','color',[(1 - colorvec(i)) 0 colorvec(i)], 'LineWidth', 2, 'MarkerSize',6)
end
% plot([0.05 0.5], [0 0], '-k', 'LineWidth', 2, 'MarkerSize',6)
hold off

% axis([0.05 0.5 -0.01 0.2])
set(gca, 'FontName', 'TimesNewRoman', 'FontSize', FontSizeAx)
set(gca, 'Units', 'normalized', 'Position', axespos);
xlabel('$$St$$','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 0,'Units', 'Normalize','Position', xlabelpos,'Interpreter', 'LaTeX');
ylabel('$$T, N$$','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 90,'Units', 'Normalize','Position', ylabelpos,'Interpreter', 'LaTeX');
% print('-depsc', '-r600','T_v_St.eps');



figure

set(gcf, 'Units', 'centimeters');
set(gcf, 'Position', afFigurePosition); % [left bottom width height]
set(gcf, 'PaperPositionMode', 'auto')
set(gcf,'DefaultAxesFontSize',FontSizeAx,'DefaultAxesFontName','TimesNewRoman','DefaultAxesGridLineStyle','-.','DefaultAxesLineWidth',2,'DefaultAxesFontWeight','Normal')

hold on
for i = 1:length(A_c_vec)
    beta = Beta(i);
    if i == 1 
        plot(2*pi*f(1:8,1)*c/U, Tnet_avg(1:8,i)./(rho*c*Sp*A(1:8,i).^2.*f(1:8,1).^2), 's-','color',[(1 - colorvec(i)) 0 colorvec(i)], 'LineWidth', 2, 'MarkerSize',6)
    else
        plot(2*pi*f(:,i)*c/U, Tnet_avg(:,i)./(rho*c*Sp*A(:,i).^2.*f(:,i).^2), 's-','color',[(1 - colorvec(i)) 0 colorvec(i)], 'LineWidth', 2, 'MarkerSize',6)
    end
end
% plot([0.05 0.5], [0 0], '-k', 'LineWidth', 2, 'MarkerSize',6)
plot([0 41], 2.25*[1 1], '--k', 'LineWidth', 2, 'MarkerSize',6)
hold off

axis([0 41 0 4])
% axis([0.05 0.5 0 4])
set(gca, 'FontName', 'TimesNewRoman', 'FontSize', FontSizeAx)
set(gca, 'Units', 'normalized', 'Position', axespos);
xlabel('$$k$$','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 0,'Units', 'Normalize','Position', xlabelpos,'Interpreter', 'LaTeX');
ylabel('$$C_T$$','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 90,'Units', 'Normalize','Position', ylabelpos,'Interpreter', 'LaTeX');
% print('-depsc', '-r600','Ct_v_k.eps');




figure

set(gcf, 'Units', 'centimeters');
set(gcf, 'Position', afFigurePosition); % [left bottom width height]
set(gcf, 'PaperPositionMode', 'auto')
set(gcf,'DefaultAxesFontSize',FontSizeAx,'DefaultAxesFontName','TimesNewRoman','DefaultAxesGridLineStyle','-.','DefaultAxesLineWidth',2,'DefaultAxesFontWeight','Normal')

hold on
for i = length(A_c_vec):-1:1
    beta = Beta(i);
    if i == 1 
            plot(2*pi*f(1:8,1)*c/U, T_avg(1:8,i)./(rho*c*Sp*A(1:8,i).^2.*f(1:8,1).^2), 's-','color',[(1 - colorvec(i)) 0 colorvec(i)], 'LineWidth', 2, 'MarkerSize',6)
    else
        plot(2*pi*f(:,i)*c/U, T_avg(:,i)./(rho*c*Sp*A(:,i).^2.*f(:,i).^2), 's-','color',[(1 - colorvec(i)) 0 colorvec(i)], 'LineWidth', 2, 'MarkerSize',6)
    end
%     plot(St_vec, T_avg(:,i)./(rho*c*2*c*A(:,i).^2.*f(:,i).^2), 's-','color',[(1 - colorvec(i)) 0 colorvec(i)], 'LineWidth', 2, 'MarkerSize',6)
end
plot([0 41], 2.25*[1 1], '--k', 'LineWidth', 2, 'MarkerSize',6)

% plot([0 50], 2.25*[1 1], '--k', 'LineWidth', 2, 'MarkerSize',6)

hold off

axis([0 40 0 4])
% axis([0.05 0.5 0 3])
set(gca, 'FontName', 'TimesNewRoman', 'FontSize', FontSizeAx)
set(gca, 'Units', 'normalized', 'Position', axespos);
xlabel('$$k$$','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 0,'Units', 'Normalize','Position', xlabelpos,'Interpreter', 'LaTeX');
ylabel('$$C_T$$','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 90,'Units', 'Normalize','Position', ylabelpos,'Interpreter', 'LaTeX');
% print('-depsc', '-r600','Ct_v_k_NoD.eps');



% figure
% 
% set(gcf, 'Units', 'centimeters');
% set(gcf, 'Position', afFigurePosition); % [left bottom width height]
% set(gcf, 'PaperPositionMode', 'auto')
% set(gcf,'DefaultAxesFontSize',FontSizeAx,'DefaultAxesFontName','TimesNewRoman','DefaultAxesGridLineStyle','-.','DefaultAxesLineWidth',2,'DefaultAxesFontWeight','Normal')
% 
% hold on
% for i = length(A_c_vec)-1:-1:1
% %     if i == 1 
% %             plot(2*pi*f(1:8,1)*c/U, T_avg(1:8,i)./(rho*c*2*c*A(1:8,i).^2.*f(1:8,1).^2), 's-','color',[(1 - colorvec(i)) 0 colorvec(i)], 'LineWidth', 2, 'MarkerSize',6)
% %     else
% %         plot(2*pi*f(:,i)*c/U, T_avg(:,i)./(rho*c*2*c*A(:,i).^2.*f(:,i).^2), 's-','color',[(1 - colorvec(i)) 0 colorvec(i)], 'LineWidth', 2, 'MarkerSize',6)
% %     end
%     plot(St_vec, T_avg(:,i)./(rho*c*2*c*A(:,i).^2.*f(:,i).^2), 's-','color',[(1 - colorvec(i)) 0 colorvec(i)], 'LineWidth', 2, 'MarkerSize',6)
% end
% plot([0 41], 2.25*[1 1], '--k', 'LineWidth', 2, 'MarkerSize',6)
% 
% % plot([0 50], 2.25*[1 1], '--k', 'LineWidth', 2, 'MarkerSize',6)
% 
% hold off
% 
% % axis([0 41 0 3])
% axis([0.05 0.8 0 3])
% set(gca, 'FontName', 'TimesNewRoman', 'FontSize', FontSizeAx)
% set(gca, 'Units', 'normalized', 'Position', axespos);
% xlabel('$$k$$','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 0,'Units', 'Normalize','Position', xlabelpos,'Interpreter', 'LaTeX');
% ylabel('$$C_T$$','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 90,'Units', 'Normalize','Position', ylabelpos,'Interpreter', 'LaTeX');
% % print('-depsc', '-r600','Ct_v_k_NoD.eps');

figure

set(gcf, 'Units', 'centimeters');
set(gcf, 'Position', afFigurePosition); % [left bottom width height]
set(gcf, 'PaperPositionMode', 'auto')
set(gcf,'DefaultAxesFontSize',FontSizeAx,'DefaultAxesFontName','TimesNewRoman','DefaultAxesGridLineStyle','-.','DefaultAxesLineWidth',2,'DefaultAxesFontWeight','Normal')

hold on
for i = 1:length(A_c_vec)
    beta = Beta(i);
%     plot(2*pi*f(:,i)*c/U, Pow_avg(:,i)./(rho*c*2*c*U^2*A(:,i).*f(:,i)), 's-','color',[(1 - colorvec(i)) 0 colorvec(i)], 'LineWidth', 2, 'MarkerSize',6)
    plot(2*pi*f(:,i)*c/U, Pow_avg(:,i), 's-','color',[(1 - colorvec(i)) 0 colorvec(i)], 'LineWidth', 2, 'MarkerSize',6)

end
% plot([0.05 0.5], [0 0], '-k', 'LineWidth', 2, 'MarkerSize',6)
% plot([0.05 0.5], 0.6*[1 1], '--k', 'LineWidth', 2, 'MarkerSize',6)
hold off

% axis([0 40 0 40])
% axis([0 80 0 120])
set(gca, 'FontName', 'TimesNewRoman', 'FontSize', FontSizeAx)
set(gca, 'Units', 'normalized', 'Position', axespos);
xlabel('$$k$$','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 0,'Units', 'Normalize','Position', xlabelpos,'Interpreter', 'LaTeX');
ylabel('$$\mathcal{P}$$','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 90,'Units', 'Normalize','Position', ylabelpos,'Interpreter', 'LaTeX');
% print('-depsc', '-r600','Ct_v_St_NoD.eps');




figure

set(gcf, 'Units', 'centimeters');
set(gcf, 'Position', afFigurePosition); % [left bottom width height]
set(gcf, 'PaperPositionMode', 'auto')
set(gcf,'DefaultAxesFontSize',FontSizeAx,'DefaultAxesFontName','TimesNewRoman','DefaultAxesGridLineStyle','-.','DefaultAxesLineWidth',2,'DefaultAxesFontWeight','Normal')

hold on
for i = 1:length(A_c_vec)
    beta = Beta(i);
%     plot(2*pi*f(:,i)*c/U, Pow_avg(:,i)./(rho*c*2*c*U^2*A(:,i).*f(:,i)), 's-','color',[(1 - colorvec(i)) 0 colorvec(i)], 'LineWidth', 2, 'MarkerSize',6)
    plot(2*pi*f(:,i)*c/U, Pow_avg(:,i)./(rho*c^2*beta*Sp*A(:,i).^2.*f(:,i).^3), 's-','color',[(1 - colorvec(i)) 0 colorvec(i)], 'LineWidth', 2, 'MarkerSize',6)

end
% plot([0.05 0.5], [0 0], '-k', 'LineWidth', 2, 'MarkerSize',6)
% plot([0.05 0.5], 0.6*[1 1], '--k', 'LineWidth', 2, 'MarkerSize',6)
hold off

% axis([0 40 0 40])
% axis([0 80 0 120])
set(gca, 'FontName', 'TimesNewRoman', 'FontSize', FontSizeAx)
set(gca, 'Units', 'normalized', 'Position', axespos);
xlabel('$$k$$','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 0,'Units', 'Normalize','Position', xlabelpos,'Interpreter', 'LaTeX');
ylabel('$$\mathcal{P}$$','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 90,'Units', 'Normalize','Position', ylabelpos,'Interpreter', 'LaTeX');
% print('-depsc', '-r600','Ct_v_St_NoD.eps');


% np = T_avg*U./Pow_avg;
% figure
% 
% set(gcf, 'Units', 'centimeters');
% set(gcf, 'Position', afFigurePosition); % [left bottom width height]
% set(gcf, 'PaperPositionMode', 'auto')
% set(gcf,'DefaultAxesFontSize',FontSizeAx,'DefaultAxesFontName','TimesNewRoman','DefaultAxesGridLineStyle','-.','DefaultAxesLineWidth',2,'DefaultAxesFontWeight','Normal')
% 
% hold on
% for i = 1:length(A_c_vec)-1
%     plot(2*pi*f(:,i)*c/U, np(:,i), 's-','color',[(1 - colorvec(i)) 0 colorvec(i)], 'LineWidth', 2, 'MarkerSize',6)
% end
% % plot([0.05 0.5], [0 0], '-k', 'LineWidth', 2, 'MarkerSize',6)
% % plot([0.05 0.5], 0.6*[1 1], '--k', 'LineWidth', 2, 'MarkerSize',6)
% hold off
% 
% axis([0 80 0 0.3])
% set(gca, 'FontName', 'TimesNewRoman', 'FontSize', FontSizeAx)
% set(gca, 'Units', 'normalized', 'Position', axespos);
% xlabel('$$k$$','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 0,'Units', 'Normalize','Position', xlabelpos,'Interpreter', 'LaTeX');
% ylabel('$$\eta$$','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 90,'Units', 'Normalize','Position', ylabelpos,'Interpreter', 'LaTeX');
% % print('-depsc', '-r600','Ct_v_St_NoD.eps');
% 
%% 
figure

set(gcf, 'Units', 'centimeters');
set(gcf, 'Position', afFigurePosition); % [left bottom width height]
set(gcf, 'PaperPositionMode', 'auto')
set(gcf,'DefaultAxesFontSize',FontSizeAx,'DefaultAxesFontName','TimesNewRoman','DefaultAxesGridLineStyle','-.','DefaultAxesLineWidth',2,'DefaultAxesFontWeight','Normal')

hold on
for i = 1:length(A_c_vec)
    beta = Beta(i);
%     plot(St_vec, np(:,i), 's-','color',[(1 - colorvec(i)) 0 colorvec(i)], 'LineWidth', 2, 'MarkerSize',6)
      plot(2*pi*f(:,i)*c/U, np(:,i), 's-','color',[(1 - colorvec(i)) 0 colorvec(i)], 'LineWidth', 2, 'MarkerSize',6)
end
% plot([0.05 0.5], [0 0], '-k', 'LineWidth', 2, 'MarkerSize',6)
% plot([0.05 0.5], 0.6*[1 1], '--k', 'LineWidth', 2, 'MarkerSize',6)
hold off

axis([0 40 0 0.5])
% axis([0.05 0.8 0 0.3])
set(gca, 'FontName', 'TimesNewRoman', 'FontSize', FontSizeAx)
set(gca, 'Units', 'normalized', 'Position', axespos);
xlabel('$$k$$','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 0,'Units', 'Normalize','Position', xlabelpos,'Interpreter', 'LaTeX');
ylabel('$$\eta$$','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 90,'Units', 'Normalize','Position', ylabelpos,'Interpreter', 'LaTeX');
% print('-depsc', '-r600','Ct_v_St_NoD.eps');




figure

set(gcf, 'Units', 'centimeters');
set(gcf, 'Position', afFigurePosition); % [left bottom width height]
set(gcf, 'PaperPositionMode', 'auto')
set(gcf,'DefaultAxesFontSize',FontSizeAx,'DefaultAxesFontName','TimesNewRoman','DefaultAxesGridLineStyle','-.','DefaultAxesLineWidth',2,'DefaultAxesFontWeight','Normal')

hold on
for i = 1:length(A_c_vec)
    beta = Beta(i);
%     plot(St_vec, np(:,i), 's-','color',[(1 - colorvec(i)) 0 colorvec(i)], 'LineWidth', 2, 'MarkerSize',6)
      plot(2*pi*f(:,i)*c/U, np_ND(:,i), 's-','color',[(1 - colorvec(i)) 0 colorvec(i)], 'LineWidth', 2, 'MarkerSize',6)
end
% plot([0.05 0.5], [0 0], '-k', 'LineWidth', 2, 'MarkerSize',6)
% plot([0.05 0.5], 0.6*[1 1], '--k', 'LineWidth', 2, 'MarkerSize',6)
hold off

axis([0 40 0 0.5])
% axis([0.05 0.8 0 0.3])
set(gca, 'FontName', 'TimesNewRoman', 'FontSize', FontSizeAx)
set(gca, 'Units', 'normalized', 'Position', axespos);
xlabel('$$k$$','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 0,'Units', 'Normalize','Position', xlabelpos,'Interpreter', 'LaTeX');
ylabel('$$\eta$$','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 90,'Units', 'Normalize','Position', ylabelpos,'Interpreter', 'LaTeX');
% print('-depsc', '-r600','Ct_v_St_NoD.eps');



figure

set(gcf, 'Units', 'centimeters');
set(gcf, 'Position', afFigurePosition); % [left bottom width height]
set(gcf, 'PaperPositionMode', 'auto')
set(gcf,'DefaultAxesFontSize',FontSizeAx,'DefaultAxesFontName','TimesNewRoman','DefaultAxesGridLineStyle','-.','DefaultAxesLineWidth',2,'DefaultAxesFontWeight','Normal')

hold on
for i = 1:length(A_c_vec)
    beta = Beta(i);
%     plot(St_vec, np(:,i), 's-','color',[(1 - colorvec(i)) 0 colorvec(i)], 'LineWidth', 2, 'MarkerSize',6)
      plot(2*pi*f(:,i)*c/U, np_ND(:,i)*beta*c.*f(:,i)/U, 's-','color',[(1 - colorvec(i)) 0 colorvec(i)], 'LineWidth', 2, 'MarkerSize',6)
end
% plot([0.05 0.5], [0 0], '-k', 'LineWidth', 2, 'MarkerSize',6)
% plot([0.05 0.5], 0.6*[1 1], '--k', 'LineWidth', 2, 'MarkerSize',6)
hold off

axis([0 40 0 0.5])
% axis([0.05 0.8 0 0.3])
set(gca, 'FontName', 'TimesNewRoman', 'FontSize', FontSizeAx)
set(gca, 'Units', 'normalized', 'Position', axespos);
xlabel('$$k$$','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 0,'Units', 'Normalize','Position', xlabelpos,'Interpreter', 'LaTeX');
ylabel('$$\eta$$','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 90,'Units', 'Normalize','Position', ylabelpos,'Interpreter', 'LaTeX');
% print('-depsc', '-r600','Ct_v_St_NoD.eps');


%% 
% figure
% 
% set(gcf, 'Units', 'centimeters');
% set(gcf, 'Position', afFigurePosition); % [left bottom width height]
% set(gcf, 'PaperPositionMode', 'auto')
% set(gcf,'DefaultAxesFontSize',FontSizeAx,'DefaultAxesFontName','TimesNewRoman','DefaultAxesGridLineStyle','-.','DefaultAxesLineWidth',2,'DefaultAxesFontWeight','Normal')
% 
% hold on
% for i = 1:length(A_c_vec)
%     plot(St_vec, Pow_avg(:,i)./(rho*c*A(:,i).^3.*f(:,i).^3), 's-','color',[(1 - colorvec(i)) 0 colorvec(i)], 'LineWidth', 2, 'MarkerSize',6)
% end
% % plot([0.05^3 0.5^3], [0 0], '-k', 'LineWidth', 2, 'MarkerSize',6)
% hold off
% 
% % axis([0 2 0 0.4])
% set(gca, 'FontName', 'TimesNewRoman', 'FontSize', FontSizeAx)
% set(gca, 'Units', 'normalized', 'Position', axespos);
% xlabel('$$St$$','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 0,'Units', 'Normalize','Position', xlabelpos,'Interpreter', 'LaTeX');
% ylabel('$$C_P$$','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 90,'Units', 'Normalize','Position', ylabelpos,'Interpreter', 'LaTeX');
% 
% figure
% 
% set(gcf, 'Units', 'centimeters');
% set(gcf, 'Position', afFigurePosition); % [left bottom width height]
% set(gcf, 'PaperPositionMode', 'auto')
% set(gcf,'DefaultAxesFontSize',FontSizeAx,'DefaultAxesFontName','TimesNewRoman','DefaultAxesGridLineStyle','-.','DefaultAxesLineWidth',2,'DefaultAxesFontWeight','Normal')
% 
% hold on
% for i = 1:length(A_c_vec)
%     plot(St_vec, Pow_avg(:,i), 's-','color',[(1 - colorvec(i)) 0 colorvec(i)], 'LineWidth', 2, 'MarkerSize',6)
% end
% % plot([0.05^3 0.5^3], [0 0], '-k', 'LineWidth', 2, 'MarkerSize',6)
% hold off
% 
% % axis([0 2 0 0.4])
% set(gca, 'FontName', 'TimesNewRoman', 'FontSize', FontSizeAx)
% set(gca, 'Units', 'normalized', 'Position', axespos);
% xlabel('$$St$$','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 0,'Units', 'Normalize','Position', xlabelpos,'Interpreter', 'LaTeX');
% ylabel('$$C_P$$','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 90,'Units', 'Normalize','Position', ylabelpos,'Interpreter', 'LaTeX');
% 
% 
% figure
% 
% set(gcf, 'Units', 'centimeters');
% set(gcf, 'Position', afFigurePosition); % [left bottom width height]
% set(gcf, 'PaperPositionMode', 'auto')
% set(gcf,'DefaultAxesFontSize',FontSizeAx,'DefaultAxesFontName','TimesNewRoman','DefaultAxesGridLineStyle','-.','DefaultAxesLineWidth',2,'DefaultAxesFontWeight','Normal')
% 
% hold on
% for i = 1:length(A_c_vec)
%     plot(St_vec, Pow_avg(:,i)./(rho*c^2*f(:,i).^3.*A(:,i).^3), 's-','color',[(1 - colorvec(i)) 0 colorvec(i)], 'LineWidth', 2, 'MarkerSize',6)
% end
% % plot([0.05^3 0.5^3], [0 0], '-k', 'LineWidth', 2, 'MarkerSize',6)
% hold off
% 
% % axis([0 2 0 0.4])
% set(gca, 'FontName', 'TimesNewRoman', 'FontSize', FontSizeAx)
% set(gca, 'Units', 'normalized', 'Position', axespos);
% xlabel('$$St$$','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 0,'Units', 'Normalize','Position', xlabelpos,'Interpreter', 'LaTeX');
% ylabel('$$C_P$$','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 90,'Units', 'Normalize','Position', ylabelpos,'Interpreter', 'LaTeX');
% 

