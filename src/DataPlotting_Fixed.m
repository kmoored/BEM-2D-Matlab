clear
close all
clc


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

num_vec = 1:6;
colorvec = -((num_vec-1)/5 - 1);

T_avg = Tnet_avg + D_avg;
% np_ND = T_avg*U./Pow_avg;

% Figure 1: Classic Thrust Coefficient
% Figure 1: Classic Thrust Coefficient
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
            end
            f_vec = k_vec*Qinf/(2*pi*c);
            St_vec = f_vec*(A_c*c)/Qinf;
            plot(k_vec, Tnet_avg(:,i,j)./(1/2*rho*Sp*Qinf^2), vec,'color',[(1 - colorvec(i)) 0 colorvec(i)], 'LineWidth', 2, 'MarkerSize',6)
%             plot((rho*Sp*Qinf^2), T_avg(:,i,j), vec,'color',[(1 - colorvec(i)) 0 colorvec(i)], 'LineWidth', 2, 'MarkerSize',6)

        end
    end
plot([0 20], [0 0], '-k', 'LineWidth', 1, 'MarkerSize',6)

hold off

% axis([0 0.09 -0.01 0.125])
axis([0 15 -0.01 10])
set(gca, 'FontName', 'TimesNewRoman', 'FontSize', FontSizeAx)
set(gca, 'Units', 'normalized', 'Position', axespos);
xlabel('$$k$$','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 0,'Units', 'Normalize','Position', xlabelpos,'Interpreter', 'LaTeX');
ylabel('$$C_t$$','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 90,'Units', 'Normalize','Position', ylabelpos,'Interpreter', 'LaTeX');
% print('-depsc', '-r600','T_v_f.eps');

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
            end
            f_vec = k_vec*Qinf/(2*pi*c);
            St_vec = f_vec*(A_c*c)/Qinf;
            plot(St_vec, Tnet_avg(:,i,j)./(1/2*rho*Sp*Qinf^2), vec,'color',[(1 - colorvec(i)) 0 colorvec(i)], 'LineWidth', 2, 'MarkerSize',6)
%             plot((rho*Sp*Qinf^2), T_avg(:,i,j), vec,'color',[(1 - colorvec(i)) 0 colorvec(i)], 'LineWidth', 2, 'MarkerSize',6)

        end
    end
plot([0 20], [0 0], '-k', 'LineWidth', 1, 'MarkerSize',6)

hold off

% axis([0 0.09 -0.01 0.125])
axis([0 1.25 -0.01 10])
set(gca, 'FontName', 'TimesNewRoman', 'FontSize', FontSizeAx)
set(gca, 'Units', 'normalized', 'Position', axespos);
xlabel('$$St$$','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 0,'Units', 'Normalize','Position', xlabelpos,'Interpreter', 'LaTeX');
ylabel('$$C_t$$','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 90,'Units', 'Normalize','Position', ylabelpos,'Interpreter', 'LaTeX');
% print('-depsc', '-r600','T_v_f.eps');



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
            end
            f_vec = k_vec*Qinf/(2*pi*c);
            St_vec = f_vec*(A_c*c)/Qinf;
            plot(k_vec, Pow_avg(:,i,j)./(1/2*rho*Sp*Qinf^3), vec,'color',[(1 - colorvec(i)) 0 colorvec(i)], 'LineWidth', 2, 'MarkerSize',6)
%             plot((rho*Sp*Qinf^2), T_avg(:,i,j), vec,'color',[(1 - colorvec(i)) 0 colorvec(i)], 'LineWidth', 2, 'MarkerSize',6)

        end
    end
plot([0 20], [0 0], '-k', 'LineWidth', 1, 'MarkerSize',6)

hold off

% axis([0 0.09 -0.01 0.125])
axis([0 15 -0.01 30])
set(gca, 'FontName', 'TimesNewRoman', 'FontSize', FontSizeAx)
set(gca, 'Units', 'normalized', 'Position', axespos);
xlabel('$$k$$','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 0,'Units', 'Normalize','Position', xlabelpos,'Interpreter', 'LaTeX');
ylabel('$$C_p$$','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 90,'Units', 'Normalize','Position', ylabelpos,'Interpreter', 'LaTeX');
% print('-depsc', '-r600','T_v_f.eps');


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
            end
            f_vec = k_vec*Qinf/(2*pi*c);
            St_vec = f_vec*(A_c*c)/Qinf;
            plot(St_vec.^3, Pow_avg(:,i,j)./(1/2*rho*Sp*Qinf^3), vec,'color',[(1 - colorvec(i)) 0 colorvec(i)], 'LineWidth', 2, 'MarkerSize',6)
%             plot((rho*Sp*Qinf^2), T_avg(:,i,j), vec,'color',[(1 - colorvec(i)) 0 colorvec(i)], 'LineWidth', 2, 'MarkerSize',6)

        end
    end
plot([0 20], [0 0], '-k', 'LineWidth', 1, 'MarkerSize',6)

hold off

% axis([0 0.09 -0.01 0.125])
axis([0 1. -0.01 30])
set(gca, 'FontName', 'TimesNewRoman', 'FontSize', FontSizeAx)
set(gca, 'Units', 'normalized', 'Position', axespos);
xlabel('$$St$$','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 0,'Units', 'Normalize','Position', xlabelpos,'Interpreter', 'LaTeX');
ylabel('$$C_p$$','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 90,'Units', 'Normalize','Position', ylabelpos,'Interpreter', 'LaTeX');
% print('-depsc', '-r600','T_v_f.eps');













%%
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
            end
            f_vec = k_vec*Qinf/(2*pi*c);
            St_vec = f_vec*(A_c*c)/Qinf;
            plot((rho*Beta*Sp*f_vec.^2*A_c.^2*c.^2), T_avg(:,i,j), vec,'color',[(1 - colorvec(i)) 0 colorvec(i)], 'LineWidth', 2, 'MarkerSize',6)
%             plot((rho*Sp*Qinf^2), T_avg(:,i,j), vec,'color',[(1 - colorvec(i)) 0 colorvec(i)], 'LineWidth', 2, 'MarkerSize',6)

        end
    end

    plot([0 6e4], 2.1*[0 6e4], '-k', 'LineWidth', 2, 'MarkerSize',6)
hold off

% axis([0 0.09 -0.01 0.125])
axis([0 6e4 -0.01 2e5])
set(gca, 'FontName', 'TimesNewRoman', 'FontSize', FontSizeAx)
set(gca, 'Units', 'normalized', 'Position', axespos);
xlabel('$$\rho S_p f^2 A^2$$','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 0,'Units', 'Normalize','Position', xlabelpos,'Interpreter', 'LaTeX');
ylabel('$$T$$','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 90,'Units', 'Normalize','Position', ylabelpos,'Interpreter', 'LaTeX');
% print('-depsc', '-r600','T_v_f.eps');


% Figure 2: New Thrust Coefficient
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
        end
        
        f_vec = (Qinf/2/pi/c)*k_vec;
        St_vec = (A_c*c/Qinf)*f_vec;
%         Beta = St_vec;
        plot(1/2/pi*k_vec, T_avg(:,i,j)./(rho*Beta*Sp*f_vec.^2*(A_c*c).^2), vec,'color',[(1 - colorvec(i)) 0 colorvec(i)], 'LineWidth', 2, 'MarkerSize',6)
    end
end
% plot([0.0007 0.08], [0 0], '-k', 'LineWidth', 2, 'MarkerSize',6)
hold off

% axis([0 0.09 -0.01 0.125])
axis([0 3 0 3])
set(gca, 'FontName', 'TimesNewRoman', 'FontSize', FontSizeAx)
set(gca, 'Units', 'normalized', 'Position', axespos);
xlabel('$$k$$','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 0,'Units', 'Normalize','Position', xlabelpos,'Interpreter', 'LaTeX');
ylabel('$$C_T$$','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 90,'Units', 'Normalize','Position', ylabelpos,'Interpreter', 'LaTeX');
% print('-depsc', '-r600','T_v_f.eps');


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
        end
        
        f_vec = (Qinf/2/pi/c)*k_vec;
        St_vec = (A_c*c/Qinf)*f_vec;
%         Beta = 1/(2*asin(A_c/2));
        plot(St_vec, Pow_avg(:,i,j)./(rho*Beta*Sp*Qinf*f_vec.^2.*(A_c*c).^2), vec,'color',[(1 - colorvec(i)) 0 colorvec(i)], 'LineWidth', 2, 'MarkerSize',6)
    end
end
% plot([0.0007 0.08], [0 0], '-k', 'LineWidth', 2, 'MarkerSize',6)
hold off

% axis([0 0.09 -0.01 0.125])
axis([0 1.25 0 100])
set(gca, 'FontName', 'TimesNewRoman', 'FontSize', FontSizeAx)
set(gca, 'Units', 'normalized', 'Position', axespos);
xlabel('$$St$$','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 0,'Units', 'Normalize','Position', xlabelpos,'Interpreter', 'LaTeX');
ylabel('$$C_P$$','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 90,'Units', 'Normalize','Position', ylabelpos,'Interpreter', 'LaTeX');
% print('-depsc', '-r600','T_v_f.eps');




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
        end
        
        f_vec = (Qinf/2/pi/c)*k_vec;
        St_vec = (A_c*c/Qinf)*f_vec;
%         Beta = 1/(2*asin(A_c/2));
        kappa = 7;
        plot((rho*Sp.*f_vec.^3.*(A_c*c).^3), Pow_avg(:,i,j), vec,'color',[(1 - colorvec(i)) 0 colorvec(i)], 'LineWidth', 2, 'MarkerSize',6)
    end
end
    plot([0 3e5], 18*[0 3e5], '-k', 'LineWidth', 2, 'MarkerSize',6)
hold off

% axis([0 0.09 -0.01 0.125])
axis([0 1e5 0 5e6])
set(gca, 'FontName', 'TimesNewRoman', 'FontSize', FontSizeAx)
set(gca, 'Units', 'normalized', 'Position', axespos);
xlabel('$$\rho S_p f^3 A^3$$','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 0,'Units', 'Normalize','Position', xlabelpos,'Interpreter', 'LaTeX');
ylabel('$$P$$','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 90,'Units', 'Normalize','Position', ylabelpos,'Interpreter', 'LaTeX');
%




% Thrust verus power
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
        end
        
        f_vec = (Qinf/2/pi/c)*k_vec;
        St_vec = (A_c*c/Qinf)*f_vec;
%         Beta = St_vec;
        plot(1/2/pi*k_vec,Tnet_avg(:,i,j)*Qinf./Pow_avg(:,i,j)*sqrt(A_c), vec,'color',[(1 - colorvec(i)) 0 colorvec(i)], 'LineWidth', 2, 'MarkerSize',6)
    end
end
% plot([0.0007 0.08], [0 0], '-k', 'LineWidth', 2, 'MarkerSize',6)
hold off

% axis([0 0.09 -0.01 0.125])
axis([0 1.25 0 0.4])
set(gca, 'FontName', 'TimesNewRoman', 'FontSize', FontSizeAx)
set(gca, 'Units', 'normalized', 'Position', axespos);
xlabel('$$k$$','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 0,'Units', 'Normalize','Position', xlabelpos,'Interpreter', 'LaTeX');
ylabel('$$\eta$$','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 90,'Units', 'Normalize','Position', ylabelpos,'Interpreter', 'LaTeX');
% print('-depsc', '-r600','T_v_f.eps');

% Thrust verus power
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
        end
        
        f_vec = (Qinf/2/pi/c)*k_vec;
        St_vec = (A_c*c/Qinf)*f_vec;
%         Beta = St_vec;
        plot(1/2/pi*k_vec,T_avg(:,i,j)*Qinf./Pow_avg(:,i,j), vec,'color',[(1 - colorvec(i)) 0 colorvec(i)], 'LineWidth', 2, 'MarkerSize',6)
    end
end
% plot([0.0007 0.08], [0 0], '-k', 'LineWidth', 2, 'MarkerSize',6)
hold off

% axis([0 0.09 -0.01 0.125])
axis([0 1.25 0 0.4])
set(gca, 'FontName', 'TimesNewRoman', 'FontSize', FontSizeAx)
set(gca, 'Units', 'normalized', 'Position', axespos);
xlabel('$$k$$','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 0,'Units', 'Normalize','Position', xlabelpos,'Interpreter', 'LaTeX');
ylabel('$$\eta_I$$','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 90,'Units', 'Normalize','Position', ylabelpos,'Interpreter', 'LaTeX');
% print('-depsc', '-r600','T_v_f.eps');


% % Thrust verus power
% figure
% 
% set(gcf, 'Units', 'centimeters');
% set(gcf, 'Position', afFigurePosition); % [left bottom width height]
% set(gcf, 'PaperPositionMode', 'auto')
% set(gcf,'DefaultAxesFontSize',FontSizeAx,'DefaultAxesFontName','TimesNewRoman','DefaultAxesGridLineStyle','-.','DefaultAxesLineWidth',2,'DefaultAxesFontWeight','Normal')
% 
% hold on
% for j = 1:length(Qinf_vec)
%     Qinf = Qinf_vec(j);
%     for i = 1:length(A_c_vec)
%         A_c = A_c_vec(i);
%         if j == 1
%             vec = 's-';
%         elseif j == 2
%             vec = 's--';
%         elseif j == 3
%             vec = 's:';
%         elseif j == 4
%             vec = 's-.';
%         else 
%             vec = 'o-';
%         end
%         
%         f_vec = (Qinf/2/pi/c)*k_vec;
%         St_vec = (A_c*c/Qinf)*f_vec;
% %         Beta = St_vec;
%         plot(k_vec,T_avg(:,i,j).*(f_vec*c*sqrt(A_c))./Pow_avg(:,i,j), vec,'color',[(1 - colorvec(i)) 0 colorvec(i)], 'LineWidth', 2, 'MarkerSize',6)
%     end
% end
% % plot([0.0007 0.08], [0 0], '-k', 'LineWidth', 2, 'MarkerSize',6)
% hold off
% 
% % axis([0 0.09 -0.01 0.125])
% axis([0 20 0 0.4])
% set(gca, 'FontName', 'TimesNewRoman', 'FontSize', FontSizeAx)
% set(gca, 'Units', 'normalized', 'Position', axespos);
% xlabel('$$St$$','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 0,'Units', 'Normalize','Position', xlabelpos,'Interpreter', 'LaTeX');
% ylabel('$$P$$','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 90,'Units', 'Normalize','Position', ylabelpos,'Interpreter', 'LaTeX');
% % print('-depsc', '-r600','T_v_f.eps');



% % Drag Coefficient
% figure
% 
% set(gcf, 'Units', 'centimeters');
% set(gcf, 'Position', afFigurePosition); % [left bottom width height]
% set(gcf, 'PaperPositionMode', 'auto')
% set(gcf,'DefaultAxesFontSize',FontSizeAx,'DefaultAxesFontName','TimesNewRoman','DefaultAxesGridLineStyle','-.','DefaultAxesLineWidth',2,'DefaultAxesFontWeight','Normal')
% 
% hold on
% for j = 1:length(Qinf_vec)
%     Qinf = Qinf_vec(j);
%     for i = 1:length(A_c_vec)
%         A_c = A_c_vec(i);
%         if j == 1
%             vec = 's-';
%         elseif j == 2
%             vec = 's--';
%         elseif j == 3
%             vec = 's:';
%         elseif j == 4
%             vec = 's-.';
%         else 
%             vec = 'o-';
%         end
%         f_vec = (Qinf/2/pi/c)*k_vec;
%         St_vec = (A_c*c/Qinf)*f_vec;
%         AoA = atan2(f_vec.*A_c*c,Qinf);
%         plot(AoA, D_avg(:,i,j)./(1/2*rho*Sp*Qinf^2), vec,'color',[(1 - colorvec(i)) 0 colorvec(i)], 'LineWidth', 2, 'MarkerSize',6)
%     end
% end
% % plot([0.0007 0.08], [0 0], '-k', 'LineWidth', 2, 'MarkerSize',6)
% hold off
% 
% % axis([0 0.09 -0.01 0.125])
% axis([0 0.8 0 0.05])
% set(gca, 'FontName', 'TimesNewRoman', 'FontSize', FontSizeAx)
% set(gca, 'Units', 'normalized', 'Position', axespos);
% xlabel('$$k$$','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 0,'Units', 'Normalize','Position', xlabelpos,'Interpreter', 'LaTeX');
% ylabel('$$C_d$$','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 90,'Units', 'Normalize','Position', ylabelpos,'Interpreter', 'LaTeX');
% % print('-depsc', '-r600','T_v_f.eps');


% figure
% 
% set(gcf, 'Units', 'centimeters');
% set(gcf, 'Position', afFigurePosition); % [left bottom width height]
% set(gcf, 'PaperPositionMode', 'auto')
% set(gcf,'DefaultAxesFontSize',FontSizeAx,'DefaultAxesFontName','TimesNewRoman','DefaultAxesGridLineStyle','-.','DefaultAxesLineWidth',2,'DefaultAxesFontWeight','Normal')
% 
% hold on
% for i = 1:length(A_c_vec)
%     plot(St_vec, T_avg(:,i), 's-','color',[(1 - colorvec(i)) 0 colorvec(i)], 'LineWidth', 2, 'MarkerSize',6)
% end
% % plot([0.05 0.5], [0 0], '-k', 'LineWidth', 2, 'MarkerSize',6)
% hold off
% 
% % axis([0.05 0.5 -0.01 0.2])
% set(gca, 'FontName', 'TimesNewRoman', 'FontSize', FontSizeAx)
% set(gca, 'Units', 'normalized', 'Position', axespos);
% xlabel('$$St$$','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 0,'Units', 'Normalize','Position', xlabelpos,'Interpreter', 'LaTeX');
% ylabel('$$T, N$$','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 90,'Units', 'Normalize','Position', ylabelpos,'Interpreter', 'LaTeX');
% % print('-depsc', '-r600','T_v_St.eps');
% 
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
%     beta = Beta(i);
%     if i == 1 
%         plot(2*pi*f(1:8,1)*c/U, Tnet_avg(1:8,i)./(rho*c*Sp*A(1:8,i).^2.*f(1:8,1).^2), 's-','color',[(1 - colorvec(i)) 0 colorvec(i)], 'LineWidth', 2, 'MarkerSize',6)
%     else
%         plot(2*pi*f(:,i)*c/U, Tnet_avg(:,i)./(rho*c*Sp*A(:,i).^2.*f(:,i).^2), 's-','color',[(1 - colorvec(i)) 0 colorvec(i)], 'LineWidth', 2, 'MarkerSize',6)
%     end
% end
% % plot([0.05 0.5], [0 0], '-k', 'LineWidth', 2, 'MarkerSize',6)
% plot([0 41], 2.25*[1 1], '--k', 'LineWidth', 2, 'MarkerSize',6)
% hold off
% 
% axis([0 41 0 4])
% % axis([0.05 0.5 0 4])
% set(gca, 'FontName', 'TimesNewRoman', 'FontSize', FontSizeAx)
% set(gca, 'Units', 'normalized', 'Position', axespos);
% xlabel('$$k$$','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 0,'Units', 'Normalize','Position', xlabelpos,'Interpreter', 'LaTeX');
% ylabel('$$C_T$$','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 90,'Units', 'Normalize','Position', ylabelpos,'Interpreter', 'LaTeX');
% % print('-depsc', '-r600','Ct_v_k.eps');
% 
% 
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
% for i = length(A_c_vec):-1:1
%     beta = Beta(i);
%     if i == 1 
%             plot(2*pi*f(1:8,1)*c/U, T_avg(1:8,i)./(rho*c*Sp*A(1:8,i).^2.*f(1:8,1).^2), 's-','color',[(1 - colorvec(i)) 0 colorvec(i)], 'LineWidth', 2, 'MarkerSize',6)
%     else
%         plot(2*pi*f(:,i)*c/U, T_avg(:,i)./(rho*c*Sp*A(:,i).^2.*f(:,i).^2), 's-','color',[(1 - colorvec(i)) 0 colorvec(i)], 'LineWidth', 2, 'MarkerSize',6)
%     end
% %     plot(St_vec, T_avg(:,i)./(rho*c*2*c*A(:,i).^2.*f(:,i).^2), 's-','color',[(1 - colorvec(i)) 0 colorvec(i)], 'LineWidth', 2, 'MarkerSize',6)
% end
% plot([0 41], 2.25*[1 1], '--k', 'LineWidth', 2, 'MarkerSize',6)
% 
% % plot([0 50], 2.25*[1 1], '--k', 'LineWidth', 2, 'MarkerSize',6)
% 
% hold off
% 
% axis([0 40 0 4])
% % axis([0.05 0.5 0 3])
% set(gca, 'FontName', 'TimesNewRoman', 'FontSize', FontSizeAx)
% set(gca, 'Units', 'normalized', 'Position', axespos);
% xlabel('$$k$$','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 0,'Units', 'Normalize','Position', xlabelpos,'Interpreter', 'LaTeX');
% ylabel('$$C_T$$','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 90,'Units', 'Normalize','Position', ylabelpos,'Interpreter', 'LaTeX');
% % print('-depsc', '-r600','Ct_v_k_NoD.eps');
% 
% 
% 
% % figure
% % 
% % set(gcf, 'Units', 'centimeters');
% % set(gcf, 'Position', afFigurePosition); % [left bottom width height]
% % set(gcf, 'PaperPositionMode', 'auto')
% % set(gcf,'DefaultAxesFontSize',FontSizeAx,'DefaultAxesFontName','TimesNewRoman','DefaultAxesGridLineStyle','-.','DefaultAxesLineWidth',2,'DefaultAxesFontWeight','Normal')
% % 
% % hold on
% % for i = length(A_c_vec)-1:-1:1
% % %     if i == 1 
% % %             plot(2*pi*f(1:8,1)*c/U, T_avg(1:8,i)./(rho*c*2*c*A(1:8,i).^2.*f(1:8,1).^2), 's-','color',[(1 - colorvec(i)) 0 colorvec(i)], 'LineWidth', 2, 'MarkerSize',6)
% % %     else
% % %         plot(2*pi*f(:,i)*c/U, T_avg(:,i)./(rho*c*2*c*A(:,i).^2.*f(:,i).^2), 's-','color',[(1 - colorvec(i)) 0 colorvec(i)], 'LineWidth', 2, 'MarkerSize',6)
% % %     end
% %     plot(St_vec, T_avg(:,i)./(rho*c*2*c*A(:,i).^2.*f(:,i).^2), 's-','color',[(1 - colorvec(i)) 0 colorvec(i)], 'LineWidth', 2, 'MarkerSize',6)
% % end
% % plot([0 41], 2.25*[1 1], '--k', 'LineWidth', 2, 'MarkerSize',6)
% % 
% % % plot([0 50], 2.25*[1 1], '--k', 'LineWidth', 2, 'MarkerSize',6)
% % 
% % hold off
% % 
% % % axis([0 41 0 3])
% % axis([0.05 0.8 0 3])
% % set(gca, 'FontName', 'TimesNewRoman', 'FontSize', FontSizeAx)
% % set(gca, 'Units', 'normalized', 'Position', axespos);
% % xlabel('$$k$$','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 0,'Units', 'Normalize','Position', xlabelpos,'Interpreter', 'LaTeX');
% % ylabel('$$C_T$$','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 90,'Units', 'Normalize','Position', ylabelpos,'Interpreter', 'LaTeX');
% % % print('-depsc', '-r600','Ct_v_k_NoD.eps');
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
%     beta = Beta(i);
% %     plot(2*pi*f(:,i)*c/U, Pow_avg(:,i)./(rho*c*2*c*U^2*A(:,i).*f(:,i)), 's-','color',[(1 - colorvec(i)) 0 colorvec(i)], 'LineWidth', 2, 'MarkerSize',6)
%     plot(2*pi*f(:,i)*c/U, Pow_avg(:,i), 's-','color',[(1 - colorvec(i)) 0 colorvec(i)], 'LineWidth', 2, 'MarkerSize',6)
% 
% end
% % plot([0.05 0.5], [0 0], '-k', 'LineWidth', 2, 'MarkerSize',6)
% % plot([0.05 0.5], 0.6*[1 1], '--k', 'LineWidth', 2, 'MarkerSize',6)
% hold off
% 
% % axis([0 40 0 40])
% % axis([0 80 0 120])
% set(gca, 'FontName', 'TimesNewRoman', 'FontSize', FontSizeAx)
% set(gca, 'Units', 'normalized', 'Position', axespos);
% xlabel('$$k$$','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 0,'Units', 'Normalize','Position', xlabelpos,'Interpreter', 'LaTeX');
% ylabel('$$\mathcal{P}$$','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 90,'Units', 'Normalize','Position', ylabelpos,'Interpreter', 'LaTeX');
% % print('-depsc', '-r600','Ct_v_St_NoD.eps');
% 
% 
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
%     beta = Beta(i);
% %     plot(2*pi*f(:,i)*c/U, Pow_avg(:,i)./(rho*c*2*c*U^2*A(:,i).*f(:,i)), 's-','color',[(1 - colorvec(i)) 0 colorvec(i)], 'LineWidth', 2, 'MarkerSize',6)
%     plot(2*pi*f(:,i)*c/U, Pow_avg(:,i)./(rho*c^2*beta*Sp*A(:,i).^2.*f(:,i).^3), 's-','color',[(1 - colorvec(i)) 0 colorvec(i)], 'LineWidth', 2, 'MarkerSize',6)
% 
% end
% % plot([0.05 0.5], [0 0], '-k', 'LineWidth', 2, 'MarkerSize',6)
% % plot([0.05 0.5], 0.6*[1 1], '--k', 'LineWidth', 2, 'MarkerSize',6)
% hold off
% 
% % axis([0 40 0 40])
% % axis([0 80 0 120])
% set(gca, 'FontName', 'TimesNewRoman', 'FontSize', FontSizeAx)
% set(gca, 'Units', 'normalized', 'Position', axespos);
% xlabel('$$k$$','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 0,'Units', 'Normalize','Position', xlabelpos,'Interpreter', 'LaTeX');
% ylabel('$$\mathcal{P}$$','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 90,'Units', 'Normalize','Position', ylabelpos,'Interpreter', 'LaTeX');
% % print('-depsc', '-r600','Ct_v_St_NoD.eps');
% 
% 
% % np = T_avg*U./Pow_avg;
% % figure
% % 
% % set(gcf, 'Units', 'centimeters');
% % set(gcf, 'Position', afFigurePosition); % [left bottom width height]
% % set(gcf, 'PaperPositionMode', 'auto')
% % set(gcf,'DefaultAxesFontSize',FontSizeAx,'DefaultAxesFontName','TimesNewRoman','DefaultAxesGridLineStyle','-.','DefaultAxesLineWidth',2,'DefaultAxesFontWeight','Normal')
% % 
% % hold on
% % for i = 1:length(A_c_vec)-1
% %     plot(2*pi*f(:,i)*c/U, np(:,i), 's-','color',[(1 - colorvec(i)) 0 colorvec(i)], 'LineWidth', 2, 'MarkerSize',6)
% % end
% % % plot([0.05 0.5], [0 0], '-k', 'LineWidth', 2, 'MarkerSize',6)
% % % plot([0.05 0.5], 0.6*[1 1], '--k', 'LineWidth', 2, 'MarkerSize',6)
% % hold off
% % 
% % axis([0 80 0 0.3])
% % set(gca, 'FontName', 'TimesNewRoman', 'FontSize', FontSizeAx)
% % set(gca, 'Units', 'normalized', 'Position', axespos);
% % xlabel('$$k$$','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 0,'Units', 'Normalize','Position', xlabelpos,'Interpreter', 'LaTeX');
% % ylabel('$$\eta$$','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 90,'Units', 'Normalize','Position', ylabelpos,'Interpreter', 'LaTeX');
% % % print('-depsc', '-r600','Ct_v_St_NoD.eps');
% % 
% %% 
% figure
% 
% set(gcf, 'Units', 'centimeters');
% set(gcf, 'Position', afFigurePosition); % [left bottom width height]
% set(gcf, 'PaperPositionMode', 'auto')
% set(gcf,'DefaultAxesFontSize',FontSizeAx,'DefaultAxesFontName','TimesNewRoman','DefaultAxesGridLineStyle','-.','DefaultAxesLineWidth',2,'DefaultAxesFontWeight','Normal')
% 
% hold on
% for i = 1:length(A_c_vec)
%     beta = Beta(i);
% %     plot(St_vec, np(:,i), 's-','color',[(1 - colorvec(i)) 0 colorvec(i)], 'LineWidth', 2, 'MarkerSize',6)
%       plot(2*pi*f(:,i)*c/U, np(:,i), 's-','color',[(1 - colorvec(i)) 0 colorvec(i)], 'LineWidth', 2, 'MarkerSize',6)
% end
% % plot([0.05 0.5], [0 0], '-k', 'LineWidth', 2, 'MarkerSize',6)
% % plot([0.05 0.5], 0.6*[1 1], '--k', 'LineWidth', 2, 'MarkerSize',6)
% hold off
% 
% axis([0 40 0 0.5])
% % axis([0.05 0.8 0 0.3])
% set(gca, 'FontName', 'TimesNewRoman', 'FontSize', FontSizeAx)
% set(gca, 'Units', 'normalized', 'Position', axespos);
% xlabel('$$k$$','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 0,'Units', 'Normalize','Position', xlabelpos,'Interpreter', 'LaTeX');
% ylabel('$$\eta$$','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 90,'Units', 'Normalize','Position', ylabelpos,'Interpreter', 'LaTeX');
% % print('-depsc', '-r600','Ct_v_St_NoD.eps');
% 
% 
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
%     beta = Beta(i);
% %     plot(St_vec, np(:,i), 's-','color',[(1 - colorvec(i)) 0 colorvec(i)], 'LineWidth', 2, 'MarkerSize',6)
%       plot(2*pi*f(:,i)*c/U, np_ND(:,i), 's-','color',[(1 - colorvec(i)) 0 colorvec(i)], 'LineWidth', 2, 'MarkerSize',6)
% end
% % plot([0.05 0.5], [0 0], '-k', 'LineWidth', 2, 'MarkerSize',6)
% % plot([0.05 0.5], 0.6*[1 1], '--k', 'LineWidth', 2, 'MarkerSize',6)
% hold off
% 
% axis([0 40 0 0.5])
% % axis([0.05 0.8 0 0.3])
% set(gca, 'FontName', 'TimesNewRoman', 'FontSize', FontSizeAx)
% set(gca, 'Units', 'normalized', 'Position', axespos);
% xlabel('$$k$$','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 0,'Units', 'Normalize','Position', xlabelpos,'Interpreter', 'LaTeX');
% ylabel('$$\eta$$','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 90,'Units', 'Normalize','Position', ylabelpos,'Interpreter', 'LaTeX');
% % print('-depsc', '-r600','Ct_v_St_NoD.eps');
% 
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
%     beta = Beta(i);
% %     plot(St_vec, np(:,i), 's-','color',[(1 - colorvec(i)) 0 colorvec(i)], 'LineWidth', 2, 'MarkerSize',6)
%       plot(2*pi*f(:,i)*c/U, np_ND(:,i)*beta*c.*f(:,i)/U, 's-','color',[(1 - colorvec(i)) 0 colorvec(i)], 'LineWidth', 2, 'MarkerSize',6)
% end
% % plot([0.05 0.5], [0 0], '-k', 'LineWidth', 2, 'MarkerSize',6)
% % plot([0.05 0.5], 0.6*[1 1], '--k', 'LineWidth', 2, 'MarkerSize',6)
% hold off
% 
% axis([0 40 0 0.5])
% % axis([0.05 0.8 0 0.3])
% set(gca, 'FontName', 'TimesNewRoman', 'FontSize', FontSizeAx)
% set(gca, 'Units', 'normalized', 'Position', axespos);
% xlabel('$$k$$','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 0,'Units', 'Normalize','Position', xlabelpos,'Interpreter', 'LaTeX');
% ylabel('$$\eta$$','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 90,'Units', 'Normalize','Position', ylabelpos,'Interpreter', 'LaTeX');
% % print('-depsc', '-r600','Ct_v_St_NoD.eps');
% 
% 
% %% 
% % figure
% % 
% % set(gcf, 'Units', 'centimeters');
% % set(gcf, 'Position', afFigurePosition); % [left bottom width height]
% % set(gcf, 'PaperPositionMode', 'auto')
% % set(gcf,'DefaultAxesFontSize',FontSizeAx,'DefaultAxesFontName','TimesNewRoman','DefaultAxesGridLineStyle','-.','DefaultAxesLineWidth',2,'DefaultAxesFontWeight','Normal')
% % 
% % hold on
% % for i = 1:length(A_c_vec)
% %     plot(St_vec, Pow_avg(:,i)./(rho*c*A(:,i).^3.*f(:,i).^3), 's-','color',[(1 - colorvec(i)) 0 colorvec(i)], 'LineWidth', 2, 'MarkerSize',6)
% % end
% % % plot([0.05^3 0.5^3], [0 0], '-k', 'LineWidth', 2, 'MarkerSize',6)
% % hold off
% % 
% % % axis([0 2 0 0.4])
% % set(gca, 'FontName', 'TimesNewRoman', 'FontSize', FontSizeAx)
% % set(gca, 'Units', 'normalized', 'Position', axespos);
% % xlabel('$$St$$','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 0,'Units', 'Normalize','Position', xlabelpos,'Interpreter', 'LaTeX');
% % ylabel('$$C_P$$','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 90,'Units', 'Normalize','Position', ylabelpos,'Interpreter', 'LaTeX');
% % 
% % figure
% % 
% % set(gcf, 'Units', 'centimeters');
% % set(gcf, 'Position', afFigurePosition); % [left bottom width height]
% % set(gcf, 'PaperPositionMode', 'auto')
% % set(gcf,'DefaultAxesFontSize',FontSizeAx,'DefaultAxesFontName','TimesNewRoman','DefaultAxesGridLineStyle','-.','DefaultAxesLineWidth',2,'DefaultAxesFontWeight','Normal')
% % 
% % hold on
% % for i = 1:length(A_c_vec)
% %     plot(St_vec, Pow_avg(:,i), 's-','color',[(1 - colorvec(i)) 0 colorvec(i)], 'LineWidth', 2, 'MarkerSize',6)
% % end
% % % plot([0.05^3 0.5^3], [0 0], '-k', 'LineWidth', 2, 'MarkerSize',6)
% % hold off
% % 
% % % axis([0 2 0 0.4])
% % set(gca, 'FontName', 'TimesNewRoman', 'FontSize', FontSizeAx)
% % set(gca, 'Units', 'normalized', 'Position', axespos);
% % xlabel('$$St$$','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 0,'Units', 'Normalize','Position', xlabelpos,'Interpreter', 'LaTeX');
% % ylabel('$$C_P$$','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 90,'Units', 'Normalize','Position', ylabelpos,'Interpreter', 'LaTeX');
% % 
% % 
% % figure
% % 
% % set(gcf, 'Units', 'centimeters');
% % set(gcf, 'Position', afFigurePosition); % [left bottom width height]
% % set(gcf, 'PaperPositionMode', 'auto')
% % set(gcf,'DefaultAxesFontSize',FontSizeAx,'DefaultAxesFontName','TimesNewRoman','DefaultAxesGridLineStyle','-.','DefaultAxesLineWidth',2,'DefaultAxesFontWeight','Normal')
% % 
% % hold on
% % for i = 1:length(A_c_vec)
% %     plot(St_vec, Pow_avg(:,i)./(rho*c^2*f(:,i).^3.*A(:,i).^3), 's-','color',[(1 - colorvec(i)) 0 colorvec(i)], 'LineWidth', 2, 'MarkerSize',6)
% % end
% % % plot([0.05^3 0.5^3], [0 0], '-k', 'LineWidth', 2, 'MarkerSize',6)
% % hold off
% % 
% % % axis([0 2 0 0.4])
% % set(gca, 'FontName', 'TimesNewRoman', 'FontSize', FontSizeAx)
% % set(gca, 'Units', 'normalized', 'Position', axespos);
% % xlabel('$$St$$','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 0,'Units', 'Normalize','Position', xlabelpos,'Interpreter', 'LaTeX');
% % ylabel('$$C_P$$','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 90,'Units', 'Normalize','Position', ylabelpos,'Interpreter', 'LaTeX');
% % 

