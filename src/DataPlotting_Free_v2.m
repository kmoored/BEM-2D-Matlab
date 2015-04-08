clear
close all
clc


%% Data
load('FreeV_2/FreeV_Pitch_Data_2.mat')
Beta = 1;
Cd_bod = 0.015;
c = 0.1;
b = 2*c;
rho = 1000;
A_bod = c*b*A_bp;
% Beta = [1 1 1 1]';
% Beta = asin(1/2*A_c_vec)./A_c_vec
Beta = A_c_vec.^(1/2);
% ct_prop = 2.25./[1 1 1 1];
ct_prop = 2.25*[1 1 1 1]';

Uc = mean(Q0(:,:,(Ncyc - 1)*Nstep + 1:end),3);
D_bod = 1/2*Cd_bod*rho*A_bod*Uc.^2;
np_1 = D_bod.*Ec;

Ct = Tnet_avg./(1/2*rho*c*b*Uc.^2);
Cp = Pow_avg./(1/2*rho*c*b*Uc.^3);
np_2 = Ct./Cp;

np_3 = Tnet_avg.*Uc./Pow_avg;
Tnet_avg;
D_bod;
np;


%% Plotting
FontSizeAx = 24;
FontSizeLb = 32;
afFigurePosition = [1 1 25 15];
ylabelpos = [-0.1 0.5];
xlabelpos = [0.5 -0.14];
axespos = [0.15 0.2 0.82 0.75];

num_vec = 1:4;
colorvec = -((num_vec-1)/3 - 1);


figure

set(gcf, 'Units', 'centimeters');
set(gcf, 'Position', afFigurePosition); % [left bottom width height]
set(gcf, 'PaperPositionMode', 'auto')
set(gcf,'DefaultAxesFontSize',FontSizeAx,'DefaultAxesFontName','TimesNewRoman','DefaultAxesGridLineStyle','-.','DefaultAxesLineWidth',2,'DefaultAxesFontWeight','Normal')

hold on
for i = 1:length(f_vec)
    f = f_vec(i);
    
    for j = 1:length(A_c_vec)
        A_c = A_c_vec(j);
        beta = Beta(j);
        Ct_prop = ct_prop(j);
        
        delT = 1/f/Nstep;
        t = ((1:Ncyc*Nstep+1) - 1)*delT;
        
        U_t = reshape(Q0(i,j,:),length(t),1);
        plot(t.*f,U_t, '-','color',[(1 - colorvec(j)) 0 colorvec(j)], 'LineWidth', 2, 'MarkerSize',6)
    end
end
hold off

axis([0 15 0 3.75])
set(gca, 'FontName', 'TimesNewRoman', 'FontSize', FontSizeAx)
set(gca, 'Units', 'normalized', 'Position', axespos);
xlabel('$$t/T$$','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 0,'Units', 'Normalize','Position', xlabelpos,'Interpreter', 'LaTeX');
ylabel('$$U(t)$$, m/s','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 90,'Units', 'Normalize','Position', ylabelpos,'Interpreter', 'LaTeX');
print('-depsc', '-r600','Ut_v_t_fA_vary.eps');


figure

set(gcf, 'Units', 'centimeters');
set(gcf, 'Position', afFigurePosition); % [left bottom width height]
set(gcf, 'PaperPositionMode', 'auto')
set(gcf,'DefaultAxesFontSize',FontSizeAx,'DefaultAxesFontName','TimesNewRoman','DefaultAxesGridLineStyle','-.','DefaultAxesLineWidth',2,'DefaultAxesFontWeight','Normal')

hold on
for i = 1:length(f_vec)
    f = f_vec(i);
    
    for j = 1:length(A_c_vec)
        A_c = A_c_vec(j);
        beta = Beta(j);
        Ct_prop = ct_prop(j);
%         M_star = 1./beta;
        
        delT = 1/f/Nstep;
        t = ((1:Ncyc*Nstep+1) - 1)*delT;
        
        U_t = reshape(Q0(i,j,:),length(t),1);
%         plot(t/(M_star/(Ct_prop*f)*sqrt(2*(1/A_bp)*Ct_prop/Cd_bod)),U_t, '-','color',[(1 - colorvec(j)) 0 colorvec(j)], 'LineWidth', 2, 'MarkerSize',6)

        plot(t/(M_star/(Ct_prop*f)*sqrt(2/A_bp*Ct_prop/Cd_bod)),U_t/(f*A_c*c*sqrt(2/A_bp*Ct_prop/Cd_bod)), '-','color',[(1 - colorvec(j)) 0 colorvec(j)], 'LineWidth', 2, 'MarkerSize',6)
    end
end
hold off

% axis([0 4 0 1.5])
set(gca, 'FontName', 'TimesNewRoman', 'FontSize', FontSizeAx)
set(gca, 'Units', 'normalized', 'Position', axespos);
xlabel('$$\hat{t}$$','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 0,'Units', 'Normalize','Position', xlabelpos,'Interpreter', 'LaTeX');
ylabel('$$\hat{U}$$','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 90,'Units', 'Normalize','Position', ylabelpos,'Interpreter', 'LaTeX');
print('-depsc', '-r600','Uhat_v_that_fA_vary.eps');


figure

set(gcf, 'Units', 'centimeters');
set(gcf, 'Position', afFigurePosition); % [left bottom width height]
set(gcf, 'PaperPositionMode', 'auto')
set(gcf,'DefaultAxesFontSize',FontSizeAx,'DefaultAxesFontName','TimesNewRoman','DefaultAxesGridLineStyle','-.','DefaultAxesLineWidth',2,'DefaultAxesFontWeight','Normal')

hold on
for j = 1:length(A_c_vec)
    A_c = A_c_vec(j);
    Ct_prop = ct_prop(j);
    
    plot(Uc(:,j),Ec(:,j),'s-','color',[(1 - colorvec(j)) 0 colorvec(j)], 'LineWidth', 2, 'MarkerSize',6) 
end
hold off

axis([0 3.5 0 5])
set(gca, 'FontName', 'TimesNewRoman', 'FontSize', FontSizeAx)
set(gca, 'Units', 'normalized', 'Position', axespos);
xlabel('$$U_c$$, m/s','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 0,'Units', 'Normalize','Position', xlabelpos,'Interpreter', 'LaTeX');
ylabel('$$\xi$$, km/kJ','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 90,'Units', 'Normalize','Position', ylabelpos,'Interpreter', 'LaTeX');
print('-depsc', '-r600','Ec_v_U_fA_vary.eps');

figure

set(gcf, 'Units', 'centimeters');
set(gcf, 'Position', afFigurePosition); % [left bottom width height]
set(gcf, 'PaperPositionMode', 'auto')
set(gcf,'DefaultAxesFontSize',FontSizeAx,'DefaultAxesFontName','TimesNewRoman','DefaultAxesGridLineStyle','-.','DefaultAxesLineWidth',2,'DefaultAxesFontWeight','Normal')

hold on
for j = 1:length(A_c_vec)
    A_c = A_c_vec(j);
    Ct_prop = ct_prop(j);
    beta = Beta(j);
    alpha_max = 1/2*(atan2(2*pi*f_vec*A_c*c,Uc(:,j)) - asin(1/2*A_c)*ones(20,1));
    
%         plot(Uc(:,j),Ec(:,j)./(f_vec*A_c*c*sqrt(2/A_bp*Ct_prop/Cd_bod)).*(rho*c^2*(2*c)*beta.*f_vec.^3*(A_c*c)^2),'s-','color',[(1 - colorvec(j)) 0 colorvec(j)], 'LineWidth', 2, 'MarkerSize',6) 
        plot(Uc(:,j),Ec(:,j).*(2*c)*c.*Uc(:,j).^2,'s-','color',[(1 - colorvec(j)) 0 colorvec(j)], 'LineWidth', 2, 'MarkerSize',6) 

%     plot(Uc(:,j),Ec(:,j).*sin(alpha_max)./(f_vec*A_c*c*sqrt(2/A_bp*Ct_prop/Cd_bod)).*(rho*c*(2*c).*f_vec.^3*(A_c*c)^3),'s-','color',[(1 - colorvec(j)) 0 colorvec(j)], 'LineWidth', 2, 'MarkerSize',6) 
%     plot(Uc(:,j),Ec(:,j)./(f_vec*A_c*c*sqrt(2/A_bp*Ct_prop/Cd_bod)).*(rho*c*(2*c).*f_vec.^3*(A_c*c)^3),'s-','color',[(1 - colorvec(j)) 0 colorvec(j)], 'LineWidth', 2, 'MarkerSize',6) 
    
%     plot(2*pi*f_vec*c./Uc(:,j),Ec(:,j)./(f_vec*A_c*c*sqrt(2/A_bp*Ct_prop/Cd_bod)).*(rho*c*(2*c).*f_vec.^3*(A_c*c)^3),'s-','color',[(1 - colorvec(j)) 0 colorvec(j)], 'LineWidth', 2, 'MarkerSize',6) 
end
hold off

axis([0 3.5 0 0.1])
set(gca, 'FontName', 'TimesNewRoman', 'FontSize', FontSizeAx)
set(gca, 'Units', 'normalized', 'Position', axespos);
xlabel('$$U_c$$, m/s','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 0,'Units', 'Normalize','Position', xlabelpos,'Interpreter', 'LaTeX');
ylabel('$$\hat{\xi}$$','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 90,'Units', 'Normalize','Position', ylabelpos,'Interpreter', 'LaTeX');
print('-depsc', '-r600','Echat_v_Uc_fA_vary.eps');
