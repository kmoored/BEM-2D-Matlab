clear
close all
clc


%% Data 

% Conventional UUVs -- [Length (m), Mass (kg), Velocity (m/s), Economy (km/kJ), COT (J/kg/m)]
Remus = [3.25 240 2.1 0.0095/0.5 0.5/0.0095/240];
Iver2 = [1.78 29.5 1.29 0.03/0.5 0.5/0.03/29.5];
Slocum = [1.5 54 0.25 5.34 1/5.34/54];

% Biology/BAUVs -- [Length (m), Mass (kg), Frequency (Hz), Velocity (m/s), Economy (km/kJ), COT (J/kg/m)]
Mantabot_UVa = [0.43 NaN NaN 0.43 0.06/0.1 2.15];
Mantabot_P = [0.44 11.86 0.375 0.33 2.34 1/2.34/11.86]';
Mantaray = [1.97*ones(4,1), 740.8*ones(4,1), [0.15 0.25 0.5 0.75]', [0.46 0.78 1.57 2.37]', [0.448 0.1956 0.056 0.027]', 1/740.8./[0.448 0.1956 0.056 0.027]'];
Stingray = [0.19*ones(2,1), 0.4*ones(2,1), [1.25 1.5]', [0.086 0.108]', [99.02 81.85]', 1/0.4./[99.02 81.85]'];
MantaMini = [1.97/sqrt(2) 261.9 1.25*sqrt(2) 1.59 0.1242 1/0.1242/(261.9)];



%% Plotting
FontSizeAx = 24;
FontSizeLb = 32;
afFigurePosition = [1 1 25 15];
ylabelpos = [-0.1 0.5];
xlabelpos = [0.5 -0.14];
axespos = [0.15 0.2 0.82 0.75];



%% Cost of Transport 

alpha = 0.003;
U = linspace(0.1,10,100);
Ec = alpha*U.^(2);

figure

set(gcf, 'Units', 'centimeters');
set(gcf, 'Position', afFigurePosition); % [left bottom width height]
set(gcf, 'PaperPositionMode', 'auto')
set(gcf,'DefaultAxesFontSize',FontSizeAx,'DefaultAxesFontName','TimesNewRoman','DefaultAxesGridLineStyle','-.','DefaultAxesLineWidth',2,'DefaultAxesFontWeight','Normal')

hold on
    plot(U,Ec,'--k', 'LineWidth', 1, 'MarkerSize',6)
%     plot(Slocum(3),Slocum(5),'or', 'LineWidth', 2, 'MarkerSize',6)
    
    plot(Iver2(3),Iver2(5),'sr', 'LineWidth', 2, 'MarkerSize',6)
%     errorbar(Iver2(3),Iver2(5),Iver2(5)*(0.4/0.5),Iver2(5)*(0.6/0.5),'^r', 'LineWidth', 2, 'MarkerSize',10)
    
    plot(Remus(3),Remus(5),'^r', 'LineWidth', 2, 'MarkerSize',6)
%     errorbar(Remus(3),Remus(5),Remus(5)*(0.4/0.5),Remus(5)*(0.6/0.5),'^r', 'LineWidth', 2, 'MarkerSize',10)

    
    plot(Mantabot_UVa(4),Mantabot_UVa(6),'pb', 'LineWidth', 2, 'MarkerSize',6)
    plot(Mantabot_P(4),Mantabot_P(6),'sb', 'LineWidth', 2, 'MarkerSize',6)
    plot(Mantaray(:,4),Mantaray(:,6),'pg', 'LineWidth', 2, 'MarkerSize',6)
    plot(Stingray(:,4),Stingray(:,6),'vg', 'LineWidth', 2, 'MarkerSize',6)
    plot(MantaMini(:,4),MantaMini(:,6),'og', 'LineWidth', 2, 'MarkerSize',6)

%     plot(Iver2(3),Iver2(5),'sr', 'LineWidth', 2, 'MarkerSize',6)
%     plot(Remus(3),Remus(5),'^r', 'LineWidth', 2, 'MarkerSize',6)

%     plot([Slocum_Ec(1);2.5;4],[Slocum_Ec(2);Iver2_Ec;Remus_Ec],'--k', 'LineWidth', 1, 'MarkerSize',6)
%     plot(1.94*Uc_4,Ec_4,'-b', 'LineWidth', 2, 'MarkerSize',6) 
%     plot(1.94*Uc_5,Ec_5,'-b', 'LineWidth', 2, 'MarkerSize',6) 
%     plot(1.94*Uc_6,Ec_6,'-b', 'LineWidth', 2, 'MarkerSize',6) 

%     plot(1.94*Uc_5,Ec_5,'-g', 'LineWidth', 2, 'MarkerSize',6) 
%     plot(1.94*Uc_est,Ec_est,'--b', 'LineWidth', 2, 'MarkerSize',6)
%     plot(Mantabot_P(1),Mantabot_P(2),'ob', 'LineWidth', 2, 'MarkerSize',10)
%     plot(0.84,Mantabot_UVa,'pb', 'LineWidth', 2, 'MarkerSize',10)
%     plot(MantaRay(:,1),MantaRay(:,2),'sb', 'LineWidth', 2, 'MarkerSize',10)
%     plot(MantaMini(:,1)*1.97,MantaMini(:,2),'pg', 'LineWidth', 2, 'MarkerSize',10)
%     plot(Slocum_Ec(1),Slocum_Ec(2),'or', 'LineWidth', 2, 'MarkerSize',10)

%     plot(Panel_2D(1),Panel_2D(2),'vb', 'LineWidth', 2, 'MarkerSize',10)
hold off

% axis([0 6 0 2.5])
axis([0.05 10 0.001 10])
set(gca, 'YScale', 'log', 'XScale','log')
set(gca, 'FontName', 'TimesNewRoman', 'FontSize', FontSizeAx)
set(gca, 'Units', 'normalized', 'Position', axespos);
xlabel('$$U_c$$, m/s','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 0,'Units', 'Normalize','Position', xlabelpos,'Interpreter', 'LaTeX');
ylabel('$$COT$$, J/kg/m','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 90,'Units', 'Normalize','Position', ylabelpos,'Interpreter', 'LaTeX');
print('-depsc', '-r600','COT_v_U.eps');


% %% Cost of Transport
% alpha = 1.2;
% U = linspace(0.1,10,100);
% Ec = alpha./U.^(2);
% 
% figure
% 
% set(gcf, 'Units', 'centimeters');
% set(gcf, 'Position', afFigurePosition); % [left bottom width height]
% set(gcf, 'PaperPositionMode', 'auto')
% set(gcf,'DefaultAxesFontSize',FontSizeAx,'DefaultAxesFontName','TimesNewRoman','DefaultAxesGridLineStyle','-.','DefaultAxesLineWidth',2,'DefaultAxesFontWeight','Normal')
% 
% hold on
% %     plot(U,Ec,'--k', 'LineWidth', 1, 'MarkerSize',6)
%     plot([Slocum_Ec(1);MantaRay(3:4,1)],[Slocum_Ec(2);MantaRay(3:4,2)],'--k', 'LineWidth', 1, 'MarkerSize',6)
% 
% %     plot([Slocum_Ec(1);2.5;4],[Slocum_Ec(2);Iver2_Ec;Remus_Ec],'--k', 'LineWidth', 1, 'MarkerSize',6)
% %     plot(1.94*Uc_4,Ec_4,'-b', 'LineWidth', 2, 'MarkerSize',6) 
% %     plot(1.94*Uc_5,Ec_5,'-b', 'LineWidth', 2, 'MarkerSize',6) 
% %     plot(1.94*Uc_6,Ec_6,'-b', 'LineWidth', 2, 'MarkerSize',6) 
% 
% %     plot(1.94*Uc_5,Ec_5,'-g', 'LineWidth', 2, 'MarkerSize',6) 
% %     plot(1.94*Uc_est,Ec_est,'--b', 'LineWidth', 2, 'MarkerSize',6)
%     errorbar(4,Remus_Ec,Remus_Ec*(0.4/0.5),Remus_Ec*(0.6/0.5),'sr', 'LineWidth', 2, 'MarkerSize',10)
%     errorbar(2.5,Iver2_Ec,Iver2_Ec*(0.4/0.5),Iver2_Ec*(0.6/0.5),'^r', 'LineWidth', 2, 'MarkerSize',10)
%     plot(Mantabot_P(1),Mantabot_P(2),'ob', 'LineWidth', 2, 'MarkerSize',10)
%     plot(0.84,Mantabot_UVa,'pb', 'LineWidth', 2, 'MarkerSize',10)
%     plot(MantaRay(:,1),MantaRay(:,2),'sb', 'LineWidth', 2, 'MarkerSize',10)
%     plot(MantaMini(:,1)*1.97,MantaMini(:,2),'pg', 'LineWidth', 2, 'MarkerSize',10)
%     plot(Slocum_Ec(1),Slocum_Ec(2),'or', 'LineWidth', 2, 'MarkerSize',10)
% 
% %     plot(Panel_2D(1),Panel_2D(2),'vb', 'LineWidth', 2, 'MarkerSize',10)
% hold off
% 
% % axis([0 6 0 2.5])
% axis([0.4 10 0.005 10])
% set(gca, 'YScale', 'log', 'XScale','log')
% set(gca, 'FontName', 'TimesNewRoman', 'FontSize', FontSizeAx)
% set(gca, 'Units', 'normalized', 'Position', axespos);
% xlabel('$$U_c$$, knots','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 0,'Units', 'Normalize','Position', xlabelpos,'Interpreter', 'LaTeX');
% ylabel('$$\xi$$, km/kJ','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 90,'Units', 'Normalize','Position', ylabelpos,'Interpreter', 'LaTeX');
% print('-depsc', '-r600','Ec_v_U.eps');
% 
% 
% 
% %% Normalizing speed by BL
% % figure
% % 
% % set(gcf, 'Units', 'centimeters');
% % set(gcf, 'Position', afFigurePosition); % [left bottom width height]
% % set(gcf, 'PaperPositionMode', 'auto')
% % set(gcf,'DefaultAxesFontSize',FontSizeAx,'DefaultAxesFontName','TimesNewRoman','DefaultAxesGridLineStyle','-.','DefaultAxesLineWidth',2,'DefaultAxesFontWeight','Normal')
% % 
% % 
% % hold on
% % %     plot(Uc_4,Ec_4*(0.5^2*0.001).*Uc_4.^2,'-b', 'LineWidth', 2, 'MarkerSize',6) 
% % %     plot(Uc_5,Ec_5*(0.001).*Uc_5.^2,'-r', 'LineWidth', 2, 'MarkerSize',6) 
% % %     plot(Uc_6,Ec_6*(1/3^2*0.001).*Uc_6.^2,'-g', 'LineWidth', 2, 'MarkerSize',6) 
% %     errorbar(4/1.94/3.25,Remus_Ec/(3.25*0.001),Remus_Ec/(3.25*0.001)*(0.4/0.5),Remus_Ec/(3.25*0.001)*(0.6/0.5),'sr', 'LineWidth', 2, 'MarkerSize',10)
% %     errorbar(2.5/1.94/1.78,Iver2_Ec/(1.78*0.001),Iver2_Ec/(1.78*0.001)*(0.4/0.5),Iver2_Ec/(1.78*0.001)*(0.6/0.5),'^r', 'LineWidth', 2, 'MarkerSize',10)
% % %     plot(4/1.94/3.25,Remus_Ec/(3.25*0.001),'ok', 'LineWidth', 2, 'MarkerSize',10)
% % %     plot(2.5/1.94/1.78,Iver2_Ec/(1.78*0.001),'^k', 'LineWidth', 2, 'MarkerSize',10)
% %     plot(Mantabot_P(1)/1.94/0.28,Mantabot_P(2)/(0.28*0.001),'ob', 'LineWidth', 2, 'MarkerSize',10)
% %     plot(1,Mantabot_UVa/(0.43*0.001),'pb', 'LineWidth', 2, 'MarkerSize',10)
% %     plot(MantaRay(1)/1.94/1.63,MantaRay(2)/(1.63*0.001),'sb', 'LineWidth', 2, 'MarkerSize',10)
% %     plot(Slocum_Ec(1)/1.7,Slocum_Ec(2)/(1.7*0.001),'or', 'LineWidth', 2, 'MarkerSize',10)
% % hold off
% % 
% % % axis([0 1.94*6 0 2.5])
% % % axis([0.3 5 0.01 10])
% % set(gca, 'YScale', 'log', 'XScale','linear')
% % set(gca, 'FontName', 'TimesNewRoman', 'FontSize', FontSizeAx)
% % set(gca, 'Units', 'normalized', 'Position', axespos);
% % xlabel('$$U_c$$, BL/s','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 0,'Units', 'Normalize','Position', xlabelpos,'Interpreter', 'LaTeX');
% % ylabel('$$\xi$$, BL/kJ','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 90,'Units', 'Normalize','Position', ylabelpos,'Interpreter', 'LaTeX');
% % print('-depsc', '-r600','Ec_v_U_BL.eps');
% % 
% % % Normalizing speed by BL
% % figure
% % 
% % set(gcf, 'Units', 'centimeters');
% % set(gcf, 'Position', afFigurePosition); % [left bottom width height]
% % set(gcf, 'PaperPositionMode', 'auto')
% % set(gcf,'DefaultAxesFontSize',FontSizeAx,'DefaultAxesFontName','TimesNewRoman','DefaultAxesGridLineStyle','-.','DefaultAxesLineWidth',2,'DefaultAxesFontWeight','Normal')
% % 
% % 
% % hold on
% % %     plot(Uc_4,Ec_4*(0.5^2*0.001).*Uc_4.^2,'-b', 'LineWidth', 2, 'MarkerSize',6) 
% % %     plot(Uc_5,Ec_5*(0.001).*Uc_5.^2,'-r', 'LineWidth', 2, 'MarkerSize',6) 
% % %     plot(Uc_6,Ec_6*(1/3^2*0.001).*Uc_6.^2,'-g', 'LineWidth', 2, 'MarkerSize',6) 
% %     errorbar(4/1.94/3.25,Remus_Ec,Remus_Ec*(0.4/0.5),Remus_Ec*(0.6/0.5),'sr', 'LineWidth', 2, 'MarkerSize',10)
% %     errorbar(2.5/1.94/1.78,Iver2_Ec,Iver2_Ec*(0.4/0.5),Iver2_Ec*(0.6/0.5),'^r', 'LineWidth', 2, 'MarkerSize',10)
% % %     plot(4/1.94/3.25,Remus_Ec/(3.25*0.001),'ok', 'LineWidth', 2, 'MarkerSize',10)
% % %     plot(2.5/1.94/1.78,Iver2_Ec/(1.78*0.001),'^k', 'LineWidth', 2, 'MarkerSize',10)
% %     plot(Mantabot_P(1)/1.94/0.28,Mantabot_P(2),'ob', 'LineWidth', 2, 'MarkerSize',10)
% %     plot(1,Mantabot_UVa,'pb', 'LineWidth', 2, 'MarkerSize',10)
% %     plot(MantaRay(1)/1.94/1.63,MantaRay(2),'sb', 'LineWidth', 2, 'MarkerSize',10)
% %     plot(Slocum_Ec(1)/1.7,Slocum_Ec(2),'or', 'LineWidth', 2, 'MarkerSize',10)
% % hold off
% % 
% % % axis([0 1.94*6 0 2.5])
% % axis([0.2 1.2 0.001 10])
% % set(gca, 'YScale', 'log', 'XScale','linear')
% % set(gca, 'FontName', 'TimesNewRoman', 'FontSize', FontSizeAx)
% % set(gca, 'Units', 'normalized', 'Position', axespos);
% % xlabel('$$U_c$$, BL/s','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 0,'Units', 'Normalize','Position', xlabelpos,'Interpreter', 'LaTeX');
% % ylabel('$$\xi$$, km/kJ','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 90,'Units', 'Normalize','Position', ylabelpos,'Interpreter', 'LaTeX');
% % print('-depsc', '-r600','Ec_v_U_BL2.eps');
% % 
% % nu = 1.004e-6;
% % L_manta = 1.97;
% % L_sting = 0.1922;
% % 
% % 
% % Mantaray = [0.15 0.46 0.448; 0.25 0.78 0.1956; 0.5 1.57 0.056; 0.75 2.37 0.027];
% % Stingray = [1.25 0.086 99.02; 1.5 0.108 81.85];
% % 
% % alpha = 4e12*0.15;
% % U = L_manta/nu*linspace(0.5,10,100);
% % Ec = alpha./U.^(2);
% % 
% % alpha1 = 7e7*0.175;
% % U1 = L_manta/nu*linspace(0.2,2,100);
% % Ec1 = alpha1./U1.^(5/4);
% % 
% % alpha2 = 4e6*4.75;
% % U2 = L_sting/nu*linspace(0.05,0.5,100);
% % Ec2 = alpha2./U2.^(5/4);
% % 
% % 
% % % Normalizing economy
% % figure
% % 
% % set(gcf, 'Units', 'centimeters');
% % set(gcf, 'Position', afFigurePosition); % [left bottom width height]
% % set(gcf, 'PaperPositionMode', 'auto')
% % set(gcf,'DefaultAxesFontSize',FontSizeAx,'DefaultAxesFontName','TimesNewRoman','DefaultAxesGridLineStyle','-.','DefaultAxesLineWidth',2,'DefaultAxesFontWeight','Normal')
% % 
% % 
% % hold on
% % %     plot(U,Ec,'--b', 'LineWidth', 1, 'MarkerSize',6)
% % %     plot(U1,Ec1,'--r', 'LineWidth', 1, 'MarkerSize',6)
% % %     plot(U2,Ec2,'--r', 'LineWidth', 1, 'MarkerSize',6)
% % 
% %     plot(L_manta/nu*Mantaray(:,2),Mantaray(:,3),'sk', 'LineWidth', 2, 'MarkerSize',10)
% %     plot(L_sting/nu*Stingray(:,2),Stingray(:,3),'ok', 'LineWidth', 2, 'MarkerSize',10)
% % hold off
% % 
% % % axis([0 1.94*6 0 2.5])
% % % axis([0.05 10 1e-2 1e3])
% % axis([1e4 1e7 1e-2 1e3])
% % set(gca, 'YScale', 'log', 'XScale','log')
% % set(gca, 'FontName', 'TimesNewRoman', 'FontSize', FontSizeAx)
% % set(gca, 'Units', 'normalized', 'Position', axespos);
% % xlabel('$$Re$$','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 0,'Units', 'Normalize','Position', xlabelpos,'Interpreter', 'LaTeX');
% % ylabel('$$\xi$$, km/kJ','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 90,'Units', 'Normalize','Position', ylabelpos,'Interpreter', 'LaTeX');
% % print('-depsc', '-r600','Ec_v_Re_rays.eps');
% % 
% % 
% % %%
% % nu = 1.004e-6;
% % L_manta = 1.97;
% % L_sting = 0.1922;
% % 
% % 
% % Mantaray = [0.15 0.46 0.448; 0.25 0.78 0.1956; 0.5 1.57 0.056; 0.75 2.37 0.027];
% % Stingray = [1.25 0.086 99.02; 1.5 0.108 81.85];
% % 
% % alpha = 0.15;
% % U = linspace(0.5,10,100);
% % Ec = alpha./U.^(2);
% % 
% % alpha1 = 0.175;
% % U1 = linspace(0.2,2,100);
% % Ec1 = alpha1./U1.^(5/4);
% % 
% % alpha2 = 4.75;
% % U2 = linspace(0.05,0.5,100);
% % Ec2 = alpha2./U2.^(5/4);
% % 
% % 
% % % Normalizing economy
% % figure
% % 
% % set(gcf, 'Units', 'centimeters');
% % set(gcf, 'Position', afFigurePosition); % [left bottom width height]
% % set(gcf, 'PaperPositionMode', 'auto')
% % set(gcf,'DefaultAxesFontSize',FontSizeAx,'DefaultAxesFontName','TimesNewRoman','DefaultAxesGridLineStyle','-.','DefaultAxesLineWidth',2,'DefaultAxesFontWeight','Normal')
% % 
% % 
% % hold on
% %     plot(U,Ec,'--b', 'LineWidth', 1, 'MarkerSize',6)
% %     plot(U1,Ec1,'--r', 'LineWidth', 1, 'MarkerSize',6)
% %     plot(U2,Ec2,'--r', 'LineWidth', 1, 'MarkerSize',6)
% %     plot(Mantaray(:,2),Mantaray(:,3),'sk', 'LineWidth', 2, 'MarkerSize',10)
% %     plot(Stingray(:,2),Stingray(:,3),'ok', 'LineWidth', 2, 'MarkerSize',10)
% % 
% % %     plot(L_manta/nu*Mantaray(:,2),Mantaray(:,3),'sk', 'LineWidth', 2, 'MarkerSize',10)
% % %     plot(L_sting/nu*Stingray(:,2),Stingray(:,3),'ok', 'LineWidth', 2, 'MarkerSize',10)
% % hold off
% % 
% % % axis([0 1.94*6 0 2.5])
% % axis([0.05 10 1e-2 1e3])
% % set(gca, 'YScale', 'log', 'XScale','log')
% % set(gca, 'FontName', 'TimesNewRoman', 'FontSize', FontSizeAx)
% % set(gca, 'Units', 'normalized', 'Position', axespos);
% % xlabel('$$U$$, m/s','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 0,'Units', 'Normalize','Position', xlabelpos,'Interpreter', 'LaTeX');
% % ylabel('$$\xi$$, km/kJ','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 90,'Units', 'Normalize','Position', ylabelpos,'Interpreter', 'LaTeX');
% % print('-depsc', '-r600','Ec_v_U_rays3.eps');
% % 
% % 
% % 
% % %%
% % nu = 1.004e-6;
% % L_manta = 1.97;
% % L_sting = 0.1922;
% % 
% % 
% % Mantaray = [0.15 0.46 0.448; 0.25 0.78 0.1956; 0.5 1.57 0.056; 0.75 2.37 0.027];
% % Stingray = [1.25 0.086 99.02; 1.5 0.108 81.85];
% % 
% % alpha = 1/15.75;
% % f = linspace(0,2,100);
% % U = alpha*f.^(4/3);
% % 
% % alpha = 1/14.5;
% % f1 = linspace(0,2,100);
% % U1 = alpha*f1;
% % 
% % 
% % % Normalizing economy
% % figure
% % 
% % set(gcf, 'Units', 'centimeters');
% % set(gcf, 'Position', afFigurePosition); % [left bottom width height]
% % set(gcf, 'PaperPositionMode', 'auto')
% % set(gcf,'DefaultAxesFontSize',FontSizeAx,'DefaultAxesFontName','TimesNewRoman','DefaultAxesGridLineStyle','-.','DefaultAxesLineWidth',2,'DefaultAxesFontWeight','Normal')
% % 
% % 
% % hold on
% % %     plot(f,U,'--r', 'LineWidth', 1, 'MarkerSize',6)
% % %     plot(f1,U1,'--b', 'LineWidth', 1, 'MarkerSize',6)
% % %     plot(U2,Ec2,'--r', 'LineWidth', 1, 'MarkerSize',6)
% % %     plot(Mantaray(:,1),Mantaray(:,2),'sk', 'LineWidth', 2, 'MarkerSize',10)
% %     plot(Stingray(:,1),Stingray(:,2),'ok', 'LineWidth', 2, 'MarkerSize',10)
% % 
% % %     plot(L_manta/nu*Mantaray(:,2),Mantaray(:,3),'sk', 'LineWidth', 2, 'MarkerSize',10)
% % %     plot(L_sting/nu*Stingray(:,2),Stingray(:,3),'ok', 'LineWidth', 2, 'MarkerSize',10)
% % hold off
% % 
% % % axis([0 1.94*6 0 2.5])
% % axis([0 2 0 0.2])
% % % set(gca, 'YScale', 'log', 'XScale','log')
% % set(gca, 'FontName', 'TimesNewRoman', 'FontSize', FontSizeAx)
% % set(gca, 'Units', 'normalized', 'Position', axespos);
% % xlabel('$$f$$, Hz','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 0,'Units', 'Normalize','Position', xlabelpos,'Interpreter', 'LaTeX');
% % ylabel('$$U$$, m/s','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 90,'Units', 'Normalize','Position', ylabelpos,'Interpreter', 'LaTeX');
% % print('-depsc', '-r600','U_v_f_Stingray3.eps');
% % 
% % 
% % 
% % alpha = 1/15.75;
% % f = linspace(0,2,100);
% % U = alpha*f.^(4/3);
% % 
% % alpha = 3.15;
% % f1 = linspace(0,1,100);
% % U1 = alpha*f1;
% % 
% % 
% % % Normalizing economy
% % figure
% % 
% % set(gcf, 'Units', 'centimeters');
% % set(gcf, 'Position', afFigurePosition); % [left bottom width height]
% % set(gcf, 'PaperPositionMode', 'auto')
% % set(gcf,'DefaultAxesFontSize',FontSizeAx,'DefaultAxesFontName','TimesNewRoman','DefaultAxesGridLineStyle','-.','DefaultAxesLineWidth',2,'DefaultAxesFontWeight','Normal')
% % 
% % 
% % hold on
% % %     plot(f,U,'--r', 'LineWidth', 1, 'MarkerSize',6)
% % %     plot(f1,U1,'--b', 'LineWidth', 1, 'MarkerSize',6)
% % %     plot(U2,Ec2,'--r', 'LineWidth', 1, 'MarkerSize',6)
% %     plot(Mantaray(:,1),Mantaray(:,2),'sk', 'LineWidth', 2, 'MarkerSize',10)
% % %     plot(Stingray(:,1),Stingray(:,2),'ok', 'LineWidth', 2, 'MarkerSize',10)
% % 
% % %     plot(L_manta/nu*Mantaray(:,2),Mantaray(:,3),'sk', 'LineWidth', 2, 'MarkerSize',10)
% % %     plot(L_sting/nu*Stingray(:,2),Stingray(:,3),'ok', 'LineWidth', 2, 'MarkerSize',10)
% % hold off
% % 
% % % axis([0 1.94*6 0 2.5])
% % axis([0 0.8 0 2.5])
% % % set(gca, 'YScale', 'log', 'XScale','log')
% % set(gca, 'FontName', 'TimesNewRoman', 'FontSize', FontSizeAx)
% % set(gca, 'Units', 'normalized', 'Position', axespos);
% % xlabel('$$f$$, Hz','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 0,'Units', 'Normalize','Position', xlabelpos,'Interpreter', 'LaTeX');
% % ylabel('$$U$$, m/s','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 90,'Units', 'Normalize','Position', ylabelpos,'Interpreter', 'LaTeX');
% % print('-depsc', '-r600','U_v_f_Mantaray2.eps');