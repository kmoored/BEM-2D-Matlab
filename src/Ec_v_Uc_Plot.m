clear
close all
clc

f = [0.1:0.01:0.75]';
c = 1/2;
Sb = 1.8;
Sp = 1/2*Sb;
Ct = 2.25;
Cd = 0.01;
rho = 1000;
A_c = 0.5;
A = A_c*c;
Beta_2 = 0.12;

load('FreeV_4/FreeV_Pitch_Data_4.mat')

Uc_4 = Ucruise;
Ec_4 = Ec;

load('FreeV_5/FreeV_Pitch_Data_5.mat')

Uc_5 = Ucruise;
Ec_5 = Ec;

load('FreeV_6/FreeV_Pitch_Data_6.mat')

Uc_6 = Ucruise;
Ec_6 = Ec;

% Uc = f*A*sqrt(2*Sb/Sp*Ct/Cd);
% Ec = Beta_2*sqrt(2*Sb/Sp*Ct/Cd)./(rho*Sp*c*f.^2*A*sqrt(A/c));

Ct = Ct;
Uc_est = f*A*sqrt(2*Sb/Sp*Ct/Cd);
Ec_est = Beta_2*sqrt(2*Sb/Sp*Ct/Cd)./(rho*Sp*c*f.^2*A*sqrt(A/c));

% Eci = interp1(Uc,Ec,2.57)

Remus_Ec = 0.0095./[0.5]';
Iver2_Ec = 0.03./[0.5]';
Mantabot_P = [0.33*1.94 2.34]';
Panel_2D = [3.62 0.2];
Mantabot_UVa = 0.06./[0.1]';
MantaRay = [1.94*0.46 0.448; 1.94*0.78 0.1956;1.94*1.57 0.056; 1.94*2.37 0.027];
MantaMini = [1.59 0.1242];
Slocum_Ec = [1.94*0.25 5.34];


%% Plotting
FontSizeAx = 24;
FontSizeLb = 32;
afFigurePosition = [1 1 25 15];
ylabelpos = [-0.1 0.5];
xlabelpos = [0.5 -0.14];
axespos = [0.15 0.2 0.82 0.75];



% figure
% 
% set(gcf, 'Units', 'centimeters');
% set(gcf, 'Position', afFigurePosition); % [left bottom width height]
% set(gcf, 'PaperPositionMode', 'auto')
% set(gcf,'DefaultAxesFontSize',FontSizeAx,'DefaultAxesFontName','TimesNewRoman','DefaultAxesGridLineStyle','-.','DefaultAxesLineWidth',2,'DefaultAxesFontWeight','Normal')
% 
% hold on
% for i = 1:length(f_vec)
%     f = f_vec(i);    
%     delT = 1/f/Nstep;
%     t = ((1:Ncyc*Nstep+1) - 1)*delT;
% 
%     U_t = reshape(Q0(i,1,:),length(t),1);
%     plot(t.*f,U_t, '-k', 'LineWidth', 2, 'MarkerSize',6)
% end
% hold off
% 
% % axis([0 15 0 3.75])
% set(gca, 'FontName', 'TimesNewRoman', 'FontSize', FontSizeAx)
% set(gca, 'Units', 'normalized', 'Position', axespos);
% xlabel('$$t/T$$','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 0,'Units', 'Normalize','Position', xlabelpos,'Interpreter', 'LaTeX');
% ylabel('$$U(t)$$, m/s','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 90,'Units', 'Normalize','Position', ylabelpos,'Interpreter', 'LaTeX');
% 

% figure
% 
% set(gcf, 'Units', 'centimeters');
% set(gcf, 'Position', afFigurePosition); % [left bottom width height]
% set(gcf, 'PaperPositionMode', 'auto')
% set(gcf,'DefaultAxesFontSize',FontSizeAx,'DefaultAxesFontName','TimesNewRoman','DefaultAxesGridLineStyle','-.','DefaultAxesLineWidth',2,'DefaultAxesFontWeight','Normal')
% 
% loglog([1.94*Uc_4; 4;2.5;Mantabot_P(1); 0.84; MantaRay(1)],[Ec_4;Remus_Ec;Iver2_Ec; Mantabot_P(2); Mantabot_UVa; MantaRay(2)],'ob', 'LineWidth', 2, 'MarkerSize',6) 
% 
% axis([0.5 10 0.01 10])
% set(gca, 'FontName', 'TimesNewRoman', 'FontSize', FontSizeAx)
% set(gca, 'Units', 'normalized', 'Position', axespos);
% xlabel('$$U_c$$, knots','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 0,'Units', 'Normalize','Position', xlabelpos,'Interpreter', 'LaTeX');
% ylabel('$$\xi$$, km/kJ','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 90,'Units', 'Normalize','Position', ylabelpos,'Interpreter', 'LaTeX');
% print('-depsc', '-r600','Ec_v_U_semilog.eps');
% 

alpha = 1.2;
U = linspace(0.1,10,100);
Ec = alpha./U.^(2);

figure

set(gcf, 'Units', 'centimeters');
set(gcf, 'Position', afFigurePosition); % [left bottom width height]
set(gcf, 'PaperPositionMode', 'auto')
set(gcf,'DefaultAxesFontSize',FontSizeAx,'DefaultAxesFontName','TimesNewRoman','DefaultAxesGridLineStyle','-.','DefaultAxesLineWidth',2,'DefaultAxesFontWeight','Normal')

hold on
%     plot(U,Ec,'--k', 'LineWidth', 1, 'MarkerSize',6)
    plot([Slocum_Ec(1);MantaRay(3:4,1)],[Slocum_Ec(2);MantaRay(3:4,2)],'--k', 'LineWidth', 1, 'MarkerSize',6)

%     plot([Slocum_Ec(1);2.5;4],[Slocum_Ec(2);Iver2_Ec;Remus_Ec],'--k', 'LineWidth', 1, 'MarkerSize',6)
%     plot(1.94*Uc_4,Ec_4,'-b', 'LineWidth', 2, 'MarkerSize',6) 
%     plot(1.94*Uc_5,Ec_5,'-b', 'LineWidth', 2, 'MarkerSize',6) 
%     plot(1.94*Uc_6,Ec_6,'-b', 'LineWidth', 2, 'MarkerSize',6) 

%     plot(1.94*Uc_5,Ec_5,'-g', 'LineWidth', 2, 'MarkerSize',6) 
%     plot(1.94*Uc_est,Ec_est,'--b', 'LineWidth', 2, 'MarkerSize',6)
    errorbar(4,Remus_Ec,Remus_Ec*(0.4/0.5),Remus_Ec*(0.6/0.5),'sr', 'LineWidth', 2, 'MarkerSize',10)
    errorbar(2.5,Iver2_Ec,Iver2_Ec*(0.4/0.5),Iver2_Ec*(0.6/0.5),'^r', 'LineWidth', 2, 'MarkerSize',10)
    plot(Mantabot_P(1),Mantabot_P(2),'ob', 'LineWidth', 2, 'MarkerSize',10)
    plot(0.84,Mantabot_UVa,'pb', 'LineWidth', 2, 'MarkerSize',10)
    plot(MantaRay(:,1),MantaRay(:,2),'sb', 'LineWidth', 2, 'MarkerSize',10)
    plot(MantaMini(:,1)*1.97,MantaMini(:,2),'pg', 'LineWidth', 2, 'MarkerSize',10)
    plot(Slocum_Ec(1),Slocum_Ec(2),'or', 'LineWidth', 2, 'MarkerSize',10)

%     plot(Panel_2D(1),Panel_2D(2),'vb', 'LineWidth', 2, 'MarkerSize',10)
hold off

% axis([0 6 0 2.5])
axis([0.4 10 0.005 10])
set(gca, 'YScale', 'log', 'XScale','log')
set(gca, 'FontName', 'TimesNewRoman', 'FontSize', FontSizeAx)
set(gca, 'Units', 'normalized', 'Position', axespos);
xlabel('$$U_c$$, knots','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 0,'Units', 'Normalize','Position', xlabelpos,'Interpreter', 'LaTeX');
ylabel('$$\xi$$, km/kJ','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 90,'Units', 'Normalize','Position', ylabelpos,'Interpreter', 'LaTeX');
print('-depsc', '-r600','Ec_v_U.eps');


%% Cost of Transport
alpha = 1.2;
U = linspace(0.1,10,100);
Ec = alpha./U.^(2);

figure

set(gcf, 'Units', 'centimeters');
set(gcf, 'Position', afFigurePosition); % [left bottom width height]
set(gcf, 'PaperPositionMode', 'auto')
set(gcf,'DefaultAxesFontSize',FontSizeAx,'DefaultAxesFontName','TimesNewRoman','DefaultAxesGridLineStyle','-.','DefaultAxesLineWidth',2,'DefaultAxesFontWeight','Normal')

hold on
%     plot(U,Ec,'--k', 'LineWidth', 1, 'MarkerSize',6)
    plot([Slocum_Ec(1);MantaRay(3:4,1)],[Slocum_Ec(2);MantaRay(3:4,2)],'--k', 'LineWidth', 1, 'MarkerSize',6)

%     plot([Slocum_Ec(1);2.5;4],[Slocum_Ec(2);Iver2_Ec;Remus_Ec],'--k', 'LineWidth', 1, 'MarkerSize',6)
%     plot(1.94*Uc_4,Ec_4,'-b', 'LineWidth', 2, 'MarkerSize',6) 
%     plot(1.94*Uc_5,Ec_5,'-b', 'LineWidth', 2, 'MarkerSize',6) 
%     plot(1.94*Uc_6,Ec_6,'-b', 'LineWidth', 2, 'MarkerSize',6) 

%     plot(1.94*Uc_5,Ec_5,'-g', 'LineWidth', 2, 'MarkerSize',6) 
%     plot(1.94*Uc_est,Ec_est,'--b', 'LineWidth', 2, 'MarkerSize',6)
    errorbar(4,Remus_Ec,Remus_Ec*(0.4/0.5),Remus_Ec*(0.6/0.5),'sr', 'LineWidth', 2, 'MarkerSize',10)
    errorbar(2.5,Iver2_Ec,Iver2_Ec*(0.4/0.5),Iver2_Ec*(0.6/0.5),'^r', 'LineWidth', 2, 'MarkerSize',10)
    plot(Mantabot_P(1),Mantabot_P(2),'ob', 'LineWidth', 2, 'MarkerSize',10)
    plot(0.84,Mantabot_UVa,'pb', 'LineWidth', 2, 'MarkerSize',10)
    plot(MantaRay(:,1),MantaRay(:,2),'sb', 'LineWidth', 2, 'MarkerSize',10)
    plot(MantaMini(:,1)*1.97,MantaMini(:,2),'pg', 'LineWidth', 2, 'MarkerSize',10)
    plot(Slocum_Ec(1),Slocum_Ec(2),'or', 'LineWidth', 2, 'MarkerSize',10)

%     plot(Panel_2D(1),Panel_2D(2),'vb', 'LineWidth', 2, 'MarkerSize',10)
hold off

% axis([0 6 0 2.5])
axis([0.4 10 0.005 10])
set(gca, 'YScale', 'log', 'XScale','log')
set(gca, 'FontName', 'TimesNewRoman', 'FontSize', FontSizeAx)
set(gca, 'Units', 'normalized', 'Position', axespos);
xlabel('$$U_c$$, knots','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 0,'Units', 'Normalize','Position', xlabelpos,'Interpreter', 'LaTeX');
ylabel('$$\xi$$, km/kJ','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 90,'Units', 'Normalize','Position', ylabelpos,'Interpreter', 'LaTeX');
print('-depsc', '-r600','Ec_v_U.eps');



%% Normalizing speed by BL
% figure
% 
% set(gcf, 'Units', 'centimeters');
% set(gcf, 'Position', afFigurePosition); % [left bottom width height]
% set(gcf, 'PaperPositionMode', 'auto')
% set(gcf,'DefaultAxesFontSize',FontSizeAx,'DefaultAxesFontName','TimesNewRoman','DefaultAxesGridLineStyle','-.','DefaultAxesLineWidth',2,'DefaultAxesFontWeight','Normal')
% 
% 
% hold on
% %     plot(Uc_4,Ec_4*(0.5^2*0.001).*Uc_4.^2,'-b', 'LineWidth', 2, 'MarkerSize',6) 
% %     plot(Uc_5,Ec_5*(0.001).*Uc_5.^2,'-r', 'LineWidth', 2, 'MarkerSize',6) 
% %     plot(Uc_6,Ec_6*(1/3^2*0.001).*Uc_6.^2,'-g', 'LineWidth', 2, 'MarkerSize',6) 
%     errorbar(4/1.94/3.25,Remus_Ec/(3.25*0.001),Remus_Ec/(3.25*0.001)*(0.4/0.5),Remus_Ec/(3.25*0.001)*(0.6/0.5),'sr', 'LineWidth', 2, 'MarkerSize',10)
%     errorbar(2.5/1.94/1.78,Iver2_Ec/(1.78*0.001),Iver2_Ec/(1.78*0.001)*(0.4/0.5),Iver2_Ec/(1.78*0.001)*(0.6/0.5),'^r', 'LineWidth', 2, 'MarkerSize',10)
% %     plot(4/1.94/3.25,Remus_Ec/(3.25*0.001),'ok', 'LineWidth', 2, 'MarkerSize',10)
% %     plot(2.5/1.94/1.78,Iver2_Ec/(1.78*0.001),'^k', 'LineWidth', 2, 'MarkerSize',10)
%     plot(Mantabot_P(1)/1.94/0.28,Mantabot_P(2)/(0.28*0.001),'ob', 'LineWidth', 2, 'MarkerSize',10)
%     plot(1,Mantabot_UVa/(0.43*0.001),'pb', 'LineWidth', 2, 'MarkerSize',10)
%     plot(MantaRay(1)/1.94/1.63,MantaRay(2)/(1.63*0.001),'sb', 'LineWidth', 2, 'MarkerSize',10)
%     plot(Slocum_Ec(1)/1.7,Slocum_Ec(2)/(1.7*0.001),'or', 'LineWidth', 2, 'MarkerSize',10)
% hold off
% 
% % axis([0 1.94*6 0 2.5])
% % axis([0.3 5 0.01 10])
% set(gca, 'YScale', 'log', 'XScale','linear')
% set(gca, 'FontName', 'TimesNewRoman', 'FontSize', FontSizeAx)
% set(gca, 'Units', 'normalized', 'Position', axespos);
% xlabel('$$U_c$$, BL/s','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 0,'Units', 'Normalize','Position', xlabelpos,'Interpreter', 'LaTeX');
% ylabel('$$\xi$$, BL/kJ','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 90,'Units', 'Normalize','Position', ylabelpos,'Interpreter', 'LaTeX');
% print('-depsc', '-r600','Ec_v_U_BL.eps');
% 
% % Normalizing speed by BL
% figure
% 
% set(gcf, 'Units', 'centimeters');
% set(gcf, 'Position', afFigurePosition); % [left bottom width height]
% set(gcf, 'PaperPositionMode', 'auto')
% set(gcf,'DefaultAxesFontSize',FontSizeAx,'DefaultAxesFontName','TimesNewRoman','DefaultAxesGridLineStyle','-.','DefaultAxesLineWidth',2,'DefaultAxesFontWeight','Normal')
% 
% 
% hold on
% %     plot(Uc_4,Ec_4*(0.5^2*0.001).*Uc_4.^2,'-b', 'LineWidth', 2, 'MarkerSize',6) 
% %     plot(Uc_5,Ec_5*(0.001).*Uc_5.^2,'-r', 'LineWidth', 2, 'MarkerSize',6) 
% %     plot(Uc_6,Ec_6*(1/3^2*0.001).*Uc_6.^2,'-g', 'LineWidth', 2, 'MarkerSize',6) 
%     errorbar(4/1.94/3.25,Remus_Ec,Remus_Ec*(0.4/0.5),Remus_Ec*(0.6/0.5),'sr', 'LineWidth', 2, 'MarkerSize',10)
%     errorbar(2.5/1.94/1.78,Iver2_Ec,Iver2_Ec*(0.4/0.5),Iver2_Ec*(0.6/0.5),'^r', 'LineWidth', 2, 'MarkerSize',10)
% %     plot(4/1.94/3.25,Remus_Ec/(3.25*0.001),'ok', 'LineWidth', 2, 'MarkerSize',10)
% %     plot(2.5/1.94/1.78,Iver2_Ec/(1.78*0.001),'^k', 'LineWidth', 2, 'MarkerSize',10)
%     plot(Mantabot_P(1)/1.94/0.28,Mantabot_P(2),'ob', 'LineWidth', 2, 'MarkerSize',10)
%     plot(1,Mantabot_UVa,'pb', 'LineWidth', 2, 'MarkerSize',10)
%     plot(MantaRay(1)/1.94/1.63,MantaRay(2),'sb', 'LineWidth', 2, 'MarkerSize',10)
%     plot(Slocum_Ec(1)/1.7,Slocum_Ec(2),'or', 'LineWidth', 2, 'MarkerSize',10)
% hold off
% 
% % axis([0 1.94*6 0 2.5])
% axis([0.2 1.2 0.001 10])
% set(gca, 'YScale', 'log', 'XScale','linear')
% set(gca, 'FontName', 'TimesNewRoman', 'FontSize', FontSizeAx)
% set(gca, 'Units', 'normalized', 'Position', axespos);
% xlabel('$$U_c$$, BL/s','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 0,'Units', 'Normalize','Position', xlabelpos,'Interpreter', 'LaTeX');
% ylabel('$$\xi$$, km/kJ','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 90,'Units', 'Normalize','Position', ylabelpos,'Interpreter', 'LaTeX');
% print('-depsc', '-r600','Ec_v_U_BL2.eps');
% 
% nu = 1.004e-6;
% L_manta = 1.97;
% L_sting = 0.1922;
% 
% 
% Mantaray = [0.15 0.46 0.448; 0.25 0.78 0.1956; 0.5 1.57 0.056; 0.75 2.37 0.027];
% Stingray = [1.25 0.086 99.02; 1.5 0.108 81.85];
% 
% alpha = 4e12*0.15;
% U = L_manta/nu*linspace(0.5,10,100);
% Ec = alpha./U.^(2);
% 
% alpha1 = 7e7*0.175;
% U1 = L_manta/nu*linspace(0.2,2,100);
% Ec1 = alpha1./U1.^(5/4);
% 
% alpha2 = 4e6*4.75;
% U2 = L_sting/nu*linspace(0.05,0.5,100);
% Ec2 = alpha2./U2.^(5/4);
% 
% 
% % Normalizing economy
% figure
% 
% set(gcf, 'Units', 'centimeters');
% set(gcf, 'Position', afFigurePosition); % [left bottom width height]
% set(gcf, 'PaperPositionMode', 'auto')
% set(gcf,'DefaultAxesFontSize',FontSizeAx,'DefaultAxesFontName','TimesNewRoman','DefaultAxesGridLineStyle','-.','DefaultAxesLineWidth',2,'DefaultAxesFontWeight','Normal')
% 
% 
% hold on
% %     plot(U,Ec,'--b', 'LineWidth', 1, 'MarkerSize',6)
% %     plot(U1,Ec1,'--r', 'LineWidth', 1, 'MarkerSize',6)
% %     plot(U2,Ec2,'--r', 'LineWidth', 1, 'MarkerSize',6)
% 
%     plot(L_manta/nu*Mantaray(:,2),Mantaray(:,3),'sk', 'LineWidth', 2, 'MarkerSize',10)
%     plot(L_sting/nu*Stingray(:,2),Stingray(:,3),'ok', 'LineWidth', 2, 'MarkerSize',10)
% hold off
% 
% % axis([0 1.94*6 0 2.5])
% % axis([0.05 10 1e-2 1e3])
% axis([1e4 1e7 1e-2 1e3])
% set(gca, 'YScale', 'log', 'XScale','log')
% set(gca, 'FontName', 'TimesNewRoman', 'FontSize', FontSizeAx)
% set(gca, 'Units', 'normalized', 'Position', axespos);
% xlabel('$$Re$$','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 0,'Units', 'Normalize','Position', xlabelpos,'Interpreter', 'LaTeX');
% ylabel('$$\xi$$, km/kJ','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 90,'Units', 'Normalize','Position', ylabelpos,'Interpreter', 'LaTeX');
% print('-depsc', '-r600','Ec_v_Re_rays.eps');
% 
% 
% %%
% nu = 1.004e-6;
% L_manta = 1.97;
% L_sting = 0.1922;
% 
% 
% Mantaray = [0.15 0.46 0.448; 0.25 0.78 0.1956; 0.5 1.57 0.056; 0.75 2.37 0.027];
% Stingray = [1.25 0.086 99.02; 1.5 0.108 81.85];
% 
% alpha = 0.15;
% U = linspace(0.5,10,100);
% Ec = alpha./U.^(2);
% 
% alpha1 = 0.175;
% U1 = linspace(0.2,2,100);
% Ec1 = alpha1./U1.^(5/4);
% 
% alpha2 = 4.75;
% U2 = linspace(0.05,0.5,100);
% Ec2 = alpha2./U2.^(5/4);
% 
% 
% % Normalizing economy
% figure
% 
% set(gcf, 'Units', 'centimeters');
% set(gcf, 'Position', afFigurePosition); % [left bottom width height]
% set(gcf, 'PaperPositionMode', 'auto')
% set(gcf,'DefaultAxesFontSize',FontSizeAx,'DefaultAxesFontName','TimesNewRoman','DefaultAxesGridLineStyle','-.','DefaultAxesLineWidth',2,'DefaultAxesFontWeight','Normal')
% 
% 
% hold on
%     plot(U,Ec,'--b', 'LineWidth', 1, 'MarkerSize',6)
%     plot(U1,Ec1,'--r', 'LineWidth', 1, 'MarkerSize',6)
%     plot(U2,Ec2,'--r', 'LineWidth', 1, 'MarkerSize',6)
%     plot(Mantaray(:,2),Mantaray(:,3),'sk', 'LineWidth', 2, 'MarkerSize',10)
%     plot(Stingray(:,2),Stingray(:,3),'ok', 'LineWidth', 2, 'MarkerSize',10)
% 
% %     plot(L_manta/nu*Mantaray(:,2),Mantaray(:,3),'sk', 'LineWidth', 2, 'MarkerSize',10)
% %     plot(L_sting/nu*Stingray(:,2),Stingray(:,3),'ok', 'LineWidth', 2, 'MarkerSize',10)
% hold off
% 
% % axis([0 1.94*6 0 2.5])
% axis([0.05 10 1e-2 1e3])
% set(gca, 'YScale', 'log', 'XScale','log')
% set(gca, 'FontName', 'TimesNewRoman', 'FontSize', FontSizeAx)
% set(gca, 'Units', 'normalized', 'Position', axespos);
% xlabel('$$U$$, m/s','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 0,'Units', 'Normalize','Position', xlabelpos,'Interpreter', 'LaTeX');
% ylabel('$$\xi$$, km/kJ','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 90,'Units', 'Normalize','Position', ylabelpos,'Interpreter', 'LaTeX');
% print('-depsc', '-r600','Ec_v_U_rays3.eps');
% 
% 
% 
% %%
% nu = 1.004e-6;
% L_manta = 1.97;
% L_sting = 0.1922;
% 
% 
% Mantaray = [0.15 0.46 0.448; 0.25 0.78 0.1956; 0.5 1.57 0.056; 0.75 2.37 0.027];
% Stingray = [1.25 0.086 99.02; 1.5 0.108 81.85];
% 
% alpha = 1/15.75;
% f = linspace(0,2,100);
% U = alpha*f.^(4/3);
% 
% alpha = 1/14.5;
% f1 = linspace(0,2,100);
% U1 = alpha*f1;
% 
% 
% % Normalizing economy
% figure
% 
% set(gcf, 'Units', 'centimeters');
% set(gcf, 'Position', afFigurePosition); % [left bottom width height]
% set(gcf, 'PaperPositionMode', 'auto')
% set(gcf,'DefaultAxesFontSize',FontSizeAx,'DefaultAxesFontName','TimesNewRoman','DefaultAxesGridLineStyle','-.','DefaultAxesLineWidth',2,'DefaultAxesFontWeight','Normal')
% 
% 
% hold on
% %     plot(f,U,'--r', 'LineWidth', 1, 'MarkerSize',6)
% %     plot(f1,U1,'--b', 'LineWidth', 1, 'MarkerSize',6)
% %     plot(U2,Ec2,'--r', 'LineWidth', 1, 'MarkerSize',6)
% %     plot(Mantaray(:,1),Mantaray(:,2),'sk', 'LineWidth', 2, 'MarkerSize',10)
%     plot(Stingray(:,1),Stingray(:,2),'ok', 'LineWidth', 2, 'MarkerSize',10)
% 
% %     plot(L_manta/nu*Mantaray(:,2),Mantaray(:,3),'sk', 'LineWidth', 2, 'MarkerSize',10)
% %     plot(L_sting/nu*Stingray(:,2),Stingray(:,3),'ok', 'LineWidth', 2, 'MarkerSize',10)
% hold off
% 
% % axis([0 1.94*6 0 2.5])
% axis([0 2 0 0.2])
% % set(gca, 'YScale', 'log', 'XScale','log')
% set(gca, 'FontName', 'TimesNewRoman', 'FontSize', FontSizeAx)
% set(gca, 'Units', 'normalized', 'Position', axespos);
% xlabel('$$f$$, Hz','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 0,'Units', 'Normalize','Position', xlabelpos,'Interpreter', 'LaTeX');
% ylabel('$$U$$, m/s','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 90,'Units', 'Normalize','Position', ylabelpos,'Interpreter', 'LaTeX');
% print('-depsc', '-r600','U_v_f_Stingray3.eps');
% 
% 
% 
% alpha = 1/15.75;
% f = linspace(0,2,100);
% U = alpha*f.^(4/3);
% 
% alpha = 3.15;
% f1 = linspace(0,1,100);
% U1 = alpha*f1;
% 
% 
% % Normalizing economy
% figure
% 
% set(gcf, 'Units', 'centimeters');
% set(gcf, 'Position', afFigurePosition); % [left bottom width height]
% set(gcf, 'PaperPositionMode', 'auto')
% set(gcf,'DefaultAxesFontSize',FontSizeAx,'DefaultAxesFontName','TimesNewRoman','DefaultAxesGridLineStyle','-.','DefaultAxesLineWidth',2,'DefaultAxesFontWeight','Normal')
% 
% 
% hold on
% %     plot(f,U,'--r', 'LineWidth', 1, 'MarkerSize',6)
% %     plot(f1,U1,'--b', 'LineWidth', 1, 'MarkerSize',6)
% %     plot(U2,Ec2,'--r', 'LineWidth', 1, 'MarkerSize',6)
%     plot(Mantaray(:,1),Mantaray(:,2),'sk', 'LineWidth', 2, 'MarkerSize',10)
% %     plot(Stingray(:,1),Stingray(:,2),'ok', 'LineWidth', 2, 'MarkerSize',10)
% 
% %     plot(L_manta/nu*Mantaray(:,2),Mantaray(:,3),'sk', 'LineWidth', 2, 'MarkerSize',10)
% %     plot(L_sting/nu*Stingray(:,2),Stingray(:,3),'ok', 'LineWidth', 2, 'MarkerSize',10)
% hold off
% 
% % axis([0 1.94*6 0 2.5])
% axis([0 0.8 0 2.5])
% % set(gca, 'YScale', 'log', 'XScale','log')
% set(gca, 'FontName', 'TimesNewRoman', 'FontSize', FontSizeAx)
% set(gca, 'Units', 'normalized', 'Position', axespos);
% xlabel('$$f$$, Hz','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 0,'Units', 'Normalize','Position', xlabelpos,'Interpreter', 'LaTeX');
% ylabel('$$U$$, m/s','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 90,'Units', 'Normalize','Position', ylabelpos,'Interpreter', 'LaTeX');
% print('-depsc', '-r600','U_v_f_Mantaray2.eps');