clear
close all
clc

N = 1000;
val = linspace(0.0001,40,N)';
k = 1./val;

[F,G] = LiftDefFac(k);

% Reproducing lift deficiency factor plot from Theodorsen 1935 (figure 4, pg 418).
FontSizeAx = 24;
FontSizeLb = 32;
afFigurePosition = [0 2 30 20];
axespos = [0.175 0.15 0.775 0.8];
ylabelpos = [-0.125 0.5];
xlabelpos = [0.5 -0.1];

figure
set(gcf, 'Units', 'centimeters','PaperPositionMode', 'auto','Position', afFigurePosition);
set(gcf,'DefaultAxesFontSize',FontSizeAx,'DefaultAxesFontName','TimesNewRoman','DefaultAxesGridLineStyle','-.','DefaultAxesLineWidth',1,'DefaultAxesFontWeight','Normal')


hold on
    plot(1./k,F,'-b','linewidth',2)
    plot(1./k,-G,'-k','linewidth',2)
hold off

axis([0 40 0 1])
grid on
set(gca, 'FontName', 'TimesNewRoman', 'FontSize', FontSizeAx)
set(gca, 'Units', 'normalized', 'Position', axespos);
legend('F','-G','Location','East')
xlabel('$$1/k$$','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 0,'Units', 'Normalize','Position', xlabelpos,'Interpreter', 'LaTeX');
ylabel('$$F$$, $$-G$$','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 0,'Units', 'Normalize','Position', ylabelpos,'Interpreter', 'LaTeX');



N = 1000;
k = linspace(0,1,N)';

[F,G] = LiftDefFac(k);

C = F + i*G;
C_mag = abs(C);
C_phi = angle(C);

% Reproducing lift deficiency factor plot from Theodorsen 1935 (figure 4, pg 418).
FontSizeAx = 24;
FontSizeLb = 32;
afFigurePosition = [0 2 30 20];
axespos = [0.175 0.15 0.775 0.8];
ylabelpos = [-0.125 0.5];
xlabelpos = [0.5 -0.1];

figure
set(gcf, 'Units', 'centimeters','PaperPositionMode', 'auto','Position', afFigurePosition);
set(gcf,'DefaultAxesFontSize',FontSizeAx,'DefaultAxesFontName','TimesNewRoman','DefaultAxesGridLineStyle','-.','DefaultAxesLineWidth',1,'DefaultAxesFontWeight','Normal')


hold on
    plot(k,C_mag,'-b','linewidth',2)
    plot(k,-C_phi,'-k','linewidth',2)
hold off

axis([0 1 0 1])
grid on
set(gca, 'FontName', 'TimesNewRoman', 'FontSize', FontSizeAx)
set(gca, 'Units', 'normalized', 'Position', axespos);
legend('C_{mag}','-C_{\Phi}','Location','East')
xlabel('$$k$$','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 0,'Units', 'Normalize','Position', xlabelpos,'Interpreter', 'LaTeX');
ylabel('$$C_{mag}$$, $$-C_{\Phi}$$','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 0,'Units', 'Normalize','Position', ylabelpos,'Interpreter', 'LaTeX');
