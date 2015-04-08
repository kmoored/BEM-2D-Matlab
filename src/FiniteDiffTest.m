clear
clc

N = 9
delt = 1/(N-1)
A = 1;

t = linspace(0,1,N)';
t_ref = linspace(0,1,1000)';
x = A/2/pi*sin(2*pi*t);
v = A*cos(2*pi*t_ref);

% Finite difference approx of velocity.
v_3 = (x(3:end) - x(1:end-2))./(t(3:end) - t(1:end-2));
v_5 = (1/3*x(1:end-4) - 8/3*x(2:end-3) + 8/3*x(4:end-1) - 1/3*x(5:end))./(t(5:end) - t(1:end-4));



% plotting
FontSizeAx = 24;
FontSizeLb = 32;
afFigurePosition = [1 1 20 15];
ylabelpos = [-0.1 0.5];
xlabelpos = [0.5 -0.15];
axespos = [0.15 0.2 0.8 0.75];

figure
set(gcf, 'Units', 'centimeters','PaperPositionMode', 'auto','Position', afFigurePosition);
set(gcf,'DefaultAxesFontSize',FontSizeAx,'DefaultAxesFontName','TimesNewRoman','DefaultAxesGridLineStyle','-.','DefaultAxesLineWidth',2,'DefaultAxesFontWeight','Normal')


hold on
    plot(t_ref,v,'-k','linewidth',3)
    plot(t(2:end-1),v_3,'ob','linewidth',3,'markersize',8)
    plot(t(3:end-2),v_5,'rx','linewidth',3,'markersize',8)
hold off

axis([0 1 -1.5*A 1.5*A])
set(gca, 'FontName', 'TimesNewRoman', 'FontSize', FontSizeAx)
set(gca, 'Units', 'normalized', 'Position', axespos);
xlabel('$$t/T$$','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 0,'Units', 'Normalize','Position', xlabelpos,'Interpreter', 'LaTeX');
ylabel('$$v$$','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 90,'Units', 'Normalize','Position', ylabelpos,'Interpreter', 'LaTeX');
