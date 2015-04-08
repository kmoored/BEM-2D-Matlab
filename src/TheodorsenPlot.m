clear
close all
clc

rho = 1;
U = 1;
f = 1;
k = pi*0.2;
St = 0.05/pi/U;
c = k*U/pi/f;
h_0 = 0.05/2/pi/f;

[Cl_Theo,Cl_c_Theo,Cl_a_Theo,Theo_t] = Theodorsen(k,St,1e4);

% Plotting Lift
FontSizeAx = 24;
FontSizeLb = 32;
afFigurePosition = [0 2 30 20];
axespos = [0.175 0.15 0.775 0.8];
ylabelpos = [-0.125 0.5];
xlabelpos = [0.5 -0.1];

figure
set(gcf, 'Units', 'centimeters','PaperPositionMode', 'auto','Position', afFigurePosition);
set(gcf,'DefaultAxesFontSize',FontSizeAx,'DefaultAxesFontName','TimesNewRoman','DefaultAxesGridLineStyle','-.','DefaultAxesLineWidth',2,'DefaultAxesFontWeight','Normal')

hold on
   
    plot(Theo_t,Cl_c_Theo*(1/2*rho*c*U^2),'-b','linewidth',3)
    plot(Theo_t,Cl_a_Theo*(1/2*rho*c*U^2),'-r','linewidth',3)
    plot(Theo_t,Cl_Theo*(1/2*rho*c*U^2),'-k','linewidth',3)

hold off

% axis([0 1 -0.04 0.04])
set(gca, 'FontName', 'TimesNewRoman', 'FontSize', FontSizeAx)
set(gca, 'Units', 'normalized', 'Position', axespos);
% legend('Steady','Unsteady','Total','Location','NorthWest')
xlabel('$$t/T$$','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 0,'Units', 'Normalize','Position', xlabelpos,'Interpreter', 'LaTeX');
ylabel('$$C_l$$','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 0,'Units', 'Normalize','Position', ylabelpos,'Interpreter', 'LaTeX');