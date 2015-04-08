clear
clc

% c = 0.12 m, Uinf = 0.06 m/s, NACA 0005, pure pitch, A/c = 0.25 (peak-to-peak)
% Case 1: d/c = infty.
% Case 2: d/c = 1/2.
% Case 3: d/c = 1/4.

DataGrd = load('SteadyGrdData');

%% Data
d_c = [10 2 1 0.9 0.85 0.8 0.75 0.7 0.65 0.6 0.55 0.5 0.45 0.4 0.35 0.3 0.25 0.2 0.15]';
d_c = d_c - 0.1305;

% Coefficient of lift
Cl = DataGrd.Cl_end';
Cl = Cl/Cl(1);

%% Experimental Data
		

%% Plotting

% Plotting coefficient of thrust.
figure
set(gcf,'DefaultAxesfontsize',20,'DefaultAxesfontname','TimesNewRoman','DefaultAxesGridLineStyle','-.')

hold on
plot(d_c,Cl,'-k','linewidth',2)
plot([0 10],[1 1],'-k','linewidth',1)


xlabel('$$h_c$$','interpreter','latex','fontname','TimesNewRoman','fontsize',20)
ylabel('$$C_l$$','interpreter','latex','fontname','TimesNewRoman','fontsize',20,'Rotation',0)

axis([0 1.5 0 2])



