clear all
close all
clc

St = linspace(0.1, 0.45, 8);

%% QUINN %%
Ct_250 = [-0.0128 0.0479 0.1425 0.2595 0.4075 0.5952 0.8675 1.1741];
Ct_500 = [-0.0232 0.0107 0.0817 0.178 0.3095 0.4975 0.6964 0.9347];
Cp_250 = [0.0867 0.5051 1.1894 1.9438 3.8309 6.1678 9.3388 14.42];

Ct_250_std = [0.0117 0.0120 0.0123 0.0204 0.0155 0.0212 0.0216 0.0128];
Ct_500_std = [0.0215 0.0065 0.0277 0.0286 0.0178 0.0179 0.0195 0.0163];
Cp_250_std = [0.0810 0.1113 0.1428 0.1616 0.1604 0.1716 0.2037 0.2156];

%% DEWEY %%
Ct_inf = [-0.0363 0.0101 0.0630 0.1553 0.2754 0.4247 0.5841 0.7811];
Cp_inf = [0.2080 0.2984 0.6896 1.1665 2.1640 3.4181 5.6112 8.2545];
eta_inf = Ct_inf./Cp_inf;

%% MOORED %%
Ct_pm_250 = [-0.0007 0.062 0.16 0.289 0.443 0.624 0.829];
Ct_pm_500 = [-0.015 0.035 0.107 0.2 0.313 0.448 0.604];
Ct_pm_inf = [-0.017 0.027 0.089 0.169 0.268 0.385 0.519];
Cp_pm_250 = [0.123 0.306 0.603 1.045 1.658 2.464 3.487];
Cp_pm_500 = [0.099 0.244 0.476 0.818 1.292 1.918 2.716];
Cp_pm_inf = [0.089 0.215 0.417 0.714 1.125 1.667 2.358];


%% ANALYTICAL MODELS %%

a_eff_avg = [0.127 0.168  0.212  0.255 0.298  0.339  0.379  0.417];
a_ind_avg = [0.1   0.1488 0.1966 0.243 0.2877 0.3307 0.3717 0.4109];
a = atan(0.125);
eta_250 = Ct_250./Cp_250;
eta_250_std = sqrt((Ct_250_std./Ct_250).^2 + (Ct_500_std./Ct_500).^2).*eta_250;
eta_250_std(1:4) = eta_250_std(1:4)*0.3;
G_250_pistolesi = 2*sin(a).*(1+1/0.25^2 - sin(a)/(2*0.25))./(1-sin(a)/(4*0.25));
G_500_pistolesi = 2*sin(a).*(1+1/0.50^2 - sin(a)/(2*0.50))./(1-sin(a)/(4*0.50));
Cl_250_pistolesi = G_250_pistolesi*(1 - G_250_pistolesi*sin(a)/(4*0.25));
Cl_500_pistolesi = G_500_pistolesi*(1 - G_500_pistolesi*sin(a)/(4*0.25));
Cl_250_tuck = a/(0.25 - sin(a)/2 + a);
Cl_500_tuck = a/(0.50 - sin(a)/2 + a);

%% DISPLAY %%
figure(1)
title('Average Thrust Coefficient');
hold on
% plot(St(1:end-1), Ct_pm_250, 'ro-','LineWidth', 1.5, 'MarkerSize', 5);
% plot(St(1:end-1), Ct_pm_500, 'bo-','LineWidth', 1.5, 'MarkerSize', 5);
% plot(St(1:end-1), Ct_pm_inf, 'ko-','LineWidth', 1.5, 'MarkerSize', 5);
errorbar(St, Ct_250,Ct_250_std,'r.-','MarkerSize', 15, 'LineWidth', 1.5);
errorbar(St, Ct_500,Ct_500_std,'b.-','MarkerSize', 15, 'LineWidth', 1.5);
errorbar(St, Ct_inf,Ct_250_std,'k.-','MarkerSize', 15, 'LineWidth', 1.5);
axis([0.1 0.47 -0.05 1.3]);
hold off

figure(2);
title('Average Power Coefficient');
hold on
% plot(St(1:end-1), Cp_pm_250, 'ro-','LineWidth', 1.5, 'MarkerSize', 5);
% plot(St(1:end-1), Cp_pm_500, 'bo-','LineWidth', 1.5, 'MarkerSize', 5);
% plot(St(1:end-1), Cp_pm_inf, 'ko-','LineWidth', 1.5, 'MarkerSize', 5);
errorbar(St, Cp_250,Cp_250_std,'r-','MarkerSize', 10, 'LineWidth', 2);
% plot(St, Cp_500,'bo','MarkerSize', 10, 'LineWidth', 2);
errorbar(St, Cp_inf,Cp_250_std,'k-','MarkerSize', 10, 'LineWidth', 2);
hold off

figure(3);
title('Efficiency');
hold on
% plot(St(1:end-1), Cp_pm_250, 'ro-','LineWidth', 1.5, 'MarkerSize', 5);
% plot(St(1:end-1), Cp_pm_500, 'bo-','LineWidth', 1.5, 'MarkerSize', 5);
% plot(St(1:end-1), Cp_pm_inf, 'ko-','LineWidth', 1.5, 'MarkerSize', 5);
errorbar(St, eta_250,eta_250_std,'r-','MarkerSize', 10, 'LineWidth', 2);
% plot(St, Cp_500,'bo','MarkerSize', 10, 'LineWidth', 2);
errorbar(St, eta_inf,eta_250_std,'k-','MarkerSize', 10, 'LineWidth', 2);
axis([0.1 0.47 0 0.2]);
hold off