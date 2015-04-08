clear
close all
clc
c=1.0;
tmax = 0.15; % percent thickness 

% Analytical solution.

[xpt,zpt,xpb,zpb,a_vdv,k,epsilon,theta1,theta2,theta,r1,r2] = VanDeVooren(c(1),200,tmax/2);

xp_a = [xpb;xpt(2:end)];
zp_a = [zpb(1:end-1)-zpb(1);0;zpt(2:end)-zpt(end)];

Aa = cos((k-1)*theta1).*cos(k*theta2) + sin((k-1)*theta1).*sin(k*theta2);
Ba = sin((k-1)*theta1).*cos(k*theta2) - cos((k-1)*theta1).*sin(k*theta2);
D0 = a_vdv*(1 - k + k*epsilon);
D1 = Aa.*(a_vdv*cos(theta) - D0) - Ba.*a_vdv.*sin(theta);
D2 = Aa.*a_vdv.*sin(theta) + Ba.*(a_vdv.*cos(theta) - D0);

u_anltc = 2*Qinf*r2.^k./r1.^(k-1).*(sin(alpha) - sin(alpha - theta))./(D1.^2 + D2.^2).*(D1.*sin(theta) + D2.*cos(theta));
w_anltc = -2*Qinf*r2.^k./r1.^(k-1).*(sin(alpha) - sin(alpha - theta))./(D1.^2 + D2.^2).*(D1.*cos(theta) - D2.*sin(theta));

Cp_anltc_t = 1 - (u_anltc.^2 + w_anltc.^2)/Qinf^2;

alpha_b = -alpha;
[xpt,zpt,xpb,zpb,a_vdv,k,epsilon,theta1,theta2,theta,r1,r2] = VanDeVooren(c(1),200,tmax_f/2);

xp_a = [xpb;xpt(2:end)];
zp_a = [zpb(1:end-1)-zpb(1);0;zpt(2:end)-zpt(end)];

Aa = cos((k-1)*theta1).*cos(k*theta2) + sin((k-1)*theta1).*sin(k*theta2);
Ba = sin((k-1)*theta1).*cos(k*theta2) - cos((k-1)*theta1).*sin(k*theta2);
D0 = a_vdv*(1 - k + k*epsilon);
D1 = Aa.*(a_vdv*cos(theta) - D0) - Ba.*a_vdv.*sin(theta);
D2 = Aa.*a_vdv.*sin(theta) + Ba.*(a_vdv.*cos(theta) - D0);

u_anltc = 2*Qinf*r2.^k./r1.^(k-1).*(sin(alpha_b) - sin(alpha_b - theta))./(D1.^2 + D2.^2).*(D1.*sin(theta) + D2.*cos(theta));
w_anltc = -2*Qinf*r2.^k./r1.^(k-1).*(sin(alpha_b) - sin(alpha_b - theta))./(D1.^2 + D2.^2).*(D1.*cos(theta) - D2.*sin(theta));

Cp_anltc_b = 1 - (u_anltc.^2 + w_anltc.^2)/Qinf^2;
