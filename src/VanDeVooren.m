% Written by Keith Moored, 12/14/10
%
%
% The values set in this file are for a 15% thick Van de Vooren airfoil.

function [xt,zt,xb,zb,a,k,epsilon,theta1,theta2,theta,r1,r2] = VanDeVooren()

tau = 12.4*(pi/180);        % Trailing edge angle.  This changes with the thickness.
epsilon = 0.075;            % Thickness parameter.
npts=101;
c=5;

k = 2 - tau/pi;
a = c*(1 + epsilon)^(k - 1)*2^-k

theta = linspace(0,pi,(npts+1)/2)'
theta1 = atan(a*sin(theta)./(a*cos(theta) - a)) + pi;
theta2 = atan2(a*sin(theta),(a*cos(theta) - epsilon*a));

r1 = sqrt((a*cos(theta) - a).^2 + a^2*sin(theta).^2);
r2 = sqrt((a*cos(theta) - epsilon*a).^2 + a^2*sin(theta).^2);

xb = r1.^k./r2.^(k-1).*(cos(k*theta1).*cos((k - 1).*theta2) + sin(k*theta1).*sin((k - 1).*theta2)) +c;
zb = -r1.^k./r2.^(k-1).*(sin(k*theta1).*cos((k - 1).*theta2) - cos(k*theta1).*sin((k - 1).*theta2));

xb(1) = c;
zb(1) = 0;

theta = linspace(pi,0,(npts+1)/2)';
theta1 = atan(a*sin(theta)./(a*cos(theta) - a)) + pi;
theta2 = atan2(a*sin(theta),(a*cos(theta) - epsilon*a));

r1 = sqrt((a*cos(theta) - a).^2 + a^2*sin(theta).^2);
r2 = sqrt((a*cos(theta) - epsilon*a).^2 + a^2*sin(theta).^2);

xt = r1.^k./r2.^(k-1).*(cos(k*theta1).*cos((k - 1).*theta2) + sin(k*theta1).*sin((k - 1).*theta2)) +c;
zt = r1.^k./r2.^(k-1).*(sin(k*theta1).*cos((k - 1).*theta2) - cos(k*theta1).*sin((k - 1).*theta2));

xt(end) = c;
zt(end) = 0;





