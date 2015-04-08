clear
close all
clc

mu = 1;
sigma = 0;
xp = [-1/2 1/2];
zp = [0 0];
t = -[1 0];
n = -[0 1];

z = linspace(-10,10,100);
% z = 1;
x = 0*z;

[dPhi_s,dPhi_d,dL] = Phi_Test(mu,sigma,x,z,xp,zp,t,n);
dL

figure
plot(z,dPhi_d,'-k','linewidth',2)
xlabel('z')
ylabel('$$\Phi$$','interpreter','latex')

