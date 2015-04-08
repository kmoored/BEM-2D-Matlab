clear
clc

syms s_p1 s_p2 s_m1 s_m2 mu_0 mu_p1 mu_p2 mu_m1 mu_m2 s

fac = 0;

%% Second-order central difference
A = [fac  fac*(s_p1 + s_p2) fac*1/2*(s_p1 + s_p2)^2  fac*1/6*(s_p1 + s_p2)^3;...
     1     (s_p1)          1/2*(s_p1)^2        1/6*(s_p1)^3;...
     1    -(s_m1)          1/2*(s_m1)^2       -1/6*(s_m1)^3;...
     fac -fac*(s_m1 + s_m2) fac*1/2*(s_m1 + s_m2)^2 -fac*1/6*(s_m1 + s_m2)^3];
 
b = [fac*mu_p2; mu_p1; mu_m1; fac*mu_m2];

s_p2 = s;
s_p1 = s;
s_m1 = s;
s_m2 = s;

phi = simple(A\b); 
dmu_ds = simple(subs(phi(2)))
pretty(dmu_ds)


%% First-order central difference
A = [1     (s_p1);...
     1    -(s_m1)];
 
b = [mu_p1; mu_m1];

s_p2 = s;
s_p1 = s;
s_m1 = s;
s_m2 = s;

phi = simple(A\b); 
dmu_ds = simple(subs(phi(2)))
pretty(dmu_ds)

%% First-order backward difference
% A = [1   -(s_m1)     ;...
%      1 -(s_m1 + s_m2)  ];
%  
% b = [mu_m1; mu_m2];
% 
% % s_p2 = s;
% % s_p1 = s;
% % s_m1 = s;
% % s_m2 = s;
% 
% phi = simple(A\b); 
% dmu_ds = simple(subs(phi(2)))
% pretty(dmu_ds)