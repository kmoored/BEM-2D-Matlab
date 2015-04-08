function [Qp] = SecondOrderCenDiff(mu,dL,Pan)

s_p2 = 1/2*dL(Pan(5)) + 1/2*dL(Pan(4));
s_p1 = 1/2*dL(Pan(4)) + 1/2*dL(Pan(3));
s_m1 = 1/2*dL(Pan(2)) + 1/2*dL(Pan(3));
s_m2 = 1/2*dL(Pan(1)) + 1/2*dL(Pan(2));

mu_p2 = mu(Pan(5));
mu_p1 = mu(Pan(4));
mu_m1 = mu(Pan(2));
mu_m2 = mu(Pan(1));

% Second-order central difference
A = [1  (s_p1 + s_p2) 1/2*(s_p1 + s_p2)^2  1/6*(s_p1 + s_p2)^3;...
     1     (s_p1)          1/2*(s_p1)^2        1/6*(s_p1)^3;...
     1    -(s_m1)          1/2*(s_m1)^2       -1/6*(s_m1)^3;...
     1 -(s_m1 + s_m2) 1/2*(s_m1 + s_m2)^2 -1/6*(s_m1 + s_m2)^3];
 
b = [mu_p2; mu_p1; mu_m1; mu_m2];

phi = A\b; 
dmu_ds = phi(2);
Qp = -dmu_ds;
