function [Qp] = SecondOrderCenForwDiff(mu,dL,Pan)

s_p2 = 1/2*dL(Pan(4)) + 1/2*dL(Pan(3));
s_p1 = 1/2*dL(Pan(3)) + 1/2*dL(Pan(2));
s_m1 = 1/2*dL(Pan(1)) + 1/2*dL(Pan(2));


mu_p2 = mu(Pan(4));
mu_p1 = mu(Pan(3));
mu_0 = mu(Pan(2));
mu_m1 = mu(Pan(1));

% Second-order central difference
A = [(s_p1 + s_p2) 1/2*(s_p1 + s_p2)^2  1/6*(s_p1 + s_p2)^3;...
        (s_p1)          1/2*(s_p1)^2        1/6*(s_p1)^3;...
        -(s_m1)          1/2*(s_m1)^2       -1/6*(s_m1)^3];
 
b = [mu_p2 - mu_0; mu_p1 - mu_0; mu_m1 - mu_0];

phi = A\b; 
dmu_ds = phi(1);
Qp = -dmu_ds;