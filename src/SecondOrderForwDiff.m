function [Qp] = SecondOrderForwDiff(mu,dL,Pan)

s_p2 = 1/2*dL(Pan(3)) + 1/2*dL(Pan(2));
s_p1 = 1/2*dL(Pan(2)) + 1/2*dL(Pan(1));

mu_p2 = mu(Pan(3));
mu_p1 = mu(Pan(2));
mu_0 = mu(Pan(1));

% Second-order backward difference
A = [ (s_p1 + s_p2) 1/2*(s_p1 + s_p2)^2;...
        (s_p1)          1/2*(s_p1)^2];
 
b = [mu_p2 - mu_0; mu_p1 - mu_0];

phi = A\b; 
dmu_ds = phi(1);
Qp = -dmu_ds;