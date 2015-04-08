% Defining the lift deficiency factor C(k) = F(k) + i G(k).
% k = pi*f*c/U.

function [F,G] = LiftDefFac(k)

% J: Bessel function of the first kind.
% Y: Bessel function of the second kind.
J0 = besselj(0,k);
J1 = besselj(1,k);

Y0 = bessely(0,k);
Y1 = bessely(1,k);

% F(k) = ( J1(J1 + Y0) + Y1(Y1 - J0) )/ ((J1 + Y0)^2 + (Y1 - J0)^2)
% G(k) = - ( Y1*Y0 + J1*J0 )/ ((J1 + Y0)^2 + (Y1 - J0)^2)
F = ( J1.*(J1 + Y0) + Y1.*(Y1 - J0) )./ ((J1 + Y0).^2 + (Y1 - J0).^2);
G = - ( Y1.*Y0 + J1.*J0 )./ ((J1 + Y0).^2 + (Y1 - J0).^2);
