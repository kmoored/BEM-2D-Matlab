function ddeltadx = BL_ode(x,delta,Ue,dUedx,s,H,nu)

f = 3.286*dUedx./Ue;
f = interp1(s,f,x);

g = 36/7*(0.3*exp(-1.33*H))./(log10(7/72*Ue/nu.*delta)).^(1.74 + 0.31*H);
g = interp1(s,g,x);

ddeltadx = -f.*delta + g; 




