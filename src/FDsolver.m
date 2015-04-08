clear
close all
clc

syms dt Q Q_p1 Q_m1 Q_m2 Q_m3

val = 3;

A = [   dt 1/2*dt^2  1/6*dt^3   1/24*dt^4;...
       -dt 1/2*dt^2 -1/6*dt^3   1/24*dt^4;...
     -2*dt   2*dt^2 -4/3*dt^3    2/3*dt^4;...
     -3*dt 9/2*dt^2 -9/2*dt^3   27/8*dt^4];
 
b = [Q_p1 - Q; Q_m1 - Q; Q_m2 - Q; Q_m3 - Q];

x = simple(A(1:val,1:val)\b(1:val));
dQ_dt = x(1)