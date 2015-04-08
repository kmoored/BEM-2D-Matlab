clear 
close all
clc

syms x0 x1 x2 x z sig 

Phi = sig/2/pi*int(log(sqrt((x - x0)^2 + z^2)),x0,x1,x2)