clear 
close all
clc

c = 1;
tmax = 0.1*c;
N = 250;

b1 = body(N,c,tmax);
b1 = TearDropShape(b1);