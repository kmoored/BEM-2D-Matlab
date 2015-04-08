clear
clc

% matlabpool open local 8

d_c_vec = [10 2 1 0.9 0.85 0.8 0.75 0.7 0.65 0.6 0.55 0.5 0.45 0.4 0.35 0.3 0.25 0.2 0.15]';
St_vec = [0.25]';

Ct_avg = zeros(length(St_vec),length(d_c_vec));
Cl_avg = zeros(length(St_vec),length(d_c_vec));
Cpow_avg = zeros(length(St_vec),length(d_c_vec));
np = zeros(length(St_vec),length(d_c_vec));

for i = 1:length(d_c_vec)
    for j = 1:length(St_vec)
        d_c = d_c_vec(i)
        St = St_vec(j)
        [Ct_avgtemp,Cl_avgtemp,Cpow_avgtemp,nptemp,Cl_endtemp,Ct_endtemp] = PanelMethod2DGrdSweep(d_c,St);
        Ct_avg(j,i) = Ct_avgtemp; 
        Cl_avg(j,i) = Cl_avgtemp; 
        Ct_end(j,i) = Ct_endtemp; 
        Cl_end(j,i) = Cl_endtemp; 
        Cpow_avg(j,i) = Cpow_avgtemp; 
        np(j,i) = nptemp; 
    end
end

% matlabpool close

save('SteadyGrdData.mat','-v7.3')
