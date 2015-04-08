matlabpool open local 8

matlabpool close force local

clear
clc

matlabpool open local 8

% Computational parameters.
N_vec = 150;
Nstep_vec = 150;
Ncyc = 8;
Nlump = 4;
TEfac_vec = 1;
val_vec = 2.5e-3;

% d_c = [1000 10 5 2 1 5/6 4/6 3/6 7/16 6/16 1/3 5/16 5/18 1/4]'

% Physical parameters.
d_c_vec = [1000 10 5 2 1 5/6 4/6 3/6 7/16 6/16]';
St_vec = [0.025:0.025:0.5]';
% d_c_vec = [7/16 6/16]'; 
% d_c_vec = [1/3 0.3 0.27 0.26 1/4]'; 
% St_vec = [0.025:0.025:0.35]';
% d_c_vec = 0;
% St_vec = [0.15:0.05:0.5]';

% Initializing matrices.
Ct_avg = zeros(length(St_vec),length(d_c_vec));
Cl_avg = zeros(length(St_vec),length(d_c_vec));
Cpow_avg = zeros(length(St_vec),length(d_c_vec));
np = zeros(length(St_vec),length(d_c_vec));
Ec = zeros(length(St_vec),length(d_c_vec));


for i = 1:length(d_c_vec)
    parfor j = 1:length(St_vec)
        d_c = d_c_vec(i)
        St = St_vec(j)
        Npanels = N_vec
        Nstep = Nstep_vec
        TEfac = TEfac_vec
        val = val_vec

        [Ct_avgt,Cl_avgt,Cpow_avgt,npt,Ect] = PanelMethod2D_v4(d_c,St,Npanels,Nstep,Ncyc,Nlump,TEfac,val); 
        Ct_avg(j,i) = Ct_avgt;
        Cl_avg(j,i) = Cl_avgt; 
        Cpow_avg(j,i) = Cpow_avgt; 
        np(j,i) = npt; 
        Ec(j,i) = Ect;
    end
end

matlabpool close force local

% save(['Data2D_Grd_N150_Nstep150_Viscous_4_refined2.mat'],'-v7.3');

exit
