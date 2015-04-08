matlabpool open local 8

matlabpool close force local

clear
clc

matlabpool open local 8

% Computational parameters.
N_vec = 150;
Nstep_vec = 150;
TEfac_vec = 1;
val_vec = 5e-3;

% Physical parameters.
% d_c_vec = [100 1 1/2 1/4]';
% St_vec = [0.1 0.25 0.5]';
d_c_vec = 0;
St_vec = 0.5;


% Initializing matrices.
Ct_avg = zeros(length(N_vec),1);
Cl_avg = zeros(length(N_vec),1);
Cpow_avg = zeros(length(N_vec),1);
np = zeros(length(N_vec),1);


% for i = 1:length(Nstep_vec)
parfor j = 1:length(N_vec)
    d_c = d_c_vec
    St = St_vec
    Npanels = N_vec(j)
    Nstep = Nstep_vec
    TEfac = TEfac_vec
    val = val_vec

    [Ct_avgt,Cl_avgt,Cpow_avgt,npt,~,~] = PanelMethod2D_Comparison(d_c,St,Npanels,Nstep,TEfac,val);
    Ct_avg(j) = Ct_avgt; 
    Cl_avg(j) = Cl_avgt; 
    Cpow_avg(j) = Cpow_avgt; 
    np(j) = npt; 
end
% end

matlabpool close force local

save(['Data2D_convergence_N_vary_Nstep80_Viscous.mat'],'-v7.3');

exit