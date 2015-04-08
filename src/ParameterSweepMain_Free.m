% Written by: Keith Moored, 1/27/11
%
% This program runs the three-dimnensional unsteady panel code over a range
% of parameters.

% Closing parallel CPUs.
% matlabpool close force local

clear
close all
clc

% % Opening parallel CPUs.
% matlabpool open local 2
% 
% % Closing parallel CPUs.
% matlabpool close force local
% 
% % Opening parallel CPUs.
% matlabpool open local 2
 
k_vec = linspace(0.1,20,25)';
A_c_vec = [0.125 0.25 0.375 0.5 0.625 0.75]';
Qinf_vec = [1 2 3 4 5]';

M_star = 1/5;
A_bp = 1;
c = 1;
Npanels = 150;
Nstep = 150;
Ncyc = 8;
Nlump = 4;
val = 2.5e-4;
TEfac = 1;

% % A_bp_vec = [1:0.25:10 50]';
% A_bp_vec = [1:0.5:10 50]';
% A_c_vec = [0.25 0.5 0.75]';

% M_star_vec = [0.5 1 2 4]';
% A_bp_vec = 5;
% f_vec = linspace(0.5,8,20)';
% A_c_vec = [0.25]';


Ct_avg = zeros(length(k_vec),length(A_c_vec),length(Qinf_vec));
Tnet_avg = zeros(length(k_vec),length(A_c_vec),length(Qinf_vec));
D_avg = zeros(length(k_vec),length(A_c_vec),length(Qinf_vec));
Cl_avg = zeros(length(k_vec),length(A_c_vec),length(Qinf_vec));
Cpow_avg = zeros(length(k_vec),length(A_c_vec),length(Qinf_vec));
Pow_avg = zeros(length(k_vec),length(A_c_vec),length(Qinf_vec));
np = zeros(length(k_vec),length(A_c_vec),length(Qinf_vec));
Ec = zeros(length(k_vec),length(A_c_vec),length(Qinf_vec));
Ucruise = zeros(length(k_vec),length(A_c_vec),length(Qinf_vec));
Q0 = zeros(length(k_vec),length(A_c_vec),length(Qinf_vec),Ncyc*Nstep+1);


for i = 1:length(Qinf_vec)   
    Qinf = Qinf_vec(i);
    
    for j = 1:length(A_c_vec)
        A_c = A_c_vec(j);

        parfor q = 1:length(k_vec)
            k = k_vec(q);

            [Ct_avgt,Tnet_avgt,D_avgt,Cl_avgt,Cpow_avgt,Pow_avgt,npt,Ect,~,~,rhot,~,Ucruiset,Q0t] = PanelMethod2D_v6(A_bp,M_star,k,A_c,c,Qinf,Npanels,Nstep,Ncyc,Nlump,TEfac,val);

            Ct_avg(q,j,i) = Ct_avgt;
            Tnet_avg(q,j,i) = Tnet_avgt;
            D_avg(q,j,i) = D_avgt;
            Cl_avg(q,j,i) = Cl_avgt;
            Cpow_avg(q,j,i) = Cpow_avgt;
            Pow_avg(q,j,i) = Pow_avgt;
            np(q,j,i) = npt;
            Ec(q,j,i) = Ect;
            rho = rhot;
            Ucruise(q,j,i) = Ucruiset;
            Q0(q,j,i,:) = Q0t;
        end
    end
end

% Closing parallel CPUs.
matlabpool close force local

%% Save Data.

save('FixedV/FixedV_Pitch_Data.mat','-v7.3')

% exit        

