% Written by: Keith Moored, 1/27/11
%
% This program runs the three-dimnensional unsteady panel code over a range
% of parameters.

% Closing parallel CPUs.
% matlabpool close force local

clear
close all
clc

% Opening parallel CPUs.
matlabpool open local 8

% Closing parallel CPUs.
matlabpool close force local

% Opening parallel CPUs.
matlabpool open local 8
 
d_c = 0;
Re_vec = [1e4]';
St_vec = [0.05:0.025:0.5]';
A_c_vec = [0.0625 0.125 0.25 0.375 0.5 0.675]';
Npanels = 150;
Nstep = 150;
Ncyc = 8;
Nlump = 4;
TEfac = 1;
val = 2.5e-3;

Ct_avg = zeros(length(St_vec),length(A_c_vec),length(Re_vec));
Tnet_avg = zeros(length(St_vec),length(A_c_vec),length(Re_vec));
D_avg = zeros(length(St_vec),length(A_c_vec),length(Re_vec));
Cl_avg = zeros(length(St_vec),length(A_c_vec),length(Re_vec));
Cpow_avg = zeros(length(St_vec),length(A_c_vec),length(Re_vec));
np = zeros(length(St_vec),length(A_c_vec),length(Re_vec));
Ec = zeros(length(St_vec),length(A_c_vec),length(Re_vec));

for i = 1:length(Re_vec)
    
    Re = Re_vec(i);
    
    for k = 1:length(A_c_vec)

        A_c = A_c_vec(k);

        parfor j = 1:length(St_vec)

            St = St_vec(j);  

            [Ct_avgt,Tnet_avgt,D_avgt,Cl_avgt,Cpow_avgt,npt,Ect,ft,A_TEt,rhot,ct,Qinft] = PanelMethod2D_v4(d_c,St,A_c,Re,Npanels,Nstep,Ncyc,Nlump,TEfac,val);

            Ct_avg(j,k,i) = Ct_avgt;
            Tnet_avg(j,k,i) = Tnet_avgt;
            D_avg(j,k,i) = D_avgt;
            Cl_avg(j,k,i) = Cl_avgt;
            Cpow_avg(j,k,i) = Cpow_avgt;
            np(j,k,i) = npt;
            Ec(j,k,i) = Ect;
            f(j,k,i) = ft;
            A(j,k,i) = A_TEt;

            rho = rhot;
            c = ct;
            U = Qinft;
        end
    end
end

% Closing parallel CPUs.
matlabpool close force local

%% Save Data.

save('FixedV/FixedV_Pitch_Data_Re_1e4.mat','-v7.3')

exit        

