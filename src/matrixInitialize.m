xw = zeros(1,Nlump*Nstep+1);
zw = zeros(1,Nlump*Nstep+1);
xl = zeros(1,2);
zl = zeros(1,2);

vtTE = zeros(1,2);
vnTE = zeros(1,2);
vtw = zeros(Nlump*Nstep,2);
vnw = zeros(Nlump*Nstep,2);
vtlw = zeros(1,2);
vnlw = zeros(1,2);

sigma = zeros(Npanels,1);
mu = zeros(Npanels,3);
muTE = zeros(1,2);
muW = zeros(1,Nlump*Nstep);
muLump = zeros(1,2);
Gamma = zeros(1,1);
GammaW = zeros(1,Nlump*Nstep+1);
GammaTot = 0;

Qp = zeros(Npanels,1);

Cw = zeros(Npanels,Nlump*Nstep);
Clw = zeros(Npanels,1);
dPhi_wd = zeros(Nlump*Nstep,Npanels);
dPhi_wd_2 = zeros(Nlump*Nstep,Npanels);

dFshear = zeros(Npanels,1);

fighand = zeros(Ncyc*Nstep+1,1);
wakeInd = 0;
Tnet_avg = 0;
D_avg = 0;
Cl_avg = 0;
Ct_avg = 0;
Pow_avg = 0;
Cpow_avg = 0;
np = 0;
Ec = 0;

% Initializing matrices for the LES calcs.
if LES == 1
    xpbod = zeros(Npanels + 1,Ncyc*Nstep + 1);
    zpbod = zeros(Npanels + 1,Ncyc*Nstep + 1);

    xpt_LES = zeros(2,Ncyc*Nstep);
    zpt_LES = zeros(2,Ncyc*Nstep);
    xpb_LES = zeros(2,Ncyc*Nstep);
    zpb_LES = zeros(2,Ncyc*Nstep);
    
    vtLEt = zeros(Ncyc*Nstep,2);
    vnLEt = zeros(Ncyc*Nstep,2);
    vtLEb = zeros(Ncyc*Nstep,2);
    vnLEb = zeros(Ncyc*Nstep,2);
    
    Vshear_t = zeros(2,Ncyc*Nstep+1);
    Vshear_b = zeros(2,Ncyc*Nstep+1);
    
    muLEt = zeros(1,Ncyc*Nstep);
    muLEb = zeros(1,Ncyc*Nstep);
    
    CLEwt = zeros(Npanels,Ncyc*Nstep-1);
    CLEwb = zeros(Npanels,Ncyc*Nstep-1);
    
    Sep_t_Store = ones(1,Ncyc*Nstep+1);
    Sep_b_Store = ones(1,Ncyc*Nstep+1);
    Sep_t = 0;
    Sep_b = 0;
else
    xpbod = [];
    zpbod = [];
    
    xpt_LES = [];
    zpt_LES = [];
    xpb_LES = [];
    zpb_LES = [];
    
    vtLEt = [];
    vnLEt = [];
    vtLEb = [];
    vnLEb = [];
    
    Vshear_t = [];
    Vshear_b = [];
    
    muLEt = [];
    muLEb = [];
    
    CLEwt = [];
    CLEwb = [];
    
    Sep_t_Store = [];
    Sep_b_Store = [];
    Sep_t = 0;
    Sep_b = 0;
end

nodalDisplacements = zeros(1,3);