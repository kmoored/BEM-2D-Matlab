if outerCorr <= 1
    fprintf('===============================================================================\n')
    fprintf(' TIME-STEP NUMBER = %i, FLOW TIME = %f\n',i_t-1,t)
    fprintf('-------------------------------------------------------------------------------\n')
end

%% Solving Rigid Body Fluid Flow Simulation
PanelMethod2D_v9_rev_wcs211
xpOld = [xp-delxp xpOld(:,1) xpOld(:,2) xpOld(:,3) xpOld(:,4) xpOld(:,5)];
zpOld = [zp-delzp zpOld(:,1) zpOld(:,2) zpOld(:,3) zpOld(:,4) zpOld(:,5)];

% Calculating the new _0 position to plug into the kinematic function
%[nodes_new(:,1),nodes_new(:,2)] = antiKinematics_HeavePitch2D(nodes(:,1),nodes(:,2),alpha_max,h_c,f,t,phi,'solid');
[ nodes_new(:,1),nodes_new(:,2) ] = antiKinematics_ZeroAoA( nodes(:,1),nodes(:,2),h_c,c,f,t,phi,Qinf,'solid',ramped,i_t );

%% Updating the kinematics (rigid body) for the structual mesh
% Updating to the new kinematics
%nodes(:,1) = (nodes_new(:,1) - nodes_new(1,1))*cos(alpha_max*sin(2*pi*f*t + phi));
%nodes(:,2) = h_c*c*sin(2*pi*f*t) + (nodes_new(:,1) - nodes_new(1,1))*sin(alpha_max*sin(2*pi*f*t + phi));
nodes(:,1) = (nodes_new(:,1) - nodes_new(1,1))*cos(atan(-2*pi*ramped(i_t)*h_c*c*f*cos(2*pi*f*t)/Qinf));
nodes(:,2) = ramped(i_t)*h_c*c*sin(2*pi*f*t) + (nodes_new(:,1) - nodes_new(1,1))*sin(atan(-2*pi*ramped(i_t)*h_c*c*f*cos(2*pi*f*t)/Qinf));

% Calculating shift in panel positions with the swimming velocity.
nodalDelxp = x_b(2)*ones(length(nodes),1);
nodalDelzp = z_b(2)*ones(length(nodes),1);

% Superposing the kinematics and swimming translations.
nodes(:,1) = nodes(:,1) + nodalDelxp;
nodes(:,2) = nodes(:,2) + nodalDelzp;

%% Calculate Wake Roolup
calcWakeRollup

%% Calculate Free-swimming
calcFreeSwimming

Cf = norm(F)/(1/2*rho*Qinf^2*c*b);
Cl = L/(1/2*rho*Qinf^2*c*b);
Ct = T/(1/2*rho*Qinf^2*c*b);
Cpow = Pow/(1/2*rho*Qinf^3*c*b);

fprintf('| Solution Information:                                                       |\n')
fprintf('|     D_visc   = %13e                                                |\n',D_visc)
fprintf('|     Cf       = %13e                                                |\n',Cf    )
fprintf('|     Cl       = %13e                                                |\n',Cl    )
fprintf('|     Ct       = %13e                                                |\n',Ct    )
fprintf('|     Cpow     = %13e                                                |\n',Cpow  )
fprintf('|     Gamma    = %13e                                                |\n',Gamma )

if i_t >= (Ncyc-1)*Nstep + 1
    Cl_avg = Cl_avg + Cl/Nstep;
    Ct_avg = Ct_avg + Ct/Nstep;
    Tnet_avg = Tnet_avg + T/Nstep;
    D_avg = D_avg + D_visc/Nstep;
    Pow_avg = Pow_avg + Pow/Nstep;
    Cpow_avg = Cpow_avg + Cpow/Nstep;

    fprintf('|                                                                             |\n')
    fprintf('| Solution Average Cycle Information:                                         |\n')
    fprintf('|     Cl_avg   = %13e                                                |\n',Cl_avg  )
    fprintf('|     Ct_avg   = %13e                                                |\n',Ct_avg  )
    fprintf('|     Tnet_avg = %13e                                                |\n',Tnet_avg)
    fprintf('|     D_avg    = %13e                                                |\n',D_avg   )
    fprintf('|     Pow_avg  = %13e                                                |\n',Pow_avg )
    fprintf('|     Cpow_avg = %13e                                                |\n',Cpow_avg)
end

fprintf('|     %3i%% Complete                                                           |\n',round((i_t/(Ncyc*Nstep+1))*100))
fprintf('+-----------------------------------------------------------------------------+\n\n')


%% Plotting
newPlotTimeStepFig
printConvergenceInfo
