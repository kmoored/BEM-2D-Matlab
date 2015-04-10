fprintf('+-----------------------------------------------------------------------------+\n')
if (fsiResidualNorm <= outerCorrTolerance)
    fprintf('| SOLUTION CONVERGED!                                                         |\n')
    fprintf('+-----------------------------------------------------------------------------+\n')
else
    fprintf('| WARNING! MAX INNER-LOOP ITERATIONS REACHED                                  |\n')
    fprintf('+-----------------------------------------------------------------------------+\n')
end

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
