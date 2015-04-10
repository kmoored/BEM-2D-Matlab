% Shifting wake strengths and dumping last one.
if (i_t > 1 && outerCorr <= 1)
    muW(1,2:end) = muW(1,1:end-1);
    muW(1,1) = muTE(1,1);
end    

%% Updating the kinematics.
if outerCorr > 1
    % Superposing the structual displacements.
    xp = xp + fluidNodeDispl(:,1) - fluidNodeDisplOld(:,1);
    zp = zp + fluidNodeDispl(:,2) - fluidNodeDisplOld(:,2);
end

% Updating to the new kinematics
if outerCorr <= 1
    %[xp,zp,Vp] = Kinematics_HeavePitch2D(xp_new,zp_new,alpha_max,h_c,f,t,phi);
    [ xp,zp,Vp,pitchingAngle(i_t),heavePos(i_t) ] = Kinematics_ZeroAoA( xp_new,zp_new,h_c,c,f,t,Qinf,ramped,i_t );
end

% Calculating locations of collocation points (mid-panel on the surface).
xc = (xp(2:end) + xp(1:end-1))/2;
zc = (zp(2:end) + zp(1:end-1))/2;

xpbod = xp;
zpbod = zp;

if outerCorr <= 1
    % Shifting body coordinate system origin and dumping last numbers (all body origin location are written to data files).
    x_b(2:end) = x_b(1:end-1);
    z_b(2:end) = z_b(1:end-1);

    % Calculating shift in panel positions with the swimming velocity. 
    delxp = x_b(2)*ones(Npanels+1,1);
    delzp = z_b(2)*ones(Npanels+1,1);
    delxc = x_b(2)*ones(Npanels,1);
    delzc = z_b(2)*ones(Npanels,1);

    % Superposing the kinematics and swimming translations.
    xp = xp + delxp;
    zp = zp + delzp;
    xc = xc + delxc;
    zc = zc + delzc;
end

if i_t > 1
    Vp(:,1) = backwardsDifference(delT,xp-delxp,xpOld(:,1),xpOld(:,2),xpOld(:,3),xpOld(:,4),xpOld(:,5),xpOld(:,6),min(i_t-1,2));
    Vp(:,2) = backwardsDifference(delT,zp-delzp,zpOld(:,1),zpOld(:,2),zpOld(:,3),zpOld(:,4),zpOld(:,5),zpOld(:,6),min(i_t-1,2));
end

if debugPlots == 1
    figure(5)
    plot(xp_new,zp_new)
    axis([min(xp_new) max(xp_new) -0.5 0.5])
    drawnow
end

Vc = 1/2*Vp(2:end,:) + 1/2*Vp(1:end-1,:);  % Collocation point velocity

% Updated normal and tangential vectors for the panels at collocation points.  Each row is a new
% collocation point.
vttemp = [(xp(2:end) - xp(1:end-1)) (zp(2:end) - zp(1:end-1))];
vt = diag(1./sqrt(vttemp(:,1).^2 + vttemp(:,2).^2))*vttemp;
vn = [-vt(:,2) vt(:,1)];  

%% Calculating the new trailing edge panels.
if (i_t > 1 && outerCorr <= 1)
    xTE_old = xTE;
    zTE_old = zTE;
end
xp_avg = xp(end-1)/2 + xp(2)/2;
zp_avg = zp(end-1)/2 + zp(2)/2;
TEvec = [xp(end) - xp_avg; zp(end) - zp_avg]/norm([xp(end) - xp_avg; zp(end) - zp_avg]);
xTE = [xp(end); xp(end) + TEvec(1)*sqrt(Q0(:,i_t)'*Q0(:,i_t))*delT*TEfac];
zTE = [zp(end); zp(end) + TEvec(2)*sqrt(Q0(:,i_t)'*Q0(:,i_t))*delT*TEfac];

dLTE = sqrt((xTE(2)-xTE(1))^2 + (zTE(2)-zTE(1))^2);

vttemp = [(xTE(2) - xTE(1)) (zTE(2) - zTE(1))];
vtTE = diag(1./sqrt(vttemp(1).^2 + vttemp(2).^2))*vttemp;
vnTE = [-vtTE(2) vtTE(1)];

%% Shedding a wake panel after first time step and specifying lumped wake.
if (i_t > 1 && outerCorr <= 1)

    muLump(2) = muLump(1);
    xl(2) = xl(1);
    zl(2) = zl(1);

    % Lumped wake doublet elements 
    if i_t > Nlump*Nstep + 1
        muLump(1) = muLump(2) - GammaW(Nlump*Nstep + 1);
        GammaTot_new = GammaTot + abs(-GammaW(Nlump*Nstep + 1));

        xl(1) = abs(GammaTot/GammaTot_new)*xl(2) + abs(-GammaW(Nlump*Nstep + 1)/GammaTot_new)*xw(Nlump*Nstep + 1);
        zl(1) = abs(GammaTot/GammaTot_new)*zl(2) + abs(-GammaW(Nlump*Nstep + 1)/GammaTot_new)*zw(Nlump*Nstep + 1);

        GammaTot = GammaTot_new;

        vttemp = [(xl(1) - xw(Nlump*Nstep))' (zl(1) - zw(Nlump*Nstep))'];
        vtlw = diag(1./sqrt(vttemp(1).^2 + vttemp(2).^2))*vttemp;
        vnlw = [-vtlw(2) vtlw(1)]; 
    end

    % Shifting wake panel coordinates
    xw(1,2:end) = xw(1,1:end-1);
    zw(1,2:end) = zw(1,1:end-1);

    % Calculating the downstream position of the newly shed wake panel
    xw(1,2) = xTE_old(2);
    zw(1,2) = zTE_old(2);       
end
if i_t > 1
    % Calculating the upstream position of the newly shed wake panel
    xw(1,1) = xTE(2);
    zw(1,1) = zTE(2);

    % Updating normal and tangential vectors of the wake panels.
    vttemp = [(xw(1,2:end) - xw(1,1:end-1))' (zw(1,2:end) - zw(1,1:end-1))'];
    vtw = diag(1./sqrt(vttemp(:,1).^2 + vttemp(:,2).^2))*vttemp;
    vnw = [-vtw(:,2) vtw(:,1)];         
end

%% Calculating source strengths. 
Vt = kron(ones(Npanels,1),Q0(:,i_t)') + Vc;   % Total body velocity at the collocation points.  The fluid velocity would be -Vt when an observer is on the body.
sigma = sum(vn.*Vt,2);          % Setting the source strength from the Dirichlet formulation to satisfy the no flux condition.


%% Moving collocation points inward for the Dirichlet formulation.
InVal = abs(CptInval*zcval);
Xc = xc - InVal.*vn(:,1);
Zc = zc - InVal.*vn(:,2);


%% Calculating influence coefficients 
% The row denotes the collocation point while the column denotes the doublet element.  

% Influence of the body panels
Cb = zeros(Npanels,Npanels);
dPhi_s_2 = zeros(Npanels,Npanels);
dPhi_d_2 = zeros(Npanels,Npanels);

[dPhi_s,dPhi_d,dL] = Phi(ones(1,Npanels),ones(1,Npanels),Xc',Zc',xp',zp',vt',vn');

if grd == 1
    [dPhi_s_2,dPhi_d_2,dLtemp2] = Phi(ones(1,Npanels),ones(1,Npanels),Xc',Zc',xp',-zp',[vt(:,1) -vt(:,2)]',[vn(:,1) -vn(:,2)]');
end

B = dPhi_s' + dPhi_s_2';
Cb(:,1:Npanels) = dPhi_d' + dPhi_d_2';


% Influence of the trailing edge panel
dPhi_d_2 = zeros(1,Npanels);
Cte = zeros(Npanels,Npanels);

[~,dPhi_d,~] = Phi(1,0,Xc',Zc',xTE',zTE',vtTE',vnTE');

if grd == 1
    [~,dPhi_d_2,~] = Phi(1,0,Xc',Zc',xTE',-zTE',[vtTE(1,1) -vtTE(1,2)]',[vnTE(1,1) -vnTE(1,2)]');
end

Cte(:,Npanels) = dPhi_d' + dPhi_d_2';
Cte(:,1) = -dPhi_d' - dPhi_d_2';

% Influence of the wake panels    
 if (i_t > 1)   
    wakeInd = min([(i_t-1) Nlump*Nstep]);
%         wakeInd = min([(find(muW == 0,1)-1) Nlump*Nstep]); 
    [~,dPhi_d,~] = Phi(ones(1,length(muW(1:wakeInd))),zeros(1,length(muW(1:wakeInd))),Xc',Zc',xw(1:wakeInd+1),zw(1:wakeInd + 1),vtw(1:wakeInd,:)',vnw(1:wakeInd,:)');

    dPhi_d_2 = 0*dPhi_d;

    if grd == 1
        [~,dPhi_d_2,~] = Phi(ones(1,length(muW(1:wakeInd))),zeros(1,length(muW(1:wakeInd))),Xc',Zc',xw(1:wakeInd+1),-zw(1:wakeInd+1),[vtw(1:wakeInd,1) -vtw(1:wakeInd,2)]',[vnw(1:wakeInd,1) -vnw(1:wakeInd,2)]');
    end

    Cw(:,1:wakeInd) = dPhi_d' + dPhi_d_2';
end

% Influence of the lumped wake panel    
 if (i_t > 1)   
    [~,dPhi_d,~] = Phi(1,0,Xc',Zc',[xw(Nlump*Nstep+1) xl(1)],[zw(Nlump*Nstep+1) zl(1)],vtlw',vnlw');

    dPhi_d_2 = 0*dPhi_d;

    if grd == 1
        [~,dPhi_d_2,~] = Phi(1,0,Xc',Zc',[xw(Nlump*Nstep+1) xl(1)],[-zw(Nlump*Nstep+1) -zl(1)],[vtlw(1,1) -vtlw(1,2)]',[vnlw(1,1) -vnlw(1,2)]');
    end

    Clw(:,1) = dPhi_d' + dPhi_d_2';
end


%% Constructing full matrix   
A = Cb + Cte;
%S = sparse(A);
%Apar = distributed(A);

% Defining Right Hand Side (RHS)
if i_t > 1
    RHS = -B*sigma - wakeinf*Cw*muW' - wakeinf*Clw*muLump(1);
else
    RHS = -B*sigma;
end
%RHSpar = distributed(RHS);

% Solving for doublet strengths.
if outerCorr <= 1
    mu(:,2:end) = mu(:,1:end-1);
end
mu(:,1) = A\RHS;

%[mu(:,1), ~] = bicgstab(A,RHS,1e-6,1000);
%[muPar(:,1), ~] = bicgstab(Apar,RHSpar,1e-6,1000);
%mu(:,1) = gather(muPar(:,1));

% Solving for TE panel strength.
if outerCorr <= 1
    muTE(1,2) = muTE(1,1);
end
muTE(1,1) = mu(Npanels,1) - mu(1,1);

% Calculating wake circulation.
if outerCorr <= 1
    GammaW(2:end) = GammaW(1:end-1);
end
if i_t == 1
    GammaW(1,1) = -muTE(1,1);
else
    GammaW(1,1) = -(muTE(1,1) - muTE(1,2));
end


%% Calculating on-body velocities and pressures. 
% Pressures are calculated from the unsteady Bernoulli's equation.  
% The flux through the body should be zero.  

% This is the order of finite difference scheme.
Order = 2; 

% Local velocity over the body due to the perturbation potential. 
Qp(1) = SecondOrderForwDiff(mu(:,1),dL,[1, 2, 3]);
Qp(2) = SecondOrderCenForwDiff(mu(:,1),dL,[1, 2, 3, 4]);

%     Qp(3:end-2,i_t) = (mu(2:end-3,i_t)- mu(4:end-1,i_t))./(1/2*dL(2:end-3,i_t) + dL(3:end-2,i_t) + 1/2*dL(4:end-1,i_t));
Pan = [1:Npanels-4; 2:Npanels-3; 3:Npanels-2; 4:Npanels-1; 5:Npanels];
for j = 1:Npanels-4  
    Qp(j+2) = CentDiff(mu(:,1),dL,Pan(:,j),Order);
end

Qp(Npanels-1) = SecondOrderCenBackDiff(mu(:,1),dL,[(Npanels - 3), (Npanels - 2), (Npanels - 1), Npanels]);
Qp(Npanels) = SecondOrderBackDiff(mu(:,1),dL,[(Npanels - 2), (Npanels - 1), Npanels]);

% Calculating the total velocity (free-stream plus the perturbation)
Qt = Qp + sum(-Vt.*vt,2);

% Calculating the pressure coefficient.
Qinf = sqrt(Q0(:,i_t)'*Q0(:,i_t));
Cp_s = 1 - Qt.^2/Qinf^2;

%%    
if i_t == 1
    Cp_us = 2/Qinf^2*mu(:,1)/delT + (sum(Vc.^2,2)/Qinf^2) - 2*Vc(:,1)/Qinf;
elseif i_t == 2
    Cp_us = 2/Qinf^2*(mu(:,1) - mu(:,2))/delT + (sum(Vc.^2,2)/Qinf^2) - 2*Vc(:,1)/Qinf;
else
    Cp_us = 2/Qinf^2*(3*mu(:,1) - 4*mu(:,2) + mu(:,3))/2/delT + (sum(Vc.^2,2)/Qinf^2) - 2*Vc(:,1)/Qinf;
end
Cp = Cp_s + Cp_us;

% Calculating the pressure coefficient.
Qp_tot = repmat(Qp,1,2).*vt + repmat(sigma,1,2).*vn;
P_s = -rho*sum(Qp_tot.^2,2)/2;

%%    
if i_t == 1
    P_us = rho*mu(:,1)/delT - rho*(sum(-Vt.*Qp_tot,2));
elseif i_t == 2
    P_us = rho*(mu(:,1) - mu(:,2))/delT - rho*(sum(-Vt.*Qp_tot,2));
else
    P_us = rho*(3*mu(:,1) - 4*mu(:,2) + mu(:,3))/2/delT - rho*(sum(-Vt.*Qp_tot,2));
end
P_b = P_s + P_us;

%% Boundary layer solver.  
if ViscDrag
    [dFi_pct,dL_pct,factor,xbc_sep,xtc_sep] = SkinFrictionSolver(xc,Qt,nu,rho,c,tmax,Cp,Qinf,0);  
    Indsep_b = 1;
    Indsep_t = 1;
    Stagpt = 75;
    dFi = dFi_pct.*kron(dL*b,ones(factor,1));
    dFshear = real(kron(eye(Npanels),ones(1,factor))*dFi);   
end

% Eliminating NaNs if present from drag calculation
dF_Ind = [];
dF_Ind = find(isnan(dFshear));
dFshear(dF_Ind) = zeros(1,length(dF_Ind));

%% Calculating the forces
% Calculating lift on the body and pitching moment about the leading
% edge.
delFp = -kron((P_b.*dL*b),[1 1]).*vn;
delFs = ViscDrag*kron(dFshear,[1 1]).*vt;
delF = delFp + delFs;
delP = sum(-delF.*Vc,2);  

F = sum(delF,1);
Fvisc = sum(delFs,1);
D_visc = Fvisc(1);
%fprintf('D_visc = %e\n',D_visc)

if BodDrag == 1
    if i_t == 1
        Fx = F(1) + 1/2*Cd_bod*rho*A_bod*Q0(1,i_t)^2;
    else
        Fx = F(1) + 1/2*Cd_bod*rho*A_bod*Q0(1,i_t-1)^2;
    end
    Fz = F(2);
else
    Fx = F(1);
    Fz = F(2);
end

Pow = sum(delP,1);

L = Fz*cos(alpha) - Fx*sin(alpha);
T = -(Fz*sin(alpha) + Fx*cos(alpha));

% Calculating non-dimensional coefficients.
Gamma = mu(Npanels,1) - mu(1,1);

Cf = norm(F)/(1/2*rho*Qinf^2*c*b);
Cl = L/(1/2*rho*Qinf^2*c*b);
Ct = T/(1/2*rho*Qinf^2*c*b);
Cpow = Pow/(1/2*rho*Qinf^3*c*b);

%% Decomposing steady and unsteady force components
delFp_s = -kron((Cp_s*1/2*rho*Qinf^2.*dL*b),[1 1]).*vn;
delFp_us = -kron((Cp_us*1/2*rho*Qinf^2.*dL*b),[1 1]).*vn;
delF_s = delFp_s;
delF_us = delFp_us;
delP_s = sum(-delF_s.*Vc,2);
delP_us = sum(-delF_us.*Vc,2);

F_s = sum(delF_s,1);
Pow_s = sum(delP_s,1);
Fx_s = F_s(1);
Fz_s = F_s(2);

F_us = sum(delF_us,1);
Pow_us = sum(delP_us,1);
Fx_us = F_us(1);
Fz_us = F_us(2);

L_s = Fz_s*cos(alpha) - Fx_s*sin(alpha);
T_s = -(Fz_s*sin(alpha) + Fx_s*cos(alpha));

Cl_s = L_s/(1/2*rho*Qinf^2*c*b);
Ct_s = T_s/(1/2*rho*Qinf^2*c*b);
Cpow_s = Pow_s/(1/2*rho*Qinf^3*c*b);

L_us = Fz_us*cos(alpha) - Fx_us*sin(alpha);
T_us = -(Fz_us*sin(alpha) + Fx_us*cos(alpha));

Cl_us = L_us/(1/2*rho*Qinf^2*c*b);
Ct_us = T_us/(1/2*rho*Qinf^2*c*b);
Cpow_us = Pow_us/(1/2*rho*Qinf^3*c*b);
