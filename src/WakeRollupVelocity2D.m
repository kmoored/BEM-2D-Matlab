function [u,w] = WakeRollupVelocity2D(X,Z,mu,muTE,muW,muLEt,muLEb,sigma,xp,zp,vt,vn,xTE,zTE,vtTE,vnTE,xw,zw,vtw,vnw,xpt_LES,zpt_LES,vtLEt,vnLEt,xpb_LES,zpb_LES,vtLEb,vnLEb,i_t,Ncyc,Nstep,LES,epSC,epBod,SC)


% Calculating flow field

% Body contribution
[u_s,w_s,u_bd,w_bd] = DubSorV(mu,sigma,X,Z,xp,zp,vt,vn,epBod,SC);



% TE contribution
[~,~,u_TEd,w_TEd] = DubSorV(muTE,0*muTE,X,Z,xTE,zTE,vtTE,vnTE,epSC,SC);



% Wake contribution
[~,~,u_wd,w_wd] = DubSorV(muW,0*muW,X,Z,xw,zw,vtw,vnw,epSC,SC);



if LES == 1 && i_t > 1
    vec = i_t:-1:2;
    u_LEtd = zeros(1,length(X),length(vec));
    u_LEbd = zeros(1,length(X),length(vec));
    w_LEtd = zeros(1,length(X),length(vec));
    w_LEbd = zeros(1,length(X),length(vec));
    
    for j = vec
        xptLE = xpt_LES(:,Ncyc*Nstep + 2 - j)';
        zptLE = zpt_LES(:,Ncyc*Nstep + 2 - j)';
        vttLE = vtLEt(Ncyc*Nstep + 2 - j,:)';
        vntLE = vnLEt(Ncyc*Nstep + 2 - j,:)';

        xpbLE = xpb_LES(:,Ncyc*Nstep + 2 - j)';
        zpbLE = zpb_LES(:,Ncyc*Nstep + 2 - j)';
        vtbLE = vtLEb(Ncyc*Nstep + 2 - j,:)';
        vnbLE = vnLEb(Ncyc*Nstep + 2 - j,:)';

        % Influence of the LE top sheet   
        [~,~,u_LEtd(1,:,j-1),w_LEtd(1,:,j-1)] = DubSorV(muLEt(Ncyc*Nstep + 2 - j)',0*muLEt(Ncyc*Nstep + 2 - j)',X,Z,xptLE,zptLE,vttLE,vntLE,epSC,SC);

        % Influence of the LE bot sheet   
        [~,~,u_LEbd(1,:,j-1),w_LEbd(1,:,j-1)] = DubSorV(muLEb(Ncyc*Nstep + 2 - j)',0*muLEb(Ncyc*Nstep + 2 - j)',X,Z,xpbLE,zpbLE,vtbLE,vnbLE,epSC,SC);

    end
end
 

if LES == 1 && i_t > 1
    u_LEtd = sum(u_LEtd,3);
    u_LEbd = sum(u_LEbd,3);
    w_LEtd = sum(w_LEtd,3);
    w_LEbd = sum(w_LEbd,3);
    u_d = u_bd + u_TEd + u_wd + u_LEtd + u_LEbd;
    w_d = w_bd + w_TEd + w_wd + w_LEtd + w_LEbd;
else
    u_d = u_bd + u_TEd + u_wd;
    w_d = w_bd + w_TEd + w_wd;
end

u = u_d + u_s;
w = w_d + w_s;


