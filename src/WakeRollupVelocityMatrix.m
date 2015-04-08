function [u,w] = WakeRollupVelocityMatrix(X,Z,mu,muTE,muW,sigma,xp,zp,vt,vn,xTE,zTE,vtTE,vnTE,xw,zw,vtw,vnw,epSC,SCw)


% Calculating flow field

% Body contribution
[u_s,w_s,u_bd,w_bd] = DubSorV(mu,sigma,X,Z,xp,zp,vt,vn);


% TE contribution
[~,~,u_TEd,w_TEd] = DubSorV(muTE,0*muTE,X,Z,xTE,zTE,vtTE,vnTE);


% Wake contribution
[~,~,u_wd,w_wd] = DubSorV(muW,0*muW,X,Z,xw,zw,vtw,vnw);


u_d = u_bd + u_TEd + u_wd;
w_d = w_bd + w_TEd + w_wd;

u = u_d + u_s;
w = w_d + w_s;


