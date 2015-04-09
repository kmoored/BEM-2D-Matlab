% Specifying the airfoil shape function (Van de Vooren airfoil) and panel 
% corner points.
% [xpt,zpt,xpb,zpb,a,k,epsilon,theta1,theta2,theta,r1,r2] = VanDeVooren(c,Npanels+1,tmax/2);

% Specifying the airfoil teardrop shape from Godoy-Diana papers. 
% [xpt,zpt,xpb,zpb] = TearDropShape(c,Npanels+1,tmax);
[xpt,zpt,xpb,zpb] = ThinPlate(c,Npanels+1,tmax);
xp = [xpb;xpt(2:end)];
xp_0 = xp - min(xp);
zp_0 = [zpb(1:end-1)-zpb(1);0;zpt(2:end)-zpt(end)];

% Calculating the meanline position
meanline_p0 = xp_0 / (max(xp_0) - min(xp_0));

% Calculating the number of elements, nodes, and elements for the structual
% solver. Also declairing the piece-wise beam thickness for each node
% except the trailing edge. The rigid portion is assumed to have the
% maximum thickness 'tmax'
% TearDropShapeSolid;
ThinPlateSolid;

xp = xp_0;      % xp will change with each timestep, but xp_0 will always remain the same.  At the first timestep, they are the same.
zp = zp_0;      % zp will change with each timestep, but zp_0 will always remain the same.  At the first timestep, they are the same.
xp_new = xp_0;
zp_new = zp_0;
nodes_new(:,1:2) = nodes_0(:,1:2);
xpOld = zeros(Npanels+1,6);
zpOld = zeros(Npanels+1,6);
xpOld(:,1) = xp;
zpOld(:,1) = zp;
nodeXOld = zeros(Nelements+1,6);
nodeZOld = zeros(Nelements+1,6);
nodeXOld(:,1) = nodes_0(:,1);
nodeZOld(:,1) = nodes_0(:,2);
nodeTOld(:,1:6) = zeros(Nelements+1,6);

% Default locations of collocation points (mid-panel on the surface).
xc_0 = (xp(2:end) + xp(1:end-1))/2;
zc_0 = (zp(2:end) + zp(1:end-1))/2;
zcval = zc_0(Npanels,1);        % To be used to move the collocations points into the body for the Dirichlet formulation.

% Calculating the collocation meanline position
meanline_c0 = xc_0 / (max(xp_0) - min(xp_0));

xc = xc_0;      % xc will change with each timestep, but xc_0 will always remain the same.  At the first timestep, they are the same.
zc = zc_0;      % zc will change with each timestep, but zc_0 will always remain the same.  At the first timestep, they are the same.

% Normal and tangential vectors for the panels at collocation points.  Each row is a new
% collocation point.
vttemp = [(xp(2:end,1) - xp(1:end-1,1)) (zp(2:end,1) - zp(1:end-1,1))];
vt(:,:,1) = diag(1./sqrt(vttemp(:,1).^2 + vttemp(:,2).^2))*vttemp;
vn(:,:,1) = [-vt(:,2,1) vt(:,1,1)];