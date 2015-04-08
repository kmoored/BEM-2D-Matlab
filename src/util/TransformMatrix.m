% (xp, zp) -- left coordinate point for each panel.
% t -- unit vector along the panel.  Running from point 1 to point 2.
% n -- unit vector normal to the panel.
% (x,z) -- point of interest.

function [Xp] = TransformMatrix(xp,zp,t,n,x,z)

K = length(xp);         % K panels
N = length(x);          % N points
 
i = reshape(repmat(1:2*K,2,1),4*K,1);
j = reshape(repmat(reshape((1:2*K)',2,K),2,1),4*K,1);
v = [t;n];

v = reshape(v,4*K,1);

Qt = sparse(i,j,v,2*K,2*K);

Xp = Qt*(repmat([x;z],K,1) - repmat(reshape([xp;zp],2*K,1),1,N));


