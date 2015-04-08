function [ub,vb] = TransformMatrixBack(t,n,up,vp)

[K,N] = size(up);           % K panels, N points.
 
i = reshape(repmat(1:2*K,2,1),4*K,1);
j = reshape(repmat(reshape((1:2*K)',2,K),2,1),4*K,1);
v = [t;n];

v = reshape(v,4*K,1);

Qt = sparse(j,i,v,2*K,2*K);

Uvec = zeros(2*K,N);
Uvec(1:2:2*K,:) = up;
Uvec(2:2:2*K,:) = vp;

Ub = Qt*Uvec;

ub = sum(Ub(1:2:2*K,:),1);
vb = sum(Ub(2:2:2*K,:),1);

