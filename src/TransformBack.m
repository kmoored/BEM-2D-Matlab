function [ub,vb,wb] = TransformBack(vc,vn,vt,up,vp,wp)

[K,~] = size(up);

Nbarpri = sparse(eye(2*K));
Nbar = sparse([Nbarpri(:,1:2:2*K) Nbarpri(:,2:2:2*K)]);

clear Nbarpri

Uvec = [up;vp;wp];
Uvec = Nbar*Uvec;
v = [t;n];
vpri = Nbar*v;

R = zeros(2*K,2*K);
for i = 1:2:2*K
    R(i:i+1,i:i+1) = vpri(i:i+1,:);
end

Ub = R'*Uvec;

ub = sum(Ub(1:2:2*K,:),1);
vb = sum(Ub(2:2:2*K,:),1);
