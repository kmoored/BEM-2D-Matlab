function [omega_y,Lambda2,Xstar,Zstar] = Lambda2Crit2D(Xf,Zf,u_p,w_p)

[Nz,Nx] = size(Xf); 

Xstar = (Xf(2:end-1,3:end) + Xf(2:end-1,1:end-2))/2;
Zstar = (Zf(3:end,2:end-1) + Zf(1:end-2,2:end-1))/2;

du_dx = (u_p(2:end-1,3:end) - u_p(2:end-1,1:end-2))./(Xf(2:end-1,3:end) - Xf(2:end-1,1:end-2));
dw_dx = (w_p(2:end-1,3:end) - w_p(2:end-1,1:end-2))./(Xf(2:end-1,3:end) - Xf(2:end-1,1:end-2));
du_dz = (u_p(3:end,2:end-1) - u_p(1:end-2,2:end-1))./(Zf(3:end,2:end-1) - Zf(1:end-2,2:end-1));
dw_dz = (w_p(3:end,2:end-1) - w_p(1:end-2,2:end-1))./(Zf(3:end,2:end-1) - Zf(1:end-2,2:end-1));
omega_y = du_dz - dw_dx;

Du_Dx = reshape(du_dx,1,1,(Nx - 2)*(Nz - 2)); 
Du_Dz = reshape(du_dz,1,1,(Nx - 2)*(Nz - 2)); 
Dw_Dx = reshape(dw_dx,1,1,(Nx - 2)*(Nz - 2)); 
Dw_Dz = reshape(dw_dz,1,1,(Nx - 2)*(Nz - 2)); 

LambdaMat = [Du_Dx.^2 + Dw_Dx.*Du_Dz, 0*Du_Dx; 0*Du_Dx,Dw_Dz.^2 + Dw_Dx.*Du_Dz];

Lambda2 = zeros((Nx - 2)*(Nz - 2),1);
for i = 1:(Nx - 2)*(Nz - 2)
    [~,vals] = eig(LambdaMat(:,:,i));   
    vals = diag(vals);
    vals = sort(vals,1,'ascend');
    Lambda2(i) = (vals(2,:) < 0);
end

Lambda2 = reshape(Lambda2,Nz - 2,Nx - 2);