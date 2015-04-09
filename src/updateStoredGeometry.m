%[xp_new,zp_new] = antiKinematics_HeavePitch2D(xp,zp,alpha_max,h_c,f,t,phi,'fluid');
%[nodes_new(:,1),nodes_new(:,2)] = antiKinematics_HeavePitch2D(nodes(:,1),nodes(:,2),alpha_max,h_c,f,t,phi,'solid');
[ xp_new,zp_new ] = antiKinematics_ZeroAoA( xp,zp,h_c,c,f,t,Qinf,'fluid',ramped,i_t );
[ nodes_new(:,1),nodes_new(:,2) ] = antiKinematics_ZeroAoA( nodes(:,1),nodes(:,2),h_c,c,f,t,Qinf,'solid',ramped,i_t );

xpOld = [xp-delxp xpOld(:,1) xpOld(:,2) xpOld(:,3) xpOld(:,4) xpOld(:,5)];
zpOld = [zp-delzp zpOld(:,1) zpOld(:,2) zpOld(:,3) zpOld(:,4) zpOld(:,5)];
nodeXOld = [nodes(:,1)-nodalDelxp nodeXOld(:,1) nodeXOld(:,2) nodeXOld(:,3) nodeXOld(:,4) nodeXOld(:,5)];
nodeZOld = [nodes(:,2)-nodalDelzp nodeZOld(:,1) nodeZOld(:,2) nodeZOld(:,3) nodeZOld(:,4) nodeZOld(:,5)];

% nodeXOld = [frameNodes(:,1) nodeXOld(:,1) nodeXOld(:,2) nodeXOld(:,3) nodeXOld(:,4) nodeXOld(:,5)];
% nodeZOld = [frameNodes(:,2) nodeZOld(:,1) nodeZOld(:,2) nodeZOld(:,3) nodeZOld(:,4) nodeZOld(:,5)];
% 
% vect1(:,1) = nodes(2:end,1)-nodes(1:end-1,1);
% vect1(:,2) = nodes(2:end,2)-nodes(1:end-1,2);
% vect1(:,3) = sqrt(vect1(:,1).^2 + vect1(:,2).^2);
% vect2(:,1) = nodes_0(2:end,1)-nodes_0(1:end-1,1);
% vect2(:,2) = nodes_0(2:end,2)-nodes_0(1:end-1,2);
% vect2(:,3) = sqrt(vect2(:,1).^2 + vect2(:,2).^2);
% 
% nodeTheta = acos((vect1(:,1).*vect2(:,1)+vect1(:,2).*vect2(:,2))./(vect1(:,3).*vect2(:,3)));
% tmpNodeTheta = zeros(Nelements+1,1);
% tmpNodeTheta(2:end,1) = nodeTheta;
% nodeTOld = [tmpNodeTheta nodeTOld(:,1) nodeTOld(:,2) nodeTOld(:,3) nodeTOld(:,4) nodeTOld(:,5)];