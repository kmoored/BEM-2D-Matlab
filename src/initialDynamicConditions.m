function [ U_n, Udot_n ] = initialDynamicConditions( Nelements, theta, nodes, nodes_0,...
    nodalDelxp, nodalDelzp, nodeXOld, nodeZOld, nodeTOld, delT, order )
%UNTITLED Summary of this function goes here
%   Initialize initial condition matricies
    U_n           = zeros(3*(Nelements+1),1);
    Udot_n        = zeros(3*(Nelements+1),1);
    
%   Transform nodes to structural beam domain positions
    nodes(:,1) = nodes(:,1) - nodalDelxp;
    nodes(:,2) = nodes(:,2) - nodalDelzp; 
    [ nodes(:,1), nodes(:,2) ] = rotatePts(nodes(:,1), nodes(:,2), -theta);

%   Calculate the initial x-displacement
    U_n(1:3:end-2,1) = nodes(:,1) - nodes_0(:,1);
    
%   Calculate the initial z-displacement
    U_n(2:3:end-1,1) = nodes(:,2) - nodes_0(:,2);
    
%   Calculate the initial deflection angle
    v(:,1) = nodes(2:end,1)-nodes(1:end-1,1);
    v(:,2) = nodes(2:end,2)-nodes(1:end-1,2);
    v(:,3) = sqrt(v(:,1).^2 + v(:,2).^2);
    
    w(:,1) = nodes_0(2:end,1)-nodes_0(1:end-1,1);
    w(:,2) = nodes_0(2:end,2)-nodes_0(1:end-1,2);
    w(:,3) = sqrt(w(:,1).^2 + w(:,2).^2);
    
    U_n(6:3:end,1) = acos((v(:,1).*w(:,1)+v(:,2).*w(:,2))./(v(:,3).*w(:,3)));
    
%   Calculate the initial x-velocity
    [ Udot_n(1:3:end-2,1) ] = backwardsDifference( delT, nodes(:,1), nodeXOld(:,1), ...
        nodeXOld(:,2), nodeXOld(:,3), nodeXOld(:,4), nodeXOld(:,5), ...
        nodeXOld(:,6), order );
    
%   Calculate the initial z-velocity
    [ Udot_n(2:3:end-1,1) ] = backwardsDifference( delT, nodes(:,2), nodeZOld(:,1), ...
        nodeZOld(:,2), nodeZOld(:,3), nodeZOld(:,4), nodeZOld(:,5), ...
        nodeZOld(:,6), order );
    
%   Calculate the initial theta-velocity
    [ Udot_n(3:3:end,1)   ] = backwardsDifference( delT, U_n(3:3:end,1), nodeTOld(:,1), ...
        nodeTOld(:,2), nodeTOld(:,3), nodeTOld(:,4), nodeTOld(:,5), ...
        nodeTOld(:,6), order );
end

