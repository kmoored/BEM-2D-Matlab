%% Updating the kinematics.
if outerCorr > 1
    % Superposing the structual displacements.
    nodes(:,1) = nodes(:,1) + nodeDispl(:,1) - nodeDisplOld(:,1);
    nodes(:,2) = nodes(:,2) + nodeDispl(:,2) - nodeDisplOld(:,2);
end

if outerCorr <= 1
    % Updating to the new kinematics
    %nodes(:,1) = (nodes_new(:,1) - nodes_new(1,1))*cos(alpha_max*sin(2*pi*f*t + phi));
    nodes(:,1) = (nodes_new(:,1) - nodes_new(1,1))*cos(atan(-2*pi*ramped(i_t)*h_c*c*f*cos(2*pi*f*t)/Qinf));
    %nodes(:,2) = h_c*c*sin(2*pi*f*t) + (nodes_new(:,1) - nodes_new(1,1))*sin(alpha_max*sin(2*pi*f*t + phi));
    nodes(:,2) = ramped(i_t)*h_c*c*sin(2*pi*f*t) + (nodes_new(:,1) - nodes_new(1,1))*sin(atan(-2*pi*ramped(i_t)*h_c*c*f*cos(2*pi*f*t)/Qinf));

    % Calculating shift in panel positions with the swimming velocity.
    nodalDelxp = x_b(2)*ones(length(nodes),1);
    nodalDelzp = z_b(2)*ones(length(nodes),1);

    % Superposing the kinematics and swimming translations.
    nodes(:,1) = nodes(:,1) + nodalDelxp;
    nodes(:,2) = nodes(:,2) + nodalDelzp; 
end

% if debugPlots == 1
%     figure(6)
%     %plot(nodes(:,1),nodes(:,2));
%     plot(nodes_new(:,1),nodes_new(:,2))
%     %axis([min(nodes_new(:,1)) max(nodes_new(:,2)) -0.5 0.5])
%     drawnow
% end

%% Determine load conditions from the fluid solver
% Calculate Pannel Lengths
lp = zeros(length(xc),1);
for i = 1:length(xc)
    lp(i,1) = sqrt((xp(i+1) - xp(i))^2 + (zp(i+1) - zp(i))^2);
end

% Calculate the force magnitude acting on the pannel due to pressure
% Then, calculate the x-z componets of this force



magPF = P_b .* lp * 1;
if ViscDrag == 1
    pF = [(magPF .* vn(:,1) * -1) + delFs(:,1) , (magPF .* vn(:,2) * -1) + delFs(:,2)];
else
    pF = [(magPF .* vn(:,1) * -1), (magPF .* vn(:,2) * -1)];
end

% Determine moment arm between top and bottom pannel points
% Collapse data to the meanline
colM = zeros(0.5*length(xc),1);
colPF = zeros(0.5*length(xc),2);
meanPt = zeros(0.5*length(xc),2);
for i = 1:0.5*length(pF)
    meanPt(i,1) = 0.5*(xc(i,1) + xc(end-i+1));
    meanPt(i,2) = 0.5*(zc(i,1) + zc(end-i+1));
    colPF(i,:) = pF(i,:) + pF(end-i+1,:);
    colM(i,1) = -1*pF(i,1)*(zc(i,1)-meanPt(i,2)) + pF(i,2)*(xc(i,1)-meanPt(i,1)) + -1*pF(end-i+1,1)*(zc(end-i+1)-meanPt(i,2)) + pF(end-i+1,2)*(xc(end-i+1) - meanPt(i,1));
end
colPF = flipud(colPF);
colM = flipud(colM);
relXp = (xp - min(xp));
relXc = (xc - min(xp));

% Interpolate the collapsed data onto the structual mesh
if interpMtd == 1
    nodalInput(:,1) = interp1(meanline_c0(0.5*length(meanline_c0)+1:length(meanline_c0),1), colPF(:,1), nodes(:,3),'linear','extrap');
    nodalInput(:,2) = interp1(meanline_c0(0.5*length(meanline_c0)+1:length(meanline_c0),1), colPF(:,2), nodes(:,3),'linear','extrap');
    nodalInput(:,6) = interp1(meanline_c0(0.5*length(meanline_c0)+1:length(meanline_c0),1), colM(:,1),  nodes(:,3),'linear','extrap');
else
    nodalInput(:,1) = spline(meanline_c0(0.5*length(meanline_c0)+1:length(meanline_c0),1), colPF(:,1), nodes(:,3));
    nodalInput(:,2) = spline(meanline_c0(0.5*length(meanline_c0)+1:length(meanline_c0),1), colPF(:,2), nodes(:,3));
    nodalInput(:,6) = spline(meanline_c0(0.5*length(meanline_c0)+1:length(meanline_c0),1), colM(:,1),  nodes(:,3));
end

%[ nodalInput(:,1), nodalInput(:,2) ] = rotatePts(nodalInput(:,1), nodalInput(:,2),-alpha_max*sin(2*pi*f*t + phi));
[ nodalInput(:,1), nodalInput(:,2) ] = rotatePts(nodalInput(:,1), nodalInput(:,2),atan(2*pi*ramped(i_t)*h_c*c*f*cos(2*pi*f*t)/Qinf));


% Create load matrix
Fload = zeros(3*(Nelements+1),1);
Fload(1:3:end-2,1) = nodalInput(:,1);
Fload(2:3:end-1,1) = nodalInput(:,2);
Fload(3:3:end,1)   = nodalInput(:,6);

% Create element area matrix
A = tBeamStruct(:,1);

% Create area moment of inertia matrix
I = 1 .* tBeamStruct(:,1).^3 ./ 12;

% Select fixed nodes
fixedNodes = fixedCounter;

% Initial element length
l_0 = c/Nelements;

% Inital displacements and velocities
% U_n           = zeros(3*(Nelements+1),1);
% Udot_n        = zeros(3*(Nelements+1),1);
% UdotDot_n     = zeros(3*(Nelements+1),1);
% U_nPlus       = zeros(3*(Nelements+1),1);
% Udot_nPlus    = zeros(3*(Nelements+1),1);
% UdotDot_nPlus = zeros(3*(Nelements+1),1);
% [ U_n, Udot_n ] = initialDynamicConditions( Nelements, alpha_max*sin(2*pi*f*t + phi), nodes, nodes_0,...
%     nodalDelxp, nodalDelzp, nodeXOld, nodeZOld, nodeTOld, delT, min(i_t-1,2) );

if i_t <= 2 && outerCorr <= 1
    U_n           = zeros(3*(Nelements+1),1);
    Udot_n        = zeros(3*(Nelements+1),1);
    UdotDot_n     = zeros(3*(Nelements+1),1);
    U_nPlus       = zeros(3*(Nelements+1),1);
    Udot_nPlus    = zeros(3*(Nelements+1),1);
    UdotDot_nPlus = zeros(3*(Nelements+1),1);
elseif (i_t > 1 && outerCorr <= 1)
    U_n           = zeros(3*(Nelements+1),1);
    Udot_n        = zeros(3*(Nelements+1),1);
    U_n(3*fixedNodes+1:end,1) = U_nPlus;
    Udot_n(3*fixedNodes+1:end,1) = Udot_nPlus;
end
% if i_t > 1
%     Udot_n(1:3:end-2,1) = backwardsDifference(delT,nodes(:,1)-nodalDelxp, nodeXOld(:,1), nodeXOld(:,2), nodeXOld(:,3), nodeXOld(:,4), nodeXOld(:,5),nodeXOld(:,6),min(i_t-1,2));
%     Udot_n(2:3:end-1,1) = backwardsDifference(delT,nodes(:,2)-nodalDelzp, nodeZOld(:,1), nodeZOld(:,2), nodeZOld(:,3), nodeZOld(:,4), nodeZOld(:,5),nodeZOld(:,6),min(i_t-1,2));
% end

