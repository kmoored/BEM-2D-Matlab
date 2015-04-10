[ U_nPlus, Udot_nPlus, UdotDot_nPlus ] = structuralSolverFSI( ...
    alpha, beta, gamma, Nelements, fixedNodes, l_0, materialDensity, YoungMod, ...
    A, I, Fload, U_n, Udot_n, delT, fracDeltaT, 'HHT');

nodeDisplacements(:,2) =  U_nPlus(2:3:end-1,1);
nodeDisplacements(:,1) =  (nodes_0(fixedNodes+1:end,2) + nodeDisplacements(:,2)).*sin(-U_nPlus(3:3:end,1) + U_nPlus(1:3:end-2,1));

tempNodes = nodes_0;
tempNodes(fixedNodes+1:end,1) = nodes_0(fixedNodes+1:end,1) + nodeDisplacements(:,1); % New x-position
tempNodes(fixedNodes+1:end,2) = nodes_0(fixedNodes+1:end,2) + nodeDisplacements(:,2); % New z-position

frameNodes = tempNodes;


%[ tempNodes(:,1), tempNodes(:,2) ] = rotatePts(tempNodes(:,1), tempNodes(:,2),alpha_max*sin(2*pi*f*t + phi));
[ tempNodes(:,1), tempNodes(:,2) ] = rotatePts(tempNodes(:,1), tempNodes(:,2),atan(-2*pi*ramped(i_t)*h_c*c*f*cos(2*pi*f*t)/Qinf));

tempNodes(:,2) = tempNodes(:,2)+ramped(i_t)*h_c*c*sin(2*pi*f*t);
% tempNodes(:,1) = (tempNodes(:,1) - tempNodes(1,1))*cos(atan(-2*pi*h_c*c*f*cos(2*pi*f*t)/Qinf));
% tempNodes(:,2) = h_c*c*sin(2*pi*f*t) + (tempNodes(:,1) - tempNodes(1,1))*sin(atan(-2*pi*h_c*c*f*cos(2*pi*f*t)/Qinf));

tempNodes(:,1) = tempNodes(:,1) + nodalDelxp;
tempNodes(:,2) = tempNodes(:,2) + nodalDelzp; 

% figure(3)
% hold on
% plot(nodes(:,1),nodes(:,2),'-o')
% plot(tempNodes(:,1),tempNodes(:,2),'-x')
% hold off
% pause(5)
% close

normal = zeros(length(elements),2);
for i = 1:length(elements)
    normal(i,1) = -1 * (tempNodes(i+1,2) - tempNodes(i,2));
    normal(i,2) =  1 * (tempNodes(i+1,1) - tempNodes(i,1));
    normal(i,:) = normal(i,:) / sqrt(normal(i,1)^2 + normal(i,2)^2);
end

topNodes(:,1) = tempNodes(1:end-1,1) + 0.5*tBeam(:,1).*normal(:,1);
topNodes(:,2) = tempNodes(1:end-1,2) + 0.5*tBeam(:,1).*normal(:,2);
topNodes(length(topNodes)+1,1:2) = tempNodes(end,1:2);
topNodes(:,3) = tempNodes(:,3);
bottomNodes(:,1) = tempNodes(1:end-1,1) - 0.5*tBeam(:,1).*normal(:,1);
bottomNodes(:,2) = tempNodes(1:end-1,2) - 0.5*tBeam(:,1).*normal(:,2);
bottomNodes(length(bottomNodes)+1,1:2) = tempNodes(end,1:2);
bottomNodes(:,3) = tempNodes(:,3);

% Interpolate the structual results back to the fluid domain
if interpMtd == 1
    bottomXp(:,1) = interp1(bottomNodes(:,3),bottomNodes(:,1),meanline_p0(1:round(0.5*length(meanline_p0))-1,1),'linear','extrap');
    bottomZp(:,1) = interp1(bottomNodes(:,3),bottomNodes(:,2),meanline_p0(1:round(0.5*length(meanline_p0))-1,1),'linear','extrap');
    topXp(:,1)    = interp1(topNodes(:,3)   ,topNodes(:,1)   ,meanline_p0(round(0.5*length(meanline_p0)):length(meanline_p0),1),'linear','extrap');
    topZp(:,1)    = interp1(topNodes(:,3)   ,topNodes(:,2)   ,meanline_p0(round(0.5*length(meanline_p0)):length(meanline_p0),1),'linear','extrap');    
else
    bottomXp(:,1) = spline(bottomNodes(:,3),bottomNodes(:,1),meanline_p0(1:round(0.5*length(meanline_p0))-1,1));
    bottomZp(:,1) = spline(bottomNodes(:,3),bottomNodes(:,2),meanline_p0(1:round(0.5*length(meanline_p0))-1,1));
    topXp(:,1)    = spline(topNodes(:,3)   ,topNodes(:,1)   ,meanline_p0(round(0.5*length(meanline_p0)):length(meanline_p0),1));
    topZp(:,1)    = spline(topNodes(:,3)   ,topNodes(:,2)   ,meanline_p0(round(0.5*length(meanline_p0)):length(meanline_p0),1));   
end

% figure(7)
% plot(bottomXp,bottomZp,topXp,topZp)
% axis equal
% pause(5)

newxp(1:round(0.5*length(meanline_p0))-1,1)                 = bottomXp;
newxp(round(0.5*length(meanline_p0)):length(meanline_p0),1) = topXp;
newzp(1:round(0.5*length(meanline_p0))-1,1)                 = bottomZp;
newzp(round(0.5*length(meanline_p0)):length(meanline_p0),1) = topZp;

% [ newxp, newzp ] = rotatePts(newxp, newzp,atan(-2*pi*h_c*c*f*cos(2*pi*f*t)/Qinf));
%newzp = newzp + h_c*c*sin(2*pi*f*t);


for i = 1:length(xp)
    if meanline_p0(i,1) <= flexionRatio
        newxp(i,1) = xp(i,1);
        newzp(i,1) = zp(i,1);
    else
        newxp(i,1) = newxp(i,1);
        newzp(i,1) = newzp(i,1);
    end
end

% figure(7)
% plot(newxp,newzp)
% axis equal
% pause(5)

clear topNodes bottomNodes
DU = [newxp - xp, newzp - zp];

