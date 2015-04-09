pivotPoint = (0.5*tmax) / (max(xp_0) - min(xp_0)); 
Nelements = elementFactor * (max(xp_0) - min(xp_0)) / (0.5*tmax);
nodes(:,1) = [min(xp_0):(max(xp_0)-min(xp_0))/Nelements:max(xp_0)]';
nodes(:,2) = 0;
nodes(:,3) = nodes(:,1) / (max(nodes(:,1)) - min(nodes(:,1)));
nodes_0 = nodes;
elements = [[1:length(nodes)-1]',[2:length(nodes)]'];
elements(:,3) = 0.5*(nodes(elements(:,1),1)+nodes(elements(:,2),1));
elements(:,4) = elements(:,3) / (max(xp_0) - min(xp_0));
tBeam = zeros(Nelements,1);
ttemp = tBeam;
tBeamStruct = tBeam;
beamCounter = 0;
fixedCounter = 0;
for i = 1:Nelements
    if nodes(i,1) <= 0.5*tmax
        tBeam(i,1) = tmax;
        tBeamStruct(i,1) = tBeam(i,1);
        ttemp(i,1) = tBeam(i,1);
        beamCounter = beamCounter + 1;
        fixedCounter = fixedCounter + 1;
    elseif nodes(i,1) >= c-0.5*tmax
        tBeam(i,1) = 2*sqrt((0.5*tmax)^2 -(nodes(i,1)-(c-0.5*tmax))^2 );
        tBeamStruct(i,1) = tBeam(i,1);
        ttemp(i,1) = tBeam(i,1);
    else
        ttemp(i,1) = tmax;
        tBeam(i,1) = tmax;
        tBeamStruct(i,1) = tBeam(i,1);
        if (constThickBeam == 1 && nodes(i,3) >= tConst)
            tBeamStruct(i,1) = tBeamStruct(i-1,1);
        else
            tBeamStruct(i,1) = tBeam(i,1);
        end
        if nodes(i,3) <= flexionRatio
            fixedCounter = fixedCounter + 1;
        end
    end
end