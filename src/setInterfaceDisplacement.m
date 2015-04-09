%fprintf('\nTime = %f, iteration: %i\n',t,outerCorr)
if outerCorr <= 1
    fprintf('===============================================================================\n')
    fprintf(' TIME-STEP NUMBER = %i, FLOW TIME = %f\n',i_t-1,t)
    fprintf('-------------------------------------------------------------------------------\n')
    fprintf('|  Iter.  |  Relax.    |  Max Displ |  Max Res  | RMS Res Norm | Max Res Norm |\n')
    fprintf('+---------+------------+------------+-----------+--------------+--------------+\n')
end

if (outerCorr < 3 || strcmp(couplingScheme,'FixedRelaxation'))

    %fprintf('Current fsi under-relaxation factor: %f\n',fsiRelaxationFactor)
    
    fluidNodeDisplOld = fluidNodeDispl;
    nodeDisplOld = nodeDispl;
    
    fluidNodeDispl = fluidNodeDispl + fsiRelaxationFactor*fsiResidual;
    nodeDispl = nodeDispl + fsiRelaxationFactor*nodeResidual;

else
    if strcmp(couplingScheme,'Aitken')
        fsiRelaxationFactor = fsiRelaxationFactor*(fsiResidualOld' * (fsiResidualOld - fsiResidual)) / (norm(fsiResidualOld - fsiResidual))^2;
        fsiRelaxationFactor = norm(fsiRelaxationFactor);
        if fsiRelaxationFactor > 1
            fsiRelaxationFactor = 1;
        end
        %fprintf('Current fsi under-relaxation factor (Aitken): %f\n',fsiRelaxationFactor)

        fluidNodeDisplOld = fluidNodeDispl;
        nodeDisplOld = nodeDispl;
        
        fluidNodeDispl = fluidNodeDispl + fsiRelaxationFactor*fsiResidual;
        nodeDispl = nodeDispl + fsiRelaxationFactor*nodeResidual;
    else
        fprintf('ERROR! Invalid coupling scheme ''%s''\nValid coupling schemes are:\n''FixedRelaxation''\n''Aitken''\n',couplingScheme)
    end
end




