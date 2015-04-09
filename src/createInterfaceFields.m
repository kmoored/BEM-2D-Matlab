fluidNodeDispl         = zeros(Npanels+1,2);
fluidNodeDisplOld      = zeros(Npanels+1,2);
solidNodeDispl         = zeros(Npanels+1,2);
nodeDispl              = zeros(length(nodes),2);
nodeDisplOld           = zeros(length(nodes),2);
fsiResidual            = zeros(Npanels+1,2);
fsiResidualOld         = zeros(Npanels+1,2);
nodeResidual           = zeros(length(nodes),2);
nodeResidualOld        = zeros(length(nodes),2);

%if firstTime == 1
    initialFsiResidualNorm = 0;
    maxInitialFsiResidualNorm = 0;
%end

fsiResidualNorm        = 0;
maxFsiResidualNorm     = 0;
