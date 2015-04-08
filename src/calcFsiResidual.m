solidNodeDispl = DU;
%fprintf('Displacements (Max, Min, Avg) = (%e, %e, %e)\n',max(sqrt(DU(:,1).^2 + DU(:,2).^2)),min(sqrt(DU(:,1).^2 + DU(:,2).^2)),mean(sqrt(DU(:,1).^2 + DU(:,2).^2)));
fsiResidualOld = fsiResidual;
nodeResidualOld = nodeResidual;
fsiResidual = solidNodeDispl - fluidNodeDispl;
nodeResidual = (tempNodes(:,1:2)-nodes(:,1:2)) - nodeDispl;
magFsiResidual = sqrt(fsiResidual(:,1).^2 + fsiResidual(:,2).^2);
%fprintf('Max fsi residual = %e\n',max(magFsiResidual))
fsiResidualNorm = norm(fsiResidual,2);
maxFsiResidualNorm = norm(fsiResidual,'inf');
%if (outerCorr == 1 && firstTime == 1)
if (outerCorr == 1 )
    initialFsiResidualNorm = fsiResidualNorm;
    maxInitialFsiResidualNorm = maxFsiResidualNorm;
    %firstTime = 0;
end
fsiResidualNorm = fsiResidualNorm / initialFsiResidualNorm;
maxFsiResidualNorm = maxFsiResidualNorm / maxInitialFsiResidualNorm;

%fprintf('Current fsi residual norm (RMS, MAX): (%e, %e)\n\n',fsiResidualNorm,maxFsiResidualNorm)

fprintf('| %7i |   %.2E |   %.2E |  %.2E |     %.2E |     %.2E |\n',outerCorr,fsiRelaxationFactor,max(sqrt(DU(:,1).^2 + DU(:,2).^2)),max(magFsiResidual),fsiResidualNorm,maxFsiResidualNorm)
