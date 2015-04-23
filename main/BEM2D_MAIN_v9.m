% Reset the system for the new run
close all
clear all
clc
pause on

% Add dependant file paths
addpath('../src/'    )
addpath('../src/geom')
addpath('../src/kine')
addpath('../src/util')
addpath('../src/fe'  )
    
% Declare program version
version = '1.1.1';
    
% Display the program title
progTitle(version);
    
% Load the input text file and interpret it
inFile = './input.dat';
loadInputData;

% Load transport properties
transportProperties;

% Load geometry, kinematic, virtual body, and panel parameters
parameters;

% Calculated simulation coditions
calculatedInput;

% Open a file for saving
saveDataFile;

% Initialize matricies
pitchingAngle = zeros(1,Ncyc*Nstep+1);
heavePos = zeros(1,Ncyc*Nstep+1);
matrixInitialize;

% Load geometry
geom;

% Initialize flight trajectory
flightTrajectory;

% Initialize the solution
initializeSolution;

for t = delT:delT:endTime
    if runTimeModifiable == 1
        loadInputData;
    end
    readFsiControls
    createInterfaceFields
    i_t = i_t + 1;
    outerCorr = 0;
    fsiResidualNorm = 1;
    if i_t > (NcycStart * Nstep + 1)
        while (fsiResidualNorm > outerCorrTolerance) && (outerCorr < nOuterCorr)
            if outerCorr == 0
                fsiResidualNorm = 0;
            end
            outerCorr = outerCorr + 1;
            % Set the interface disdplacement
            setInterfaceDisplacement

            % Move the fluid pannels
            % Solve the fluid domain
            PanelMethod2D_v9_rev_wcs211

            % Set the interface force
            setInterfaceForceFSI

            % Solve the solid domain
            solidSolveFSI_v2

            % Calculate the FSI residual
            calcFsiResidual

            if (fsiResidualNorm <= outerCorrTolerance) || (outerCorr >= nOuterCorr)
                printConvergenceInfo
                calcWakeRollup
                calcFreeSwimming
                updateStoredGeometry
                newPlotTimeStepFig
		saveDataFiles
                break
            end
        end
    else
        rigidBodySimulation
    end
end
