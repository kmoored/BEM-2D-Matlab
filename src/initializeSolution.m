t = 0;
outerCorr = 0;
i_t = 1;
createInterfaceFields
fprintf('\nInitializing the flow solution for FLOW TIME = %f\n\n',t)
PanelMethod2D_v9_rev_wcs211
newPlotTimeStepFig

nodes(:,1) = (nodes_new(:,1) - nodes_new(1,1))*cos(atan(-2*pi*ramped(i_t)*h_c*c*f*cos(2*pi*f*t)/Qinf));
nodes(:,2) = ramped(i_t)*h_c*c*sin(2*pi*f*t) + (nodes_new(:,1) - nodes_new(1,1))*sin(atan(-2*pi*ramped(i_t)*h_c*c*f*cos(2*pi*f*t)/Qinf));
