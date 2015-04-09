function [xt,zt,xb,zb] = ThinPlate(c,npts,D)
    % Stepping through each spanwise position to calculate the positions of the 
    % fin neutral plane at the given time step.
    xb = linspace(c,0,(npts+1)/2)';
    xt = linspace(0,c,(npts+1)/2)';
    zb = -0.5*D*ones((npts+1)/2,1);
    zt =  0.5*D*ones((npts+1)/2,1);

    % Circle origins ( z is assumed to be zero)
    oF = 0.5*D;
    oB = c - 0.5*D;

    % Calculate new theta positions for points on the rounded ends
    count1 = length(xb(xb>=oB));
    count2 = length(xb(xb<=oF));
    count3 = length(xt(xt<=oF));
    count4 = length(xt(xt>=oB));

    thetab = linspace(0,pi,count1+count2)';
    thetat = linspace(pi,2*pi,count3+count4)';

    % Calculate transform leading and trailing edge points
    x1 = oB + 0.5*D*cos(thetab(1:count1));
    z1 = 0  - 0.5*D*sin(thetab(1:count1));
    x2 = oF + 0.5*D*cos(thetab(end:-1:end-count1+1));
    z2 = 0  - 0.5*D*sin(thetab(end:-1:end-count1+1));
    x3 = oF + 0.5*D*cos(thetat(1:count3));
    z3 = 0  - 0.5*D*sin(thetat(1:count3));
    x4 = oB + 0.5*D*cos(thetat(end:-1:end-count3+1));
    z4 = 0  - 0.5*D*sin(thetat(end:-1:end-count3+1));

    % Replace x and z transformed points
    xb(1:count1) = x1;
    xb(end:-1:end-count2+1) = x2;
    xt(1:count3) = x3;
    xt(end:-1:end-count4+1) = x4;

    zb(1:count1) = z1;
    zb(end:-1:end-count2+1) = z2;
    zt(1:count3) = z3;
    zt(end:-1:end-count4+1) = z4;
end

