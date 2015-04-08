classdef body
    properties 
        N
        c
        tmax
        xp
        zp
        xp_0
        zp_0
    end
    methods 
        function obj = body(N,c,tmax)
            obj.N = N;
            obj.c = c;
            obj.tmax = tmax;
        end
        function obj = TearDropShape(obj)
            npts = obj.N + 1;
            D = obj.tmax;
            
            % Stepping through each spanwise position to calculate the positions of the 
            % fin neutral plane at the given time step.
            xb = linspace(pi,0,(npts+1)/2)';
            xt = linspace(0,pi,(npts+1)/2)';

            % Slopes and intersects for the line segments
            m = -D/2/(obj.c - D/2);
            b = D/2 + D^2/4/(obj.c - D/2);

            % Tear drop shape equation.
            x_c = 1/2*(1 - cos(xb));
            xb = x_c*obj.c;
            xb1 = xb(xb <= D/2);
            xb2 = xb(xb > D/2);

            zb2 = -m*xb2 - b;
            zb1 = -sqrt((D/2)^2 - (xb1 - D/2).^2);
            zb = [zb2; zb1];

            % Tear drop shape equation.
            x_c = 1/2*(1 - cos(xt));
            xt = x_c*obj.c;
            xt1 = xt(xt <= D/2);
            xt2 = xt(xt > D/2);

            zt1 = sqrt((D/2)^2 - (xt1 - D/2).^2);
            zt2 = m*xt2 + b;
            zt = [zt1; zt2];
            
            xpt = [xb;xt(2:end)];
            obj.xp_0 = xpt - min(xpt);
            obj.zp_0 = [zb(1:end-1)-zb(1);0;zt(2:end)-zt(end)];

            obj.xp = obj.xp_0;      % xp will change with each timestep, but xp_0 will always remain the same.  At the first timestep, they are the same.
            obj.zp = obj.zp_0;      % zp will change with each timestep, but zp_0 will always remain the same.  At the first timestep, they are the same.

        end
    end
end