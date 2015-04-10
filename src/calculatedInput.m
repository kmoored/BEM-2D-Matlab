Re = Qinf*c/nu;               % Reynolds number.
AoA_max = abs(atan2(pi*f*A_c*c,Qinf)...
    + alpha_max)*180/pi;      % Maximum angle of attack, degrees.

% Determine time-step
if f == 0
    delT = 1/0.001/Nstep;     % s 
else
    delT = 1/f/Nstep;         % s
end

% Free-stream velocity
Uinf = Qinf*cos(alpha);
Winf = Qinf*sin(alpha);

% Simulation end time
endTime = Nstep * Ncyc * delT;

% Determine ramping constants if applicable
if rampedStart == 1
    ramped = ramp(0:delT:endTime,1/(Nstep*delT),0) + ...
        ramp(0:delT:endTime,-1/(Nstep*delT),Nstep*delT);
else
    ramped = ones(numel(0:delT:endTime),1);
end

% Print calculations
fprintf('     alpha_max = %f\n',alpha_max);
fprintf('     Re = %e\n',Re);
fprintf('     AoA_max = %f\n',AoA_max);
fprintf('     delT = %f\n',delT);