function [ r ] = ramp( t, slope, startTime )
%RAMP: This function can generate a ramp signal based on the following
%inputs:
%   t         - vector of time samples
%   slope     - slope of the ramp signal
%   startTime - location where the ramp turns on

%   Get the number of samples in output signal
    N = numel(t);
    
%   Initialize the ramp signal
    r = zeros(N,1);
    
%   Find the index where the ramp turns on
    if median(diff(t)) > 0
        startInd = min(find(t>=startTime));
        popInd   = startInd:N;
    elseif median(diff(t)) < 0
        % Time-reversed ramp
        startTime = -startTime;
        startInd  = max(find(t>=startTime));
        popInd    = 1:startInd;
        slope     = -slope;
    end
    
%   For indices greater than the start time, compute proper signal value
%   based on slope
    r(popInd) = slope.*(t(popInd) + startTime) - 2*startTime*slope;
end

