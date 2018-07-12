function [ x, y ] = get_snapshotPairs(params)
%get_snapshotPairs_SMD: Generate a bunch of snapshot pairs for the
%autonomous spring-mass-damper system (SMD)
%   Detailed explanation goes here

[Ts, numTrials, x0ub, x0lb, tf] = deal(params.Ts, params.numTrials, params.x0ub, params.x0lb, params.tf);

% initialize output arrays
x = []; y = [];
tspan = [0, tf];
for i = 1:numTrials
    x0 = (x0ub - x0lb)*rand + x0lb; % random initial state
    
    [t, v] = ode45(@(t,x) vf_real(x), tspan, x0);   % simulate plant response
    
    tq = 0:Ts:tf;
    vq = interp1(t,v,tq);   % interpolate results to get samples at sampling interval Ts
    
    for j = 1:length(vq)-1
        xk(j,:) = vq(j,:);
        yk(j,:) = vq(j+1,:);
    end
    
    % append snapshot pairs from this trial onto set of all pairs
    x = [x; xk];
    y = [y; yk];
end

end