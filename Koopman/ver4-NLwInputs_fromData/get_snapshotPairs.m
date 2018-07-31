function [ x, y ] = get_snapshotPairs(params)
%get_snapshotPairs: Generate a bunch of snapshot pairs from the data by
%simulating the plant system dynamics at a bunch of initial conditions and
%inputs
%   Detailed explanation goes here

[Ts, numTrials, x0ub, x0lb, tf] = deal(params.Ts, params.numTrials, params.x0ub, params.x0lb, params.tf);

% initialize output arrays
x = []; y = [];
tspan = [0, tf];
for i = 1:numTrials
    x0 = (x0ub - x0lb)*rand + x0lb; % random initial state
    
    [t, v] = ode45(@(t,x) vf_real(x, get_u(t, params)), tspan, x0);   % simulate plant response input (defined below)
    u = get_u(t, params);
    
    tq = 0:Ts:tf;
    vq = interp1(t,v,tq);   % interpolate results to get samples at sampling interval Ts
    uq = interp1(t,u,tq)';
    
    % inject noise with standard deviation 0.01
    mean = 0;   % mean offset
    sigma = 0.01;   % standard deviation
    noise = sigma .* randn(size(vq)) + mean;
    vq = vq + noise;
    
    for j = 1:length(vq)-1
        xk(j,:) = [ vq(j,:), uq(j,:) ];
        yk(j,:) = [ vq(j+1,:), uq(j,:) ];
    end
    
    % append snapshot pairs from this trial onto set of all pairs
    x = [x; xk];
    y = [y; yk];
end

end