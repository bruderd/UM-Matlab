function [ x, y ] = get_snapshotPairs_SMD(Ts, numTrials, x0max, tf, m, b, k)
%get_snapshotPairs_SMD: Generate a bunch of snapshot pairs for the
%autonomous spring-mass-damper system (SMD)
%   Detailed explanation goes here

% initialize output arrays
x = []; y = [];

for i = 1:numTrials
    [t, v] = sim_SpringMassDamper(x0max*rand, tf, m, b, k);
    
    tq = 0:Ts:tf;
    vq = interp1(t,v,tq);
    
    for j = 1:length(vq)-1
        xk(j,:) = vq(j,:);
        yk(j,:) = vq(j+1,:);
    end
    
    % append snapshot pairs from this trial onto set of all pairs
    x = [x; xk];
    y = [y; yk];
end

end

