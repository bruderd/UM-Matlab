function snapshotPairs = gen_snapshotPairs( numTrials )
%gen_snapshoPairs: Runs a bunch of simulations with randomized intial
%                  conditions and compiles a set of snapshot pairs.
%   Detailed explanation goes here

%% Define system parameters (USER EDIT SECTION)
params = struct;

% initial conditions
params.phi1                = 0;
params.dtphi1              = 0;
params.phi2                = 0;
params.dtphi2              = 0;

% physical parameters
params.g                   = 9.81; 
params.m1                  = 1; 
params.m2                  = 1; 
params.l1                  = 1; 
params.l2                  = 1;

% "measurement" parameters
params.Ts                  = 1/30;
params.mean                = 0;     % mean of noise 
params.sigma               = 0;     % standard dev of noise
params.duration            = 5;   % in seconds

%% simulate a bunch of trials with randomized initial conditions
x = []; y = [];
for i = 1:numTrials

    % randomize initial conditions
    params.phi1     = pi*rand - pi/2;
    params.phi2     = (2/3)*pi*rand - pi/3;
    params.dtphi1   = 0;
    params.dtphi2   = 0;
    
    data = gen_data(params);
    
    xk = zeros(length(data.t), size(data.x,2) + size(data.u,2));
    yk = zeros(length(data.t), size(data.x,2) + size(data.u,2));
    for j = 1:length(data.t)-1
        xk(j,:) = [ data.x(j,:), data.u(j,:) ];
        yk(j,:) = [ data.x(j+1,:), data.u(j,:) ];
    end
    
    % append snapshot pairs from this trial onto set of all pairs
    x = [x; xk];
    y = [y; yk];
    
end

% define output struct
snapshotPairs = struct;
snapshotPairs.x = x;
snapshotPairs.y = y;

end

