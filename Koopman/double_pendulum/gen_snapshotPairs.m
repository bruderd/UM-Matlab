function snapshotPairs = gen_snapshotPairs( numTrials, dataFileName )
%gen_snapshoPairs: Runs a bunch of simulations with randomized intial
%                  conditions and compiles a set of snapshot pairs.
%   Detailed explanation goes here

%% Define system parameters (USER EDIT SECTION)
params = struct;
progress = waitbar(0,'Initializing parameters...');

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
params.sigma               = 0.01;     % standard dev of noise
params.duration            = 10;   % in seconds

%% simulate a bunch of trials with randomized initial conditions
x = []; y = [];
for i = 1:numTrials
    for k = 1:10
        
        % generate parameters for input sinusoid
        params.amp = 5 * (k / 10);
        params.freq = 10 * (i / numTrials);    
    
        % randomize initial conditions
        params.phi1     = (pi/2)*rand - pi/4;
        params.phi2     = (pi/2)*rand - pi/4;
        params.dtphi1   = 0;
        params.dtphi2   = 0;
        
        % generate data from one simulation
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
    waitbar(i/numTrials, progress, 'Running Simulations...');
end

% define output struct
snapshotPairs = struct;
snapshotPairs.x = x;
snapshotPairs.y = y;

if exist('dataFileName', 'var')
    save([dataFileName, '.mat'], 'snapshotPairs');
end

waitbar(1, progress, 'Done.');

end

