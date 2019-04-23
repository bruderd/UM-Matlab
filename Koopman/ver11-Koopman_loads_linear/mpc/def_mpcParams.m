% def_mpcParams

%% define MPC parameters
mpcParams = struct;
mpcParams.name = 'armloadSim_horizon50_ubounded';

mpcParams.Np      = 50;                 % prediction horizon
mpcParams.nc      = 2;                % number of input constraints (usually 2 times number of inputs, because uppper and lower limit)

% input constraints
mpcParams.umin = -pi;     % min value that can be taken by any input
mpcParams.umax = pi;    % max value that can be taken by any input

% slope constraints
mpcParams.slope = ( mpcParams.umax - mpcParams.umin ) / 5;  % limits the absolute value of slope

% smoothness constraint
mpcParams.smooth = 0.04;    % larger = smoother

% input costs
mpcParams.rc = 0.1;    % running cost weight
mpcParams.tc = 100;    % terminal cost weight

% horizon for the load estimator
mpcParams.Nw = 50;

%% save parameter file
save(['paramFiles' , filesep , mpcParams.name , '.mat'] , 'mpcParams');