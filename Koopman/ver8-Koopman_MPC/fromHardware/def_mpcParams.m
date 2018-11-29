% def_mpcParams

%% define MPC parameters
mpcParams = struct;
mpcParams.name = 'larm_8v_horizon100_ubounded';

mpcParams.Np      = 100;                 % prediction horizon
mpcParams.nc      = 6;                % number of input constraints

% input constraints
mpcParams.umin = 0;     % min value that can be taken by any input
mpcParams.umax = 8;    % max value that can be taken by any input

%% save parameter file
save(['mpc' , filesep ,'paramFiles' , filesep , mpcParams.name , '.mat'] , 'mpcParams');