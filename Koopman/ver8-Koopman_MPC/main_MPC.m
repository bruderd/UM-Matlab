% main_MPC
%
%
%
%
%

%% load in system model
[model_file , model_path] = uigetfile('./models/');
model = load( [model_path , model_file] );

%% set MPC parameters
mpc = struct;

mpc.Ts      = model.params.Ts;      % sampling time must be the same as discrete model
mpc.tf      = ;                     % total length of MPC simulation
mpc.tspan   = 0 : mpc.Ts : mpc.tf;  % time vector for MPC simulation
mpc.Np      = ;                 % prediction horizon


%% define reference trajectory


%% define cost function


%% run MPC simulation