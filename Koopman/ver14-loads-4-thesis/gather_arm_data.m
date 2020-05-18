% gather_arm_data
%
% Perform simulations of an Arm model under various rampNhold inputs and 
% loading conditions.
%
% Generates the following objects
%   train - cell array containing data for each training trial
%   val - cell array containing data for each validation trial
%   data4sysid - struct containing all the data that can be used by Ksysid
%                class

%% Construct Arm class

Arm_setup;  % edit arm parameters in Arm_setup.m

%% Set experimental parameters

% number of trials
num_train = 10;
num_val = 3;

% size of len arrays determines the number of each type of trial
train_len = ones(1,num_train) * 60;     % length of each training trial (s)
val_len = ones(1,num_val) * 20;       % length of each valiation trial (s)

% ramp length bounds
ramps_lb = 0.25;
ramps_ub = 3;

% load condition bounds
% w_lb = [ 0 , -pi/2 ];
% w_ub = [ 1 , pi/2 ];
w_lb = [0 0];   % unloaded
w_ub = [0 0];   % unloaded

% input bounds
u_lb = -0.9*pi * ones( 1 , Arm.params.Nmods );
u_ub = 0.9*pi * ones( 1 , Arm.params.Nmods );

%% conduct training trials

train = cell( size(train_len) );
for i = 1 : length( train_len )
    % get random ramp signal for the load condition
    w_ramp_len = ( ramps_ub - ramps_lb ) * rand + ramps_lb;
    w_in = Arm.get_rampNhold( train_len(i) , w_ramp_len , w_lb , w_ub );
    
    % get random ramp signal for the input
    u_ramp_len = ( ramps_ub - ramps_lb ) * rand + ramps_lb;
    [ u_in , t_in ] = Arm.get_rampNhold( train_len(i) , u_ramp_len , u_lb , u_ub );
    
    % simulate the system
    train{i} = Arm.simulate( t_in , u_in , w_in );
end

%% conduct validation trials

val = cell( size(val_len) );
for i = 1 : length( val_len )
    % get random ramp signal for the load condition
    w_ramp_len = ( ramps_ub - ramps_lb ) * rand + ramps_lb;
    w_in = Arm.get_rampNhold( val_len(i) , w_ramp_len , w_lb , w_ub );
    
    % get random ramp signal for the input
    u_ramp_len = ( ramps_ub - ramps_lb ) * rand + ramps_lb;
    [ u_in , t_in ] = Arm.get_rampNhold( val_len(i) , u_ramp_len , u_lb , u_ub );
    
    % simulate the system
    val{i} = Arm.simulate( t_in , u_in , w_in );
end

%% put data together into a single struct

file_name = 'exp3';
data4sysid = Data.get_data4sysid( train , val , true , file_name );





