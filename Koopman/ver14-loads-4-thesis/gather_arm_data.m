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
num_val = 5;

% size of len arrays determines the number of each type of trial
train_len = ones(1,num_train) * 60;     % length of each training trial (s)
val_len = ones(1,num_val) * 20;       % length of each valiation trial (s)

% input ramp length bounds
ramps_lb = 0.25;
ramps_ub = 3;

% load ramp length bounds
w_ramps_lb = 4;
w_ramps_ub = 10;

% load condition bounds
% w_lb = [ 0 , -pi/2 ];
% w_ub = [ 10 , pi/2 ];
% w_lb = [0 0];   % endeff load only
% w_ub = [10 0];   % endeff load only
% w_lb = [ 0 , -pi/2 ];   % gravity load only
% w_ub = [ 0 , pi/2 ];    % gravity load only
w_lb = [0 0];   % unloaded
w_ub = [0 0];   % unloaded 

% input bounds
u_lb = -0.9*pi * ones( 1 , Arm.params.Nmods );
u_ub = 0.9*pi * ones( 1 , Arm.params.Nmods );

%% conduct training trials

train = cell( 1 , num_train );
for i = 1 : num_train
    % get random ramp signal for the load condition
    w_ramp_len = ( w_ramps_ub - w_ramps_lb ) * rand + w_ramps_lb;
    w_in = Arm.get_rampNhold( train_len(i) , w_ramp_len , w_lb , w_ub );
    
%     get constant load signal
%     w_in = w_lb + ( i / num_train ) * ( w_ub - w_lb );
    
    % get random ramp signal for the input
    u_ramp_len = ( ramps_ub - ramps_lb ) * rand + ramps_lb;
    [ u_in , t_in ] = Arm.get_rampNhold( train_len(i) , u_ramp_len , u_lb , u_ub );
    
    % simulate the system
    train{i} = Arm.simulate( t_in , u_in , w_in );
end

%% conduct validation trials

val = cell( 1 , num_val );
for i = 1 : num_val
    % get random ramp signal for the load condition
    w_ramp_len = ( w_ramps_ub - w_ramps_lb ) * rand + w_ramps_lb;
    w_in = Arm.get_rampNhold( val_len(i) , w_ramp_len , w_lb , w_ub );
    
%     % get constant load signal
%     w_in = w_lb + ( i / num_val ) * ( w_ub - w_lb );
    
    % get random ramp signal for the input
    u_ramp_len = ( ramps_ub - ramps_lb ) * rand + ramps_lb;
    [ u_in , t_in ] = Arm.get_rampNhold( val_len(i) , u_ramp_len , u_lb , u_ub );
    
    % simulate the system
    val{i} = Arm.simulate( t_in , u_in , w_in );
end

%% put data together into a single struct

file_name = 'debug-1link-angles-noload';
data4sysid = Data.get_data4sysid( train , val , true , file_name );





