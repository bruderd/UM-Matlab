% gather_arm_data_mesh
%
% Perform simulations of an Arm model under various rampNhold inputs and 
% loading conditions.
%
% Differs from gather_arm_data in the way training loads are selected.
% It covers a whole grid of constant loads rather than ramping them.
%
% Generates the following objects
%   train - cell array containing data for each training trial
%   val - cell array containing data for each validation trial
%   data4sysid - struct containing all the data that can be used by Ksysid
%                class

%% Construct Arm class

Arm_setup;  % edit arm parameters in Arm_setup.m

%% Set experimental parameters

% size of len arrays determines the number of each type of trial
train_len = 40;     % length of each training trial (s)
val_len = 20;       % length of each valiation trial (s)

% input ramp length bounds
ramps_lb = 0.25;
ramps_ub = 3;

% % load ramp length bounds
% w_ramps_lb = 5;
% w_ramps_ub = 10;

% load condition bounds
% w_lb = [ 0 , -pi/2 ];
% w_ub = [ 10 , pi/2 ];
% w_lb = [0 0];   % endeff load only
% w_ub = [10 0];   % endeff load only
w_lb = [ 0 , -pi/2 ];   % gravity load only
w_ub = [ 1 , pi/2 ];    % gravity load only
% w_lb = [0 0];   % unloaded
% w_ub = [0 0];   % unloaded 

% load condition meshgrid
numsteps_w_train = [ 2 , 11 ];
numsteps_w_val = [ 2 , 5 ];
[ W1_train , W2_train ] = meshgrid( linspace( w_lb(1) , w_ub(1) , numsteps_w_train(1) ) ,...
                                linspace( w_lb(2) , w_ub(2) , numsteps_w_train(2) ) );
[ W1_val , W2_val ] = meshgrid( linspace( w_lb(1) , w_ub(1) , numsteps_w_val(1) ) ,...
                                linspace( w_lb(2) , w_ub(2) , numsteps_w_val(2) ) );

% input bounds
u_lb = -0.9*pi * ones( 1 , Arm.params.Nmods );
u_ub = 0.9*pi * ones( 1 , Arm.params.Nmods );

% number of trials
num_train = numel(W1_train);
num_val = numel(W1_val);

%% conduct training trials

train = cell( 1 , num_train );
for i = 1 : num_train
    
    % get constant load signal
%     w_in = w_lb + ( (i-1) / (num_train-1) ) * ( w_ub - w_lb );
    w_in = [ W1_train(i) , W2_train(i) ];
    
    % get random ramp signal for the input
    u_ramp_len = ( ramps_ub - ramps_lb ) * rand + ramps_lb;
    [ u_in , t_in ] = Arm.get_rampNhold( train_len , u_ramp_len , u_lb , u_ub );
    
    % simulate the system
    train{i} = Arm.simulate( t_in , u_in , w_in );
end

%% conduct validation trials

val = cell( 1 , num_val );
for i = 1 : num_val
    
    % get constant load signal
%     w_in = w_lb + ( (i-1) / (num_val-1) ) * ( w_ub - w_lb );
    w_in = [ W1_val(i) , W2_val(i) ];
    
    % get random ramp signal for the input
    u_ramp_len = ( ramps_ub - ramps_lb ) * rand + ramps_lb;
    [ u_in , t_in ] = Arm.get_rampNhold( val_len , u_ramp_len , u_lb , u_ub );
    
    % simulate the system
    val{i} = Arm.simulate( t_in , u_in , w_in );
end

%% put data together into a single struct

file_name = 'thesis-3link-markers-grav-endload-01';
data4sysid = Data.get_data4sysid( train , val , true , file_name );