% def_trajectory: Defines the refeference trajecory for the trial

ref = struct;

ref.name = 'setpoint_1p5pi_10sec';

ref.T = 10;    % total time of trajectory (s)

%% define shape of refrence trajectory
addpath(['trajectory' , filesep , 'functions']);
ref.y = [ 3*pi/2 , 0 ];
% ref.y = get_blockM([0,0], 8 , 8);   % collection of points that defines the shape of the trajectory
% ref.y = -[ ref.y(:,2) , -ref.y(:,1) ];    % include after the block M to turn it into sigma
% ref.y = get_circle([0,1] , 1);
% ref.y = get_reachPoint([-0.4,-0.6]);
% ref.y = get_reachLine([0,-6]);
% ref.y = get_polygon( ( [1 2 ; 2 0 ; 4 0 ; 2.5 -1.5 ; 3 -3.5 ; 1 -2.5 ; -1 -3.5 ; -0.5 -1.5 ; -2 0 ; 0 0] + [0,0] ).* (9/6) ); % STAR
% ref.y = get_polygon( [0 0; 1 1 ; 2 1.5 ; 3 1 ; 4 0 ; 4 -1 ; 3 -2 ; 2 -3 ; 1 -4 ; 0 -5 ;...
%                     -1 -4; -2 -3; -3 -2; -4 -1; -4 0; -3 1; -2 1.5; -1 1; 0 0] ); % HEART
% ref.y = get_sinusoid( [ 1 1 ] , [ 4 4 ] , 30 );
% ref.y = get_pacman( [0,0] , 3 );
rmpath(['trajectory' , filesep , 'functions']);

%% ensure trajectory starts from resting configuration of system
% ref.y = [ 0, 2 ; ref.y ];   % planar manipulator
ref.y = [ 0, 0 ; ref.y ];   % laser tracker

%% define time vector
ref.t = linspace( 0 , ref.T , size(ref.y,1) )'; % timestep must be the same as model.params.Ts


%% save reference trajectory struct
save(['trajectory' , filesep , 'files' , filesep , ref.name , '.mat'] , 'ref');