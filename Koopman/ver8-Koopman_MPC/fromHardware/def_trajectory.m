% def_trajectory: Defines the refeference trajecory for the trial

ref = struct;

ref.name = 'blockM_larm_8x8_3min';

ref.T = 180;    % total time of trajectory (s)

%% define shape of refrence trajectory
addpath(['trajectory' , filesep , 'functions']);
ref.y = get_blockM([-2,0], 8 , 8);   % collection of points that defines the shape of the trajectory
rmpath(['trajectory' , filesep , 'functions']);


%% define time vector
ref.t = linspace( 0 , ref.T , size(ref.y,1) )'; % timestep must be the same as model.params.Ts


%% save reference trajectory struct
save(['trajectory' , filesep , 'files' , filesep , ref.name , '.mat'] , 'ref');