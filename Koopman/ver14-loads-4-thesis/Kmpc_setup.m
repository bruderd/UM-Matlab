% Kmpc_setup.m
%
% Creates an instance of the Kmpc class bases on user arguments.
% Must have a Ksysid class in the workspace for this to work

Kmpc = Kmpc( Ksysid ,...
        'horizon' , 10 ,...
        'input_bounds' , [ -Inf , Inf ],... %[ -7*pi/8 , 7*pi/8 ] ,...
        'input_slopeConst' , [1e-1],... %1e-1 ,...
        'input_smoothConst' , [],... %[1e-1] ,...
        'state_bounds' , [] ,...
        'cost_running' , 1 ,...   % 0.1
        'cost_terminal' , 100 ,...  % 100
        'cost_input' , 0 ,...    % 1e-1
        'projmtx' , Ksysid.model.C(end-1:end,:) ,...  % just end effector
        'load_obs_horizon' , 50);   % only needed for loaded models