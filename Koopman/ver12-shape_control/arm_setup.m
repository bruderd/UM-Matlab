% arm_setup.m
%
% Sets the values of parameters and derives equations of motion for a 
% hyper-redundant planar manipulator with the following properties:
%   -Every link has the same mass, inertia, length
%   -The input is groups of joint-torques
%   -Every module is identical (same number of joints/links)
%

%% Define parameters
params = struct;

params.sysName = 'arm_3-mods_1-links_01-Ts';

params.Nmods = 3;   % number of modules (actuated sections)
params.nlinks = 1;      % number of links in each module
params.Nlinks = params.Nmods * params.nlinks;   % total number of links in robot

% general system parameters (make sure to include these an any system class)
params.nx = params.Nlinks * 2;   % dimension of the full state (joint angles and joing velocities)
param.ny = 2 * (params.Nlinks ) + 2;   % dimension of measured output (mocap coordinates + end effector orientation)
params.nu = params.Nlinks;  % dimension of the input (reference angle at each joint)

% manipulator parameters
params.L = 0.3;    % total length of robot (m)
params.l = params.L / params.Nlinks;
params.k = -0.00001;    % stiffness at each joint
params.d = 1e-4;    % viscous damping at each joint
params.m = 0.0001;   % mass of each link (kg)
params.i = (1/3) * params.m * params.l^2;   % inertia of each link
params.g = 9.81;    % gravity constant (m/s^2)

% mocap parameters
params.markerPos = ( ( 0 : params.Nmods ) * params.l * params.nlinks ) / params.L;  % position of mocap markers along the arm

% input parameters
params.ku = 1e-3; % effective input stiffness

% simulation parameters
params.Ts = 1e-2;   % sampling time
params.umax = pi/2; % maximum input value (scalar for all modules, vector for different limits per module)

%% Derive the equations of motion
% EOM = arm_set_EOM(params);
% 
% % Save functions handles for functions associated with EOM
% fcns = EOM.fcns;
% 
%% Create class for this system

% arm = arm( params , fcns );
arm = arm( params );

% save this system for later use
dirname = [ 'systems' , filesep , params.sysName ];
unique_dirname = auto_rename( dirname , '(0)' );
arm.params.sysName = erase( unique_dirname , ['systems', filesep] ) ;    % modify the system name parameter    

% create directory for the system, and save the class
mkdir( unique_dirname );
mkdir( [ unique_dirname , filesep , 'simulations' ] );  % make simulation subfolder
fname = [ unique_dirname , filesep , params.sysName, '.mat' ];
save( fname , 'arm' );

