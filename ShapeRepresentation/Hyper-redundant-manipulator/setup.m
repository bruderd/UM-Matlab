% setup.m
%
% Sets the values of parameters and derives equations of motion for a 
% hyper-redundant planar manipulator with the following properties:
%   -Every link has the same mass, inertia, length
%   -The input is groups of joint-torques
%   -Every module is identical (same number of joints/links)
%

%% Define parameters
params = struct;

params.sysName = 'arm_3-mods_1-links_taylor';

params.Nmods = 3;   % number of modules (actuated sections)
params.nlinks = 1;      % number of links in each module
params.Nlinks = params.Nmods * params.nlinks;   % total number of links in robot

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
params.Ts = 1e-3;   % sampling time
params.umax = pi/2; % maximum input value (scalar for all modules, vector for different limits per module)

%% Derive the equations of motion
EOM = set_EOM(params);

% Save functions handles for functions associated with EOM
fcns = EOM.fcns;

%% Create class for this system

arm = arm( params , fcns );

% save this system for later use
fname = [ 'systems' , filesep , params.sysName , '.mat' ];
unique_fname = auto_rename( fname , '(0)' );
save( unique_fname , 'arm' );

