%main_sysid
%main_sysid: A generic "main" function for learning model from data
%   Performs linear system identification of nonlinear systems using a
%   lifting technique based on the Koopman operator projected onto a finite
%   monomial basis.
%
%   INPUTS:
%       getData is a string that specifies whether data will be generated
%       by simulation, from experimental measurements, or loaded from a 
%       file. Possible values are 'sim', 'exp', and 'file'.
%
%   OUTPUTS:
%       koopman is a struct containing the value of the learned Koopman
%       operator and other values related to the koopman sysid. See the
%       koopmanSysid function for more details.


%% Define system parameters (USER EDIT SECTION)
params = struct;

params.getData = 'exp';            % ('exp, 'file', or 'sim')
params.basis = 'fourier';   % ('fourier' or 'poly')

% Koopman Sysid parameters
params.n = 3;   % dimension of state space (including state derivatives)
params.p = 1;   % dimension of input
params.naug = params.n + params.p; % dimension of augmented state (DNE)

% select maximum "degree" for basis elements (NOTE: m1 = 1)
params.maxDegree = 1;   % maximum degree of vector field monomial basis
params.m1 = 1;  % maximum degree of observables to be mapped through Lkj (DNE)

% define lifting function and basis
disp('Defining basis of observables...')
if strcmp(basis, 'fourier')
    params = def_fourierLift(params);  % creates fourier lifting function, fourierLift;
elseif strcmp(basis, 'poly')
    params = def_polyLift(params);  % creates polynomial lifting function, polyLift
end
disp('Done.')

% Another Koopman fitting parameter to penalize model complexity
params.t = (1/params.N) * params.N^2; % penalty on model complexity

% parameters for reading in data
params.numTrials        = 6;        % numer of sysid trials
params.numVals          = 1;        % number of validation trials
params.Ts               = 0.02;     % sampling period
params.K                = 5000;     % numer of snapshotPairs to take
params.numericalDerivs  = false;    % choose whether or not to take numerical derivatives of states (boolean)

params.systemName          = 'snake_5000pts_scale1_fourierBasis_allData';  % name of current system
params.filterWindow        = floor( [1/params.Ts, 1/params.Ts] );  % if taking numerical derivatives, specifies the moving mean window before and after derivatives taken.

% output parameters
params.compareon           = true;  % boolean to decide whether to compare model simulation to validation data
params.ploton              = true;  % boolean to turn error plot on or off

%% Get training data
[data, some_snapshotPairs] = get_trainingData(params);

%% Learn the approximate Koopman operator and corresponding NL system
koopman = koopmanSysid_CG(some_snapshotPairs, params);

%% Simulate the results and compare to validation trial(s)

