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

%% Former input 
getData = 'exp';       % (exp, file, or sim)
basis = 'poly';      % (fourier or poly)

%% Define system parameters (USER EDIT SECTION)
params = struct;
progress = waitbar(0,'Initializing parameters...');

params.basis = fourier;

% Koopman Sysid parameters
params.n = 3;   % dimension of state space (including state derivatives)
params.p = 1;   % dimension of input
params.naug = params.n + params.p; % dimension of augmented state (DNE)

% select maximum degrees for monomial bases (NOTE: m1 = 1)
params.maxDegree = 1;   % maximum degree of vector field monomial basis
params.m1 = 1;  % maximum degree of observables to be mapped through Lkj (DNE)

% define lifting function and basis
if strcmp(basis, 'fourier')
    params = def_fourierLift(params);  % creates fourier lifting function, fourierLift;
elseif strcmp(basis, 'poly')
    params = def_polyLift(params);  % creates polynomial lifting function, polyLift
end

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
params.ploton              = true;  % boolean to turn error plot on or off

%% Generate or load data from file
waitbar(.33,progress,'Generating data...');

addpath('generateData');

if strcmp(getData, 'sim')
    data = gen_data_fromSim( params );
elseif strcmp(getData, 'exp')
%     data = gen_data_fromExp_xonly( params );
    data = gen_data_fromExp( params );
elseif strcmp(getData, 'file')
    % Prompt user to identify data file
    [data_file,data_path] = uigetfile;
    matcontents = load([data_path, data_file]); % must be a .mat file
    data = matcontents.data;
end

rmpath('generateData')

%% Use Koopman operator to perform sysid
waitbar(.5,progress,'Performing Koopman based system identification...');

% take random subset of snapshot pairs
some_snapshotPairs = get_randsnapshotPairs(params.K, data.snapshotPairs);

koopman = koopmanSysid_CG(some_snapshotPairs, params);
% if strcmp(basis, 'fourier')
%     koopman = koopmanSysid_fourier(some_snapshotPairs, params);  % creates fourier lifting function, fourierLift;
% elseif strcmp(basis, 'poly')
%     koopman = koopmanSysid(some_snapshotPairs, params);  % creates polynomial lifting function, polyLift
% end


%% error
waitbar(0.75,progress,'Comparing to validation data set...');

% [error, xkoop] = koopmanSimulation( data.validation, params, koopman ); % only uses koopman transpose, no ODE
if strcmp(basis, 'fourier')
    [error, koopsim] = koopmanValidation_fourier( data, params, koopman );
elseif strcmp(basis, 'poly')
    [error, koopsim] = koopmanValidation( data, params, koopman );
end


%% compare koopman results to those from sysid toolbox
waitbar(0.85,progress,'Preparing data for Matlab SysId toolbox...');

% convert data to a format matlabs sysid toolbox can use
[zsysid_merged, zval_merged, zsysid, zval] = prep_iddata(data);

% save in struct for output
data4sysid = struct;
data4sysid.sysid_merged = zsysid_merged;
data4sysid.val_merged = zval_merged;
data4sysid.val = zval;
data4sysid.sysid = zsysid;
for k = 1:params.numVals
    valID = ['val', num2str(k)];
    zID = ['z', num2str(k)];
    data4sysid.valkoop.(zID) = iddata(koopsim.(valID).x, data.(valID).u, data.valparams.Ts, 'Name', 'Koopman');
end

% show comparison of Koopman system verses ground truth
if params.ploton
    for k = 1: params.numVals
        zID = ['z', num2str(k)];
        figure
        compare(data4sysid.val.(zID), data4sysid.valkoop.(zID));
    end
end

waitbar(1,progress,'Done.');
close(progress);