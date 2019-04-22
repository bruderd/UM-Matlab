%main_sysid_loaded
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

if ~exist('params' ,'var')  % recycle struct from previous run 
    params = struct;
end
params.getData = 'file';            % ('exp' or 'file')
params.basisID = 'poly';   % ('fourier' or 'poly' or 'fourier_sparser' or 'thinplate' or 'gaussian' or 'hermite')

% parameters for reading in data (these affect how shapshot pairs built from raw data).
params.numTrials        = 10;        % numer of sysid trials
params.numVals          = 10;        % number of validation trials
params.Ts               = 0.02;     % sampling period
params.K                = 100*50*10;     % numer of snapshotPairs to take
params.numericalDerivs  = false;    % choose whether or not to take numerical derivatives of states (boolean)
params.scale            = 0.9;      % scale down all state to be in range [-scale , scale]
params.nd               = 1;        % number of delays to include in the snapshot pairs

params.name             = 'arm_loaded_load0to9';    % user defined custom part of name
params.system              = [ params.name , '_' , num2str(params.K) , 'sp_' , num2str(params.numTrials) , 'sysid_' ,  num2str(params.numVals) , 'val_' , num2str(params.nd), 'del']; % name of the current system
% params.filterWindow        = floor( [1/params.Ts, 1/params.Ts] );  % if taking numerical derivatives, specifies the moving mean window before and after derivatives taken.
params.filterWindow        = floor( [1, 1] );  % no smoothing
% params.filterWindow        = floor( [6, 6] );  % if taking numerical derivatives, specifies the moving mean window before and after derivatives taken.

% Koopman Sysid parameters
params.n = 2;   % dimension of state space (including state derivatives)
params.p = 1;   % dimension of input
params.ny = 2;  % dimension of output
params.nw = 1;  % dimension of load
params.naug = params.n + params.p; % dimension of augmented state (DNE)
params.nzeta = params.n + params.nd * (params.naug);    % dimensinon of zeta (DNE)

% select maximum "degree" for basis elements
params.maxDegree = 4;   % maximum degree of vector field monomial basis

% system name imcludeing the degree 
params.systemName = [ params.system , '_' , params.basisID , num2str(params.maxDegree) , '_'];

% only do this if the Basis is not already defined. Will need to clear before running with a different basis or maxDegree
if ~isfield(params , 'Basis')   
    % define lifting function and basis
    disp('Defining basis of observables...')
    if strcmp(params.basisID, 'fourier')
        params = def_fourierLift(params);  % creates fourier lifting function, fourierLift;
    elseif strcmp(params.basisID, 'poly')
        params = def_polyLift_loaded(params);  % creates polynomial lifting function, polyLift
    elseif strcmp(params.basisID, 'fourier_sparser')
        params = def_fourierLift_sparser(params);
    elseif strcmp(params.basisID, 'thinplate')
        params = def_thinplateLift(params);
    elseif strcmp(params.basisID, 'gaussian')
        params = def_gaussianLift(params);
    elseif strcmp(params.basisID, 'hermite')
        params = def_hermiteLift(params);
    end
    disp('Done.')
end


% Koopman sysid tuning parameters
params.t        = 1 * params.N; % penalty on model complexity
params.epsilon  = 1; % model accuracy tolerance (larger value = less accurate)
params.percSat  = 0.75;  % percentage of snapshot pairs that must satisfy accuracy tolerance

% output parameters
params.validateon          = true;  % boolean to decide whether to validate the model
params.ploton              = true;  % boolean to turn validation/error/compare plot on or off
params.compareon           = true;  % boolean to decide whether to convert to iddata and compare to validation data with Matlab compare function


%% Get training data
[data, some_snapshotPairs, params] = get_trainingData(params);

%% Learn the approximate Koopman operator and corresponding NL system
[U , koopData ] = get_Koopman_loaded( some_snapshotPairs, params );
% statespace      = sysid_statespaceSys( U, some_snapshotPairs, params );
lifted          = sysid_liftedSys_loaded( U , params , koopData );

%% Simulate the results and compare to validation trial(s)
if params.validateon
    disp('Comparing to validation data set...');
%     [error, koopsim] = koopmanValidation( data, params, statespace, lifted );
    [error, koopsim] = val_liftedSys_loaded(data, lifted);
    disp('Done.')
end

%% Convert to iddata and compare to validation data
if params.compareon
    disp('Converting to iddata format...');
    data4sysid = get_data4sysid( data , koopsim , params );
    disp('Done.')
end

%% Save the model file if the user agrees
saveModel = questdlg('Would you like to save this model?');
if strcmp( saveModel , 'Yes' )
    % save the dynamics
    save_model(lifted);
    
    % save the lifting function
    liftFunDest = ['liftingFunctions' , filesep , 'lift_' , params.systemName, '.m'];
    copyfile('stateLift.m' , liftFunDest);
    
    % save the other lifting function
    WliftFunDest = ['liftingFunctions' , filesep , 'Wlift_' , params.systemName, '.m'];
    copyfile('Wlift.m' , WliftFunDest);
end

