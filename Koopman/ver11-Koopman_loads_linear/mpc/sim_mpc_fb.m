% sim_mpc_fb: simulate an mpc trial on hardware using model for prediction
%   The "fb" suffix denotes that this version of calls get_MPCinput_fb
%   NOTE: Only compatible with models that have 0-1 delays.

clear all;  % this guarantees that the persistent variables have been cleared from getMPCinput function

%% initializations
% Load model and other system parameters into workspace
[modelFile, modelPath] = uigetfile(['..' , filesep , 'models'] , 'Choose model file...');    % load model parameters
[mpcParams, mpcParamsPath] = uigetfile(['paramFiles'] , 'Choose mpc parameter file...'); 
[trajectory, trajectoryPath] = uigetfile(['trajectory' , filesep , 'files'], 'Choose reference trajectory file...');

model = load([modelPath , filesep , modelFile]);    % load model file
load([mpcParamsPath , filesep , mpcParams]);   % load mpc params
load([trajectoryPath , filesep , trajectory]);    % load reference trajectory

Ts = model.params.Ts;    % sampling time of model

% resample reference trajectory
tr = 0 : Ts : ref.T;
yr_unsc = interp1(ref.t , ref.y , tr);

% scale reference trajectory for mpc
yr = yr_unsc * diag(model.params.xScaleFactor);

% Define mpc matrices
mpc = set_mpc( model , mpcParams );

% initialize waypoint index
now = 1;

% initialize mpc data struct
mpcData = struct;
mpcData.T = [];
mpcData.U = [];
mpcData.Y = [];
mpcData.K = [];
mpcData.Yr = [];
mpcData.Psi = [];
mpcData.load = [];

%% run simulation

% get handle for the lifting function
cd( ['..' , filesep , 'liftingFunctions'] );
stateLift = str2func( [ 'lift_' , model.params.systemName ] );
cd(['..' , filesep , 'mpc']);

mpcData.Y = zeros( 1 , model.params.n ); %[0,0];
mpcData.U = zeros( 2 , model.params.p ); % [0 0 0 ; 0 0 0];
mpcData.K = 0;
mpcData.T = 0;
W = zeros(model.params.N * mpc.params.Nw , 1 + model.params.nw);

i = 1;
while now < length(yr)  % stop when all trajectory points have been targeted
% for i = 1 : length(yr)
    
    % current time
    t = (i-1) * model.params.Ts;
    
    % current step
    k = i;
    
    % Calculates current output (i.e. laser position)
    yk = mpcData.Y(i,:)';     % current measured output
    
    % scale down laser data for calculating input
    ysck = diag(model.params.xScaleFactor , 0) * yk;
    
    % Update the current waypoint
    now = now + 1; wp = yr( now , :);   % increment always
%     [ wp , now ] = get_waypoint( ysck , yr , now , 0.025 );   % increment only if within epsilon ball
    
    % scale down the previous state and input
    ysckm1 = diag( model.params.xScaleFactor ) * mpcData.Y(end,:)' ; 
    ucskm1 = diag( model.params.uScaleFactor ) * mpcData.U(end-1:end,:)' ;
    
    % estimate the load
    [ load , W ] = get_loadEstimate( yk , W , mpcData , mpc , model);    % accepts unscaled measurements
    
    % caluclate input
    usck = get_MPCinput_fb( k , ysck , ysckm1, ucskm1, load , wp, mpc , model );
    uk = diag(model.params.uScaleFactor .^ (-1) , 0) * usck;  % scale input back up
    
    %% simulate the system from k to k+1
    if i <= model.params.nd     % all stuff before trial starts is zero
        xk = ysck;
        xdk = zeros( model.params.n * model.params.nd , 1 );
        udk = zeros( model.params.p * model.params.p , 1 );
    else
        xk = ysck;
        xdk = reshape( flipud( mpcData.Y( i-model.params.nd : i-1 , : ) * diag(model.params.xScaleFactor) )' , [model.params.n * model.params.nd , 1] );
        udk = reshape( flipud( mpcData.U( i-model.params.nd : i-1 , : ) * diag(model.params.uScaleFactor) )' , [model.params.p * model.params.nd , 1] );
    end
    zetak = [xk; xdk; udk];
    psik = stateLift(zetak,0.5);  % NEED TO INPUT THE REAL LOAD HERE!
    psikp1 = model.Asim * psik + model.Bsim * usck;
    xkp1 = diag(model.params.xScaleFactor.^(-1)) * ( model.C * psikp1 );    % isolate state and scale back up
    
    % update values of output parameters
    mpcData.T = [ mpcData.T ; t ];
    mpcData.Y = [ mpcData.Y ; xkp1' ];
    mpcData.U = [ mpcData.U ; uk' ];
    mpcData.K = [ mpcData.K ; k ];
    mpcData.Yr = [ mpcData.Yr ; yr_unsc(now , :) ];    % this is the unscaled trajectory
    mpcData.Psi = [ mpcData.Psi ; psik' ];   % the lifted state at each timestep
    mpcData.load = [ mpcData.load ; load' ];
    
    i = i + 1;      % increment step counter
end

%% plot results

figure; plot(mpcData.T, mpcData.U)
figure; plot(mpcData.Y(:,1), mpcData.Y(:,2))
hold on
plot(mpcData.Yr(:,1),mpcData.Yr(:,2),'*')
hold off



