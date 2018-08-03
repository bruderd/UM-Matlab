function data = run_experiment(dataFileName)
%RUN_EXPERIMENT Summary of this function goes here
%   Detailed explanation goes here
clearvars -except dataFileName;

%% Define system parameters
params = struct;

% initial conditions
params.phi1                = 0;
params.dtphi1              = 0;
params.phi2                = 0;
params.dtphi2              = 0;

% physical parameters
params.g                   = 9.81; 
params.m1                  = 1; 
params.m2                  = 1; 
params.l1                  = 1; 
params.l2                  = 1;

% "measurement" parameters
params.Ts                  = 1/30;
params.mean                = 0;     % mean of noise 
params.sigma               = 0;     % standard dev of noise

% animation parameters
params.duration            = 3600;   % in seconds
params.fps                 = 30;
params.movie               = true;

%% Simulate system

x0 = [params.phi1; params.phi2; params.dtphi1; params.dtphi2];
sol = ode45(@(t, x) double_pendulum_ODE(t, x, get_input(t,x,params), params),[0 params.duration], x0);
tout = sol.x;
usol = get_input(tout);       % get inputs corresponding to output points

%% Get state "measurements"

% get the state at sampled times
t = (0 : params.Ts : params.duration)';
y = deval(sol, t)';
u = interp1(tout, usol, t);

% define observed states (this may be different than simulation states, 
%   e.g. could track only angular positions, or the xy position of end eff)
xobs = y(:,1:4);    % full state

% inject noise to simulate measurement noise
noise = params.sigma .* randn(size(xobs)) + params.mean;
x = xobs + noise;

% define output
data = struct;
data.t = t;
data.x = x;
data.u = u;

%% Save experimental data as .csv and .mat file
if exist('dataFileName', 'var')
    csvwrite([dataFileName, '.csv'], [data.t, data.x, data.u]);
    save([dataFileName, '.mat'], 't', 'x', 'u');
end

end



%% Define the system input as a function of time (and state if desired)

function u = get_input(t,x,params)
%   will want to parametrize in terms of some params later...

% u = 10*sin( (1/(2*pi)) * t) .* sin( 3*t - 1.5*cos(t) );
u = 10*sin(0.1*t) + cos(t);
% u = 20 * rand  - 10;


end

