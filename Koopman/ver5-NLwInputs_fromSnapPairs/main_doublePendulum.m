function main_doublePendulum(snapshotPairs)
% main_doublePendulum: Run the Koopman sysid process on data collected
% from doublePendulum simulation
%   snapshotPairs :     struct contianing fields "x" and "y" which are the
%                       snapshot pairs for the Koopman operator calculation
%
%   DNE = "do not edit this line"

% clearvars -except dataFileName;

%% USER EDIT SECTION: set parameter values

params = struct;  % DNE
progress = waitbar(0,'Initializing parameters...');

% Koopman Sysid parameters
params.n = 4;   % dimension of state space (including state derivatives)
params.p = 1;   % dimension of input
params.naug = params.n + params.p; % dimension of augmented state (DNE)

% select maximum degrees for monomial bases (NOTE: m1 = 1)
params.maxDegree = 1;   % maximum degree of vector field monomial basis
params.m1 = 1;  % maximum degree of observables to be mapped through Lkj (DNE)

% define lifting function and basis
% params = def_polyLift(params);  % creates polynomial lifting function, polyLift
params.interval = [-50, 50];     % interval (in all states) upon which dynamics are defined
params = def_sinLift(params);   % creates sinusoidal lifting function, sinLift

% choose whether or not to take numerical derivatives of states ('on' or 'off')
numericalDerivs = 'off';

params.Ts = 1/30;   % sampling period

% animation parameters
params.duration            = 5;
params.fps                 = 30;
params.movie               = true;

% double pendulum parameters
params.phi1                = pi/2;
params.dtphi1              = 0;
params.phi2                = pi/2;
params.dtphi2              = 0;
params.g                   = 9.81; 
params.m1                  = 1; 
params.m2                  = 1; 
params.l1                  = 1; 
params.l2                  = 1;


%% Simulate and find Koopman operator from "measurements" (DNE)
waitbar(.33,progress,'Calculating Koopman Operator...');

[x,y] = deal(snapshotPairs.x, snapshotPairs.y);
U = get_Koopman(x,y, params);

%% Calculate the infiniesimal generator as funtion of coeffients, and from data (DNE)
Ldata = get_Ldata(U, params);   % infinitesimal generator from data
Ldata_hat = Ldata(:, 1:params.N1);  % N x N1 version of Ldata (projection onto the polyBasis)
vecLdata = Ldata_hat(:);    % vectorized version of Ldata matrix

vecstackL = zeros(params.N*params.N1, params.N*params.n);
for k = 1:params.N
    for j = 1:params.n
        Lkj = get_Lkj(k,j,params);
        
        % convert all the Lkj's into vectors and stack them horizontally
        vecLkj = Lkj(:);
        vecstackL(:, (k-1)*params.n + j) = vecLkj;
    end
end

%% solve for the coefficients, i.e. Eq. (18) from Mauroy and Gonclaves (DNE)

W = pinv(vecstackL) * vecLdata;

% matrix of coefficents of monomials
w = reshape(W, [params.n, params.N]);

% dynamics (gives symbolic expression in terms of state and input)
% vf2 = w * params.polyBasis; 
vf2 = w * params.sinBasis; 
matlabFunction(vf2, 'File', 'vf_sysid', 'Vars', {params.x, params.u});

%% Run simulatio of sysId'd system and compare results to real system (DNE)
waitbar(.67,progress,'Simulating dynamics...');

tspan = [0, params.duration];    
x0sim = [params.phi1; params.phi2; params.dtphi1; params.dtphi2]; % same initial state as data initial state
sol_sysid = ode45(@(t,x) vf_sysid(x, get_input(t, x, params)), tspan, x0sim);
sol_real = ode45(@(t,x) vf_doublePendulum(x, get_input(t, x, params), params), tspan, x0sim);

[tsysid, xsysid] = deal(sol_sysid.x, sol_sysid.y);
[treal, xreal] = deal(sol_real.x, sol_real.y);

% plot the results
figure
subplot(2,1,1)
plot(treal, xreal)
title('Real system')
subplot(2,1,2)
plot(tsysid, xsysid)
title('Identified system')

% animate the results
animate_doublePendulum(sol_real, sol_sysid, params);

waitbar(1,progress,'Done.');
end

%% USER EDIT: Define input for the simulation

function u = get_input(t,x,params)
%   will want to parametrize in terms of some params later...

% u = 4*sin( (1/(2*pi)) * t) .* sin( 3*t - 1.5*cos(t) );
% u = 3*(1 - exp(-0.1*t));
u = 0;
% u = 10*sin(0.1*t) + cos(t);

end

