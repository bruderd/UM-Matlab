% main: Run the Koopman sysid process on a user defined nonlinear system
%
%   DNE = "do not edit this line"
clear 

%% USER EDIT SECTION

% define parameters
params = struct;  % DNE

params.n = 2;   % dimension of state space
params.maxDegree = 3;   % maximum degree of monomial basis
params = def_polyLift(params);  % creates the lifting function, polyLift

params.Ts = 0.01;   % sampling period
params.numTrials = 10;  % number of "experiments" to run from different initial conditions
params.x0ub = [1, 1]';  % upper bound of initial state condition
params.x0lb = [-1, -1]'; % lower bound of initial state condition
params.tf = 10; % simulation time for each "experiment"


% define dynamics
x = params.x;    % DNE

vf = [x(2);...
      (1 - x(1)^2) * x(2) - x(1)];
 
matlabFunction(vf, 'File', 'vf_real', 'Vars', {x}); % create plant dynamics func. DNE


%% Simulate and find Koopman operator from "measurements" (DNE)

[x,y] = get_snapshotPairs(params);
U = get_Koopman(x,y);

%% Calculate the infiniesimal generator as funtion of coeffients, and from data (DNE)
Ldata = get_Ldata(U, params);   % infinitesimal generator from data
vecLdata = Ldata(:);    % vectorized version of Ldata matrix

vecstackL = zeros(params.N^2, params.N*params.n);
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

% dynamics (gives symbolic expression in terms of state)
vf2 = w * params.polyBasis; 
matlabFunction(vf2, 'File', 'vf_sysid', 'Vars', {params.x});

%% Run simulatio of sysId'd system and compare results to real system (DNE)

tspan = [0, params.tf];
x0sim = (params.x0ub - params.x0lb)*rand + params.x0lb; % random initial state
[tsysid, xsysid] = ode45(@(t,x) vf_sysid(x), tspan, x0sim);
[treal, xreal] = ode45(@(t,x) vf_real(x), tspan, x0sim);

% plot the results
figure
subplot(2,1,1)
plot(treal,xreal)
title('Real system')
subplot(2,1,2)
plot(tsysid,xsysid)
title('Identified system')


