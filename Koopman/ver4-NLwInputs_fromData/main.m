% main: Run the Koopman sysid process on a data collected from a real sys
%
%   DNE = "do not edit this line"
clear all;

%% USER EDIT SECTION

% define parameters
params = struct;  % DNE

params.n = 4;   % dimension of state space (including state derivatives)
params.p = 1;   % dimension of input
params.naug = params.n + params.p; % dimension of augmented state (DNE)

% select maximum degrees for monomial bases (NOTE: m1 = 1 )
params.maxDegree = 1;   % maximum degree of vector field monomial basis
params.m1 = 1;  % maximum degree of observables to be mappee through Lkj (DNE)
params = def_polyLift(params);  % creates the lifting function, polyLift

% choose whether or not to take numerical derivatives of states ('on' or 'off')
numericalDerivs = 'on';

params.Ts = 0.2;   % sampling period


%% Read in user data from csv

% Option 2nd argument for whether numerical derivative should be taken
data = get_data(params, numericalDerivs);


%% Simulate and find Koopman operator from "measurements" (DNE)

[x,y] = get_snapshotPairs(data, params);
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
vf2 = w * params.polyBasis; 
matlabFunction(vf2, 'File', 'vf_sysid', 'Vars', {params.x, params.u});

%% Run simulatio of sysId'd system and compare results to real system (DNE)

simlen = 500;   % simulate for first 500 time points of data
tspan = [0, data.t(simlen)];    
x0sim = data.x(1,:)'; % same initial state as data initial state
[tsysid, xsysid] = ode45(@(t,x) vf_sysid(x, get_simInput(t, data, params)), tspan, x0sim);

% plot the results
figure
subplot(2,1,1)
plot(data.t(1:simlen), data.x(1:simlen,:))
title('Real system')
subplot(2,1,2)
plot(tsysid, xsysid)
title('Identified system')


